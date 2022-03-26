import itertools
from typing import Any, List, Tuple

import click
import pysam
import pipe

Segment = Tuple[str, int, int]
CigarSegment = Tuple[int, int]
Read = Any


def get_next_segment(current_segment, next_cigar_segment) -> Segment:
    # coordinates of the next CIGAR segment in genome coordinates
    genome_shift = (
        next_cigar_segment[1] if next_cigar_segment[0] in [0, 2, 3] else 0
    )  # shift if code is M, D, N
    return (
        current_segment[0],
        current_segment[2],
        current_segment[2] + genome_shift,
        next_cigar_segment[0],
    )


def junctions_from_single_read(r1) -> List[Segment]:
    read_start = (
        r1.header.get_reference_name(r1.reference_id),
        r1.reference_start,
        r1.reference_start,
        None,
    )
    segments = (
        itertools.accumulate(r1.cigartuples, func=get_next_segment, initial=read_start)
        | pipe.where(lambda x: x[3] == 3)  # code is N
        | pipe.map(lambda seg: (seg[0], seg[1] + 1, seg[2], seg[3]))
    )  # return from pysam to sam 1-based coordinates
    return list(segments)


def extract_juctions_from_read(r1):
    r1_junctions = junctions_from_single_read(r1)
    strand = "-" if r1.is_reverse else "+"
    return [
        (chrom, start, end)
        for chrom, start, end, _ in r1_junctions
    ]


def parse_SJDB(SJDB_file):
    with open(SJDB_file) as f:
        return set(
            f.readlines()
            | pipe.where(lambda l: (len(l.split()) >= 3))
            | pipe.map(lambda l: l.split())
            | pipe.map(lambda columns: (columns[0], int(columns[1]), int(columns[2])))
        )


def is_known_junction(junction, SJDB_set):
    return junction in SJDB_set


def output_junction(junction):
    print("\t".join(str(c) for c in junction))

def filter_length(read, max_length):
    if max_length == -1:
        return True
    else:
        return read.reference_length <= max_length

@click.command()
@click.argument("sam")
@click.argument("sjdb", nargs=-1)
@click.option("--max-length", default=-1, type=int)
def main(sjdb, sam, max_length: int):
    SJDB_set = set.union(*[parse_SJDB(s) for s in sjdb])
    # print(SJDB_set)
    samfile = pysam.AlignmentFile(sam)
    outfile = pysam.AlignmentFile("-", "w", template=samfile)
    reads_w_new_junctions = (
        samfile
        | pipe.where(lambda read: 'N' in read.cigarstring)
        | pipe.where(lambda read: not all(is_known_junction(j, SJDB_set) for j in extract_juctions_from_read(read)))
        | pipe.where(lambda read: filter_length(read, max_length))
    )
    for r in reads_w_new_junctions:
        outfile.write(r)
    # print(list(SJDB_set)[:100])
    # for s in samfile:
    #     junctions = extract_juctions_from_read(s)
    #     # print(junctions)
    #     # has_gap = ('N' in s.cigarstring)
    #     # junctions_are_known = list(is_known_junction(j, SJDB_set) for j in junctions)
    #     # print(('N' in s.cigarstring), list(is_known_junction(j, SJDB_set) for j in junctions),
    #     #  not all(is_known_junction(j, SJDB_set) for j in junctions) and ('N' in s.cigarstring))
    #     if not all(is_known_junction(j, SJDB_set) for j in junctions) and ('N' in s.cigarstring):
    #         outfile.write(s)

if __name__ == "__main__":
    main()
