import itertools
from typing import Any, List, Tuple
import sys

import click
import pipe


def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.islice(zip(a, b), None, None, 1)


def is_new_junction(bed12_line, parsed_junctions):
    fields = bed12_line.split('\t')
    sizes = [int(i) for i in fields[10].rstrip(',').split(',')]
    starts = [int(i) for i in fields[11].rstrip(',').split(',')]
    coord1 = int(fields[1]) + starts[0] + sizes[0] + 1
    coord2 = int(fields[1]) + starts[1]
    # site1 = (fields[0], str(int(fields[1]) + starts[0] + sizes[0] + 1))
    # site2 = (fields[0], str(int(fields[1]) + starts[1] - 1))
    
    return not (((fields[0], str(coord1)) in parsed_junctions) or (((fields[0], str(coord2)) in parsed_junctions))) and (coord2 - coord1 > 0)


def junction_to_bed_entry(bed12_line):
    fields = bed12_line.split('\t')
    sizes = [int(i) for i in fields[10].rstrip(',').split(',')]
    starts = [int(i) for i in fields[11].rstrip(',').split(',')]
    site1 = str(int(fields[1]) + starts[0] + sizes[0])
    site2 = str(int(fields[1]) + starts[1])
    fields[1], fields[2] = site1, site2
    return '\t'.join(fields[:9]) 


def parse_SJDB(SJDB_file):
    with open(SJDB_file) as f:
        return set(
            f.readlines()
            | pipe.where(lambda l: (len(l.split()) >= 3))
            | pipe.map(lambda l: l.split())
            | pipe.map(lambda columns: [(columns[0], columns[1]), (columns[0], columns[2])])
            | pipe.chain
        )


@click.command()
@click.option("--sjdb", help="BED3 file with junctions")       
def main(sjdb):
    parsed_junctions = parse_SJDB(sjdb)
    filtered_entries = (
        sys.stdin 
        | pipe.where(lambda l: is_new_junction(l, parsed_junctions)) 
        | pipe.map(junction_to_bed_entry)
    )
    for junction in filtered_entries:
        print(junction.strip())


if __name__ == "__main__":
    main()