import itertools
from typing import Any, List, Tuple
import sys

import click
import pipe
from intervaltree import Interval, IntervalTree


def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.islice(zip(a, b), None, None, 1)


def parse_SJDB(SJDB_file, radius):
    tree_dict = {}
    with open(SJDB_file) as f:
        for line in f.readlines():
            fields = line.split()
            if len(fields) < 3:
                continue
            l_chrom, l_start, l_end = fields[0], int(fields[1]), int(fields[2])
            interval1, interval2  = Interval(l_start - radius, l_start + radius), Interval(l_end - radius, l_end + radius)
            if l_chrom in tree_dict:
                tree_dict[l_chrom].add(interval1)
                tree_dict[l_chrom].add(interval2)
            else:
                tree_dict[l_chrom] = IntervalTree([interval1, interval2])
    return tree_dict

def is_new_junction(l, parsed_junctions):
    fields = l.split()
    l_chrom, l_start, l_end = fields[0], int(fields[1]), int(fields[2])
    if l_chrom not in parsed_junctions:
        return False
    return (not parsed_junctions[l_chrom].overlaps(l_start)) and (not parsed_junctions[l_chrom].overlaps(l_end))

@click.command()
@click.option("--sjdb", help="BED3 file with junctions")     
@click.option("--radius", help="radius from known junction", type=int, default=25)  
def main(sjdb, radius):
    parsed_junctions = parse_SJDB(sjdb, radius)
    filtered_entries = (
        sys.stdin 
        | pipe.where(lambda l: is_new_junction(l, parsed_junctions)) 
    )
    for junction in filtered_entries:
        print(junction.strip())


if __name__ == "__main__":
    main()