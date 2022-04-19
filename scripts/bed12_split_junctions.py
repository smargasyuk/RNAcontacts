import itertools
from typing import Any, List, Tuple
import sys


def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.islice(zip(a, b), None, None, 1)


def get_new_lines(line: str):
    orig_fields = line.split('\t')
    orig_fields[10] = [int(i) for i in orig_fields[10].rstrip(',').split(',')]
    orig_fields[11] = [int(i) for i in orig_fields[11].rstrip(',').split(',')]
    for block0, block1 in pairwise(zip(orig_fields[10], orig_fields[11])):
        new_line = orig_fields.copy()
        old_start = int(new_line[1])
        new_line[2] = old_start + block1[0] + block1[1]
        new_line[1] = old_start + block0[1]
        new_line[6] = new_line[1]
        new_line[7] = new_line[2]
        new_line[9] = 2
        new_line[10] = ','.join(map(str, [block0[0], block1[0]]))
        new_line[11] = ','.join(map(str, [0, block1[1] - block0[1]]))
        yield '\t'.join(map(str, new_line))
        

def main():
    new_entries = map(get_new_lines, sys.stdin)
    for read in new_entries:
        for junction in read:
            print(junction)


if __name__ == "__main__":
    main()