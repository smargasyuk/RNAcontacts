import sys

import click
from colour import Color

def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)

arm_len = 20
uplim = 10
colors = list(Color("blue").range_to(Color("red"), uplim))

def color_map(score):
    return ",".join(map(lambda x: str(int(x * 255)), colors[clamp(score, 0, uplim - 1)].rgb))

def generate_view(line):
    fields = line.strip().split("\t")
    fields = fields[:6]
    fields[3] = f"{fields[3]};count={fields[4]}"
    fields += [fields[1], fields[2], color_map(int(fields[4])), "2"]
    fields.append(f"{arm_len},{arm_len}")
    fields.append(f"0,{int(fields[2]) - int(fields[1]) + arm_len}")
    fields[1], fields[2] = str(int(fields[1]) - arm_len), str(int(fields[2]) + arm_len)
    return "\t".join(fields)


@click.command()
def main():
    for junction in sys.stdin:
        print(generate_view(junction))


if __name__ == "__main__":
    main()