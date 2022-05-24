import sys


arm_len = 20
color = "255,255,255"


def generate_view(line):
    fields = line.strip().split("\t")
    fields = fields[:6]
    fields[3] = f"{fields[3]};count={fields[4]}"
    fields += [fields[1], fields[2], color, "2"]
    fields.append(f"{arm_len},{arm_len}")
    fields.append(f"0,{int(fields[2]) - int(fields[1]) + arm_len}")
    fields[1], fields[2] = str(int(fields[1]) - arm_len), str(int(fields[2]) + arm_len)
    return "\t".join(fields)


def main():
    for junction in sys.stdin:
        print(generate_view(junction))


if __name__ == "__main__":
    main()
