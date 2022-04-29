import sys


def main():
    for idx, line in enumerate(sys.stdin):
        fields = line.strip().split("\t")
        fields[3] = f"id={idx},count={fields[4]}"
        print("\t".join(fields))


if __name__ == "__main__":
    main()
