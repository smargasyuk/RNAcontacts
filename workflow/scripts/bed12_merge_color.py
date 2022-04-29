import click
import distinctipy


@click.command()
@click.argument("files", nargs=-1)
def main(files):

    colors = distinctipy.get_colors(len(files), pastel_factor=0.7)
    color_strings = [
        ",".join([str(int(i * 255)) for i in rgb_frac_tuple])
        for rgb_frac_tuple in colors
    ]

    for file_idx, bed in enumerate(files):
        with open(bed) as f:
            for line in f.readlines():
                fields = line.strip().split("\t")
                fields[8] = color_strings[file_idx]
                print("\t".join(fields))


if __name__ == "__main__":
    main()
