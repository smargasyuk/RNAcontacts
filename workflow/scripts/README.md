Scripts for RNAcontacts pipeline

+ `neo.pl` : extracts junctions from stdin in SAM format for mate `ARGV[1]` and prints junctions not present in file `ARGV[0]`;

+ `bed12_merge_color.py`: merges bed12 files and assigns different colors to entries from different files;

+ `bedj_filter_junctions.py`: reads junctions in bed3 format from stdin and prints junctions with both ligation sites outside `radius` from ligation sites from `sjdb` file with junctions;

+ `bedj_to_bed12.py`: reads junctions in bed3 format from stdin and prints junctions in bed12 format, where 'intron' is an actual junction and 'exons' are 20nt flanking regions;

+ `bed_enumerate.py`: reads junctions in bed format from stdin and prints them with names containing line idx and read support






