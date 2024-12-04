#!/usr/bin/perl 

# This script reads a list of junctions from $ARGV[0] and parses 
# a SAM from STDIN to extract and print neojunctions in STAR chimeric 
# format: chr1 pos1 str1 chr2 pos2 str2 read_id


$BAM_FREAD1 = 0x40;
$BAM_FREAD2 = 0x80;
$BAM_FREVERSE = 0x10;

open FILE, $ARGV[0];
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end) = split /\t/, $line;
    $intron{$chr}{$beg}{$end}++;
}
close FILE;

@STRAND = ("+","-");
while(<STDIN>) {
    ($id, $flag, $ref, $pos, $qual, $cigar) = split /\t/;
    next if($flag & 0x100); # skip secondary alignments

    $rev = ($flag & $BAM_FREVERSE) ? 1 : 0; # reverse complemented yes no
    $str = $STRAND[$rev]; # report reverse flag from BAM, do not flip it

    while($cigar=~/(\d+)(\w)/g) {
        $increment = $1;
        $operation = $2;
        if($operation eq 'M') {  
                $pos += $increment;
        }
        if($operation eq 'D') {  
                $pos += $increment;
        }
        if($operation eq 'N') {
                $beg = $pos;
                $end = $pos + $increment - 1;
                print join("\t", $ref, $beg, $str, $ref, $end, $str, $id), "\n" unless($intron{$ref}{$beg}{$end});
                $pos += $increment; 
        }
    }
}
