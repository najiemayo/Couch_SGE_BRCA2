# ~/python/3.7.5/bin/python


import argparse

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment


def main():
    args = parse_args()
    run(args.fastq, args.sequence, args.complement)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-f', dest='fastq',
                        help='test name of the panel')
    parser.add_argument('-c', dest='complement', action='store_true',
                        help='to map the reverse complement of the reads')
    parser.add_argument('-r', dest='reference',
                        default='~/GRCh38_full_analysis_set_plus_decoy_hla.fa',
                        help='(optional/required if non-human) reference fasta of organism')
    parser.add_argument('-b', dest='bam',
                        help='bam file as input')
    parser.add_argument('-s', dest='sequence', required=True, help='Target sequence to use')
    args = parser.parse_args()
    return args


def run(fastq, target_seq, complement):
    #    if fastq:
    if not complement:
        parse_fastq(fastq, target_seq)
    else:
        parse_fastq(fastq, target_seq, complement)


#    elif bam:


def parse_fastq(filename, seq, complement=False):
    with open(filename, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            seq1 = record.seq
            seq2 = Seq(seq)
            if not complement:
                alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
                for alignment in alignments:
                    print(format_alignment(*alignment))
            else:
                new_seq2 = seq2.reverse_complement()
                alignments = pairwise2.align.globalms(seq1, new_seq2, 2, -1, -0.5, -0.1)
                for alignment in alignments:
                    print(format_alignment(*alignment))


# def parse_bam

if __name__ == "__main__":
    main()
