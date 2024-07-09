#!~/python/3.7.5/bin/python

import argparse
import json


def main():
    args = parse_args()
    run(args.chrom,
        args.out_dir,
        args.position,
        args.sequence,
        args.reference, args.json_temp)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', dest='chrom', required=True,
                        help='chromosome with chr prefix')
    parser.add_argument('-o', dest='out_dir', required=True,
                        help='path to output dir')
    parser.add_argument('-p', dest='position', required=True,
                        help='start position of the sequence')
    parser.add_argument('-s', dest='sequence', required=True,
                        help='sequence to create all possible mutations')
    parser.add_argument('-r', dest='reference',
                        default='~/GRCh38_full_analysis_set_plus_decoy_hla.fa',
                        help='(optional/required if non-human) reference fasta of organism')
    parser.add_argument('-j', dest='json_temp', default='wdl/input.json')
    args = parser.parse_args()
    return args


def run(chrom, out_dir, position, sequence, reference, json_dump):
    parse_json(chrom, out_dir, position, sequence, reference, json_dump)


def parse_json(chrom, out_dir, position, sequence, reference, json_dump):
    f = open(json_dump, 'r')
    data = json.load(f)
    modify_dict(data, chrom, out_dir, position, sequence, reference)


def modify_dict(data, chrom, out_dir, position, sequence, reference):
    new_dict = data
    sample = "SSM_Allpossible"
    new_dict['SeqToAnno.CHR'] = chrom
    new_dict['SeqToAnno.START'] = position
    new_dict['SeqToAnno.SEQ'] = sequence
    new_dict['SeqToAnno.SAMPLENAME'] = sample
    json_outname = out_dir + "/" + sample + ".seqAnno.input.json"
    with open(json_outname, "w") as outfile:
        json.dump(new_dict, outfile, indent=4)


if __name__ == "__main__":
    main()
