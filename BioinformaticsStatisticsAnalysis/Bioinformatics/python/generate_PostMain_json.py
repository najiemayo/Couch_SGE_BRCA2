#!~/python/3.7.5/bin/python

import argparse
import json


def main():
    args = parse_args()
    run(args.sample,
        args.vcf,
        args.out_dir,
        args.reference,
        args.json_temp)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', dest='sample', required=True,
                        help='name of the samples')
    parser.add_argument('-v', dest='vcf', required=True,
                        help='path to vcf')
    parser.add_argument('-o', dest='out_dir', required=True,
                        help='path to output dir')
    parser.add_argument('-r', dest='reference',
                        default='~/GRCh38_full_analysis_set_plus_decoy_hla.fa',
                        help='(optional/required if non-human) reference fasta of organism')
    parser.add_argument('-j', dest='json_temp', default='wdl/PostMain.json')
    args = parser.parse_args()
    return args


def run(sample, vcf, out_dir, reference, json_temp):
    parse_json(sample, vcf, out_dir, reference, json_temp)


def parse_json(sample, vcf, out_dir, reference, json_dump):
    f = open(json_dump, 'r')
    data = json.load(f)
    modify_dict(data,sample, vcf, out_dir, reference)


def modify_dict(data, sample, vcf, out_dir, reference):
    new_dict = data
    sample = sample
    new_dict['SeqToAnno.INPUT_VCF'] = vcf
    new_dict['SeqToAnno.SAMPLENAME'] = sample
    json_outname = out_dir + "/" + sample + ".PostMain.input.json"
    with open(json_outname, "w") as outfile:
        json.dump(new_dict, outfile, indent=4)


if __name__ == "__main__":
    main()
