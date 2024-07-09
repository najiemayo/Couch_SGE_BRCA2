#!~/python/3.7.5/bin/python

import argparse
import json
import os
from collections import defaultdict


def main():
    args = parse_args()
    run(args.fastq_dir,
        args.out_dir, args.cutadapt_flag, args.seqprep_flag,
        args.reference, args.json_temp)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-f', dest='fastq_dir', required=True,
                        help='test name of the panel')
    parser.add_argument('-o', dest='out_dir', required=True,
                        help='path to output dir')
    parser.add_argument('-c', dest='cutadapt_flag', required=True,
                        help='string for cutadapt flag')
    parser.add_argument('-s', dest='seqprep_flag', required=True,
                        help='string for seqprep flag')
    parser.add_argument('-r', dest='reference',
                        default='~/GRCh38_full_analysis_set_plus_decoy_hla.fa',
                        help='(optional/required if non-human) reference fasta of organism')
    parser.add_argument('-j', dest='json_temp', default='wdl/input.json')
    args = parser.parse_args()
    return args


def run(fastq_dir, out_dir, cutadapt_flag, seqprep_flag, reference, json_dump):
    parse_json(fastq_dir,
               out_dir, cutadapt_flag, seqprep_flag,json_dump)


def parse_json(fastq_dir, out_dir, cutadapt_flag, seqprep_flag, json_temp):
    f = open(json_temp, 'r')
    data = json.load(f)
    modify_dict(data, fastq_dir, out_dir, cutadapt_flag, seqprep_flag)


def modify_dict(data, fastq_dir, out_dir, cutadapt_flag, seqprep_flag):
    sample_fastq_dict = get_sample_name(fastq_dir)
    print(sample_fastq_dict)
    for sample, fastq in sample_fastq_dict.items():
        outdir = out_dir + '/' + sample

        # Read 1 is not always first
        if 'L1_R1' in fastq[0]:
            R1_fastq = fastq[0]
            R2_fastq = fastq[1]
        else:
            R1_fastq = fastq[1]
            R2_fastq = fastq[0]
        new_dict = data
        new_dict['FastqToSam.FASTQ1'] = R1_fastq
        new_dict['FastqToSam.FASTQ2'] = R2_fastq
        new_dict['FastqToSam.sortsam.OUTDIR'] = outdir
        new_dict['FastqToSam.OUTDIR'] = outdir
        new_dict['FastqToSam.CUTADAPT_OPTIONS'] = cutadapt_flag
        new_dict['FastqToSam.SEQPREP_OPTIONS'] = seqprep_flag
        json_outname = out_dir + "/" + sample + ".input.json"
        with open(json_outname, "w") as outfile:
            json.dump(new_dict, outfile, indent=4)


def get_sample_name(fastq_dir):
    fastq_dict = defaultdict(list)
    for fname in os.listdir(fastq_dir):
        if fname.endswith("fastq.gz"):
            path = os.path.join(fastq_dir, fname)
            sample_name = os.path.basename(path).split(".")[0]
            fastq_dict[sample_name].append(path)
    return fastq_dict


if __name__ == "__main__":
    main()
