# ~/python3_env/bin/python

"""
This scripts takes screenshots of 
the variant in a bam file

Author: Rohan David Gnanaolivu
"""

import argparse
import logging
import os
import subprocess
import sys

# IGVTtools works with Java 11 or greater
JAVA = "~/java/jdk-11.0.9+11/bin/java"


def main():
    logging.basicConfig(stream=sys.stderr, level=logging.ERROR,
                        format='%(name)s (%(levelname)s): %(message)s')
    args = parse_input_args()
    read_vcf(args.BAMLIST,
             args.IGVTOOLS, args.out_dir,
             args.bed)


def parse_input_args():
    parser = argparse.ArgumentParser(description='Create IGV Snapshots for variants of interest from vcf files',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=""" 
            NOte: this only works if you have X11 forwarding on!
            uses a vcf file as input
        """)
    parser.add_argument('-d', help='output directory', dest='out_dir', type=str, default=os.getcwd())
    parser.add_argument('-j', help='igvtools path', default='~/igv/igv-2.8.9/lib/', dest='IGVTOOLS')
    parser.add_argument('-f', type=str, help='The file to parse', required=True, dest='bed')
    parser.add_argument('-b', type=str, help='A file containing a list of BAM files', required=True, dest='BAMLIST')
    parser.add_argument('-l', '--verbose', dest='verbose_count', action='count', default=0,
                        help='increases log verbosity for each occurence.')
    args = parser.parse_args()
    return args


def read_bam_list(bamlist):
    sample_tuple_list = []
    with open(bamlist, 'r') as f:
        for line in f:
            new_line = line.strip()
            sample_name = os.path.basename(new_line).split(".")[0]
            sample_tuple_list.append((sample_name, line))
    return sample_tuple_list


def create_script(out_dir, bam, pos,
                  sample_id, variant, bed):
    """
    This function create a igv batch script to run with
    igvtools
    to modify how to view the image modify this batch code
    """
    lines = ['new',
             'genome ~/GRCh38_full_analysis_set_plus_decoy_hla.fa',
             'load ' + bam, 'snapshotDirectory ' + out_dir, 'goto ' + pos, 'sort base', 'squish', 'load ' + bed,
             'snapshot ' + sample_id + '_' + variant + '.png', 'exit']
    #    lines.append('genome hg38')
    line = '\n'.join(lines)
    return line


def read_vcf(bamlist, igvtools, out_dir, bed):
    """
    This function parses the vcf and adds 10bp to the variant
    to view on IGV
    also i add in 4g memory which works for IDT bams
    """
    sample_tuple = read_bam_list(bamlist)
    with open(bed, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                line = line.strip().split('\t')
                start_pos = int(line[1]) - 50
                end_pos = int(line[2]) + 50
                pos = str(line[0]) + ':' + str(start_pos) + '-' + str(end_pos)
                region = str(line[0]) + '_' + line[1] + '_' + line[2]
                for some_tuple in sample_tuple:
                    bam_file = some_tuple[1]
                    cmd = create_script(out_dir, bam_file, pos, some_tuple[0], region, bed)
                    with open('tmp.batch', 'w') as f:
                        f.write(cmd)
                        module_path = '--module-path=' + igvtools
                    subprocess.call([JAVA, module_path, '-Xmx4g', '@~/igv/igv-2.8.9/igv.args',
                                     '--module=org.igv/org.broad.igv.ui.Main', '-b', 'tmp.batch'])


if __name__ == "__main__":
    main()
