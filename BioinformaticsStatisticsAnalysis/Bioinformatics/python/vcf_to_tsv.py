import argparse
import logging
import sys

from cyvcf2 import VCF


# noinspection DuplicatedCode
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--vcf",
                        dest='vcf_file',
                        required=True,
                        help="VCF file to extract columns from")

    parser.add_argument("-f", "--fields",
                        dest='fields_to_extract',
                        nargs='+',
                        default={'RefAA', 'AltAA', 'AApos', 'RefCodon', 'AltCodon', 'AltCodon', 'EventCount',
                                 'EventType', 'SpliceAI_ALLELE', 'SpliceAI_GENE', 'SpliceAI_DS_AG', 'SpliceAI_DS_AL',
                                 'SpliceAI_DS_DG', 'SpliceAI_DS_DL', 'SpliceAI_DP_AG', 'SpliceAI_DP_AL',
                                 'SpliceAI_DP_DG', 'SpliceAI_DP_DL'},
                        help="INFO fields to extract from VCF file")

    parser.add_argument("-o", "--out",
                        dest='tsv_out',
                        default=None,
                        help="TSV file to create")

    parser.add_argument("-V", "--verbose",
                        dest="logLevel",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default="INFO",
                        help="Set the logging level")

    args = parser.parse_args()
    logging.basicConfig(stream=sys.stderr, level=args.logLevel,
                        format='%(name)s (%(levelname)s): %(message)s')

    logger = logging.getLogger(__name__)
    logger.setLevel(args.logLevel)
    for x in args._get_kwargs():
        logger.info(x)

    return args, logger


if __name__ == "__main__":
    args, logger = parse_args()
    args = args.__dict__

    # Create an empty dictionary
    info_dict = dict(zip(sorted(args['fields_to_extract']), '.' * args['fields_to_extract'].__len__()))
    vcf = VCF(args['vcf_file'])

    # Write the header
    if args['tsv_out'] is None:
        o = open(args['vcf_file'].replace('.vcf', '.tsv'), 'w')
    else:
        o = open(args['tsv_out'], 'w')
    header = "#CHROM\tPOS\tREF\tALT\t"
    header += '\t'.join(sorted(info_dict.keys()))
    o.write(header + '\n')

    for variant in vcf:
        line = [variant.CHROM, variant.start, variant.REF, ','.join(variant.ALT)]
        for i in info_dict:
            try:
                line.append(variant.INFO[i])
            except KeyError:
                line.append('.')
        line = '\t'.join(str(x) for x in line)
        o.write(line + '\n')
    logger.info("Conversion complete!")
