import argparse
import logging
import sys


# noinspection DuplicatedCode
def parse_args():
    """
    target_dict = {'chrom': 'chr17', 'from': 56774036,
                   'seq': 'TCATCTTTCTGTTGACAGTATGCAGTTGGCAGTAGATGTGCAGATACCAGAATGTTTTGGAGGAGTGGCAGGTGAAGCAGTTTTTATTGATACAGAGGG'}
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--chrom",
                        dest='chrom',
                        required=True,
                        help="Chromosome of target")

    parser.add_argument("-f", "--from",
                        dest='from',
                        type=int,
                        required=True,
                        help="Beginning coordinate of sequence")

    parser.add_argument("-s", "--seq",
                        dest='seq',
                        required=True,
                        help="Sequence of target")

    parser.add_argument("-o", "--outfile",
                        dest='outfile',
                        required=True,
                        help="File name of output")

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


class makeAltVCF:

    def __init__(self, target_dict):
        self.target_dict = target_dict

    @staticmethod
    def _altGenerator(base):
        bases = ["A", "C", "T", 'G']
        for b in bases:
            if b != base:
                yield b

    def getline(self):
        k = self.target_dict.get('seq')
        k = list(k).__len__()
        for i in range(k):
            chrom = self.target_dict.get('chrom')
            pos = int(self.target_dict.get('from')) + i
            ref_base = self.target_dict.get('seq')
            ref_base = ref_base[i]
            for alt in self._altGenerator(ref_base):
                yield chrom, pos, '.', ref_base, alt

    @staticmethod
    def write_vcf(outfile_name):
        """
        Writes the VCF file of all possible alternates
        :param outfile_name:  Name of VCF file to create
        :return: None
        """
        with open(outfile_name, 'w') as r:
            r.write("##fileformat=VCFv4.1\n")
            r.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            r.write("##contig=<ID=chr1,length=248956422>\n")
            r.write("##contig=<ID=chr2,length=242193529>\n")
            r.write("##contig=<ID=chr3,length=198295559>\n")
            r.write("##contig=<ID=chr4,length=190214555>\n")
            r.write("##contig=<ID=chr5,length=181538259>\n")
            r.write("##contig=<ID=chr6,length=170805979>\n")
            r.write("##contig=<ID=chr7,length=159345973>\n")
            r.write("##contig=<ID=chr8,length=145138636>\n")
            r.write("##contig=<ID=chr9,length=138394717>\n")
            r.write("##contig=<ID=chr10,length=133797422>\n")
            r.write("##contig=<ID=chr11,length=135086622>\n")
            r.write("##contig=<ID=chr12,length=133275309>\n")
            r.write("##contig=<ID=chr13,length=114364328>\n")
            r.write("##contig=<ID=chr14,length=107043718>\n")
            r.write("##contig=<ID=chr15,length=101991189>\n")
            r.write("##contig=<ID=chr16,length=90338345>\n")
            r.write("##contig=<ID=chr17,length=83257441>\n")
            r.write("##contig=<ID=chr18,length=80373285>\n")
            r.write("##contig=<ID=chr19,length=58617616>\n")
            r.write("##contig=<ID=chr20,length=64444167>\n")
            r.write("##contig=<ID=chr21,length=46709983>\n")
            r.write("##contig=<ID=chr22,length=50818468>\n")
            r.write("##contig=<ID=chrX,length=156040895>\n")
            r.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            for j in q.getline():
                j = list(j)
                j.extend(['.', '.', '.', 'GT', "0/1"])
                j = '\t'.join([str(x) for x in j])
                r.write(j + "\n")


if __name__ == '__main__':
    args, logger = parse_args()
    q = makeAltVCF(args.__dict__)
    q.write_vcf(args.outfile)
