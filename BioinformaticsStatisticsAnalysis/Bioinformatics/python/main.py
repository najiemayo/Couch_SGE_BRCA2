"""Console script for Site saturation mutagenesis experiments."""
import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from CountReads import CountReads
from parse_cigar import split_cigar, filter_matches
import dna_functions as dnaf
from output import VCFWriter, TSVWriter, MetricsWriter
from deduper import add_to_dict


# noinspection DuplicatedCode
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--sam",
                        dest='samfile',
                        required=True,
                        help="SAM file containing reads")

    parser.add_argument("-c", "--chrom",
                        dest='chrom',
                        required=True,
                        help="Chromosome of target")

    parser.add_argument("-p", "--positions",
                        dest='positions',
                        required=True,
                        type=int,
                        nargs='+',
                        help="locations of the required changes")

    parser.add_argument("-m", "--cigar_min",
                        dest='cigar_min',
                        default=60,
                        type=int,
                        help="required length contiguous bases in cigar string")

    parser.add_argument("-b", "--bases",
                        dest='bases',
                        required=True,
                        choices=['A', 'C', 'G', 'T'],
                        nargs='+',
                        help="Non reference bases that must be located at positions")

    parser.add_argument("-r", "--ref",
                        dest='seq',
                        required=True,
                        help="Target sequence (NOT containing mutations!")

    parser.add_argument("-P", "--target_bounds",
                        dest='target_bounds',
                        required=True,
                        nargs=2,
                        type=int,
                        help="Start and stop [low, high] of the full codon region to assess "
                             "(Do not include partial codons here!)")

    parser.add_argument("-z", "--strand",
                        dest='strand',
                        required=True,
                        choices=['F', 'R'],
                        help="Forward (F) or Reverse (R) strand")

    parser.add_argument("-d", "--codon_distance_up",
                        dest='codon_distance_up',
                        default=0,
                        type=int,
                        help="Number of bases in incomplete codon in the upstream direction (smaller positions)")

    parser.add_argument("-D", "--codon_distance_down",
                        dest='codon_distance_down',
                        default=0,
                        type=int,
                        help="Number of bases in incomplete codon in the downstream direction (larger positions)")

    parser.add_argument("-S", "--read_start",
                        dest='read_start',
                        required=True,
                        type=int,
                        help="Where the read begins (relative to the forward strand)")

    parser.add_argument("-x", "--strict",
                        dest='strict',
                        default=True,
                        type=bool,
                        help="Set this flag if you want to require only 1 required base change instead of all")

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


class SamParser:

    def __init__(self, args, log_level='INFO'):
        self.args = args
        self.log_level = log_level
        self.logger = None
        self.read_start = args['read_start']
        self.target_dict = self.build_target_dict(args)
        self.indicator_dict = self.build_indicator_dict(args)
        self.validate_args()
        self.sample_name = os.path.basename(args['samfile']).replace(".sam", "")
        self.set_logger()
        self.mm_count = 0
        self.invalid_count = 0
        self.total_count = 0
        self.usable_count = 0
        self.frameshift_count = 0
        self.nochange_count = 0
        self.synonymous_count = 0
        self.missense_count = 0
        self.stop_count = 0
        self.cigar_count = 0
        self.off_target = 0
        self.mixed_count = 0
        self.multi_synonymous_count = 0
        self.noncoding_count = 0
        self.cigar_min = args['cigar_min']
        self.strict = args['strict']

    def validate_args(self):
        assert os.path.exists(self.args['samfile']), "Must provide args.samfile"
        assert 'chrom' in self.args.keys(), "Must provide args.chrom"
        assert 'positions' in self.args.keys(), "Must provide args.positions"
        for k in self.args['positions']:
            assert type(k) == int, "Positions must be integers"
        assert 'bases' in self.args.keys(), "Must provide args.bases"
        for k in self.args['bases']:
            assert k.upper() in ['A', 'C', 'T', 'G'], "args.bases must only be of A, C, T, or G"
        assert 'seq' in self.args.keys()
        for k in self.args['seq']:
            assert k.upper() in ['A', 'C', 'T', 'G'], "args.seq must only be of A, C, T, or G"
        return True

    def set_logger(self):
        logging.basicConfig(stream=sys.stderr, level=self.log_level,
                            format='%(name)s (%(levelname)s): %(message)s')

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(self.log_level)

    @staticmethod
    def build_target_dict(args):
        """
        Build a target dictionary

        :param args:
            {
                'samfile': SAM file containing reads,
                'chrom': chromsome of target,
                'positions': list of positions where mutations should be,
                'bases': list of bases corresponding to each position,
                'seq': sequence of your target (with desired mutation),
                'target_bounds': coordinates of your target,
                'logLevel': level of logging
            }
        :return: {
                    'chrom': 'chr17',
                    'from': 56774036,
                    'to': 56774134,
                    'size': 99,
                    'mmc': 2,
                    'seq': 'TCATCTTTCTGTTGACAGTATGCAGTTGGCAGT
                                'AGATGTGCAGATACCAGAATGTTTTGGAGGAGT'
                                'CGCAGGTGAAGCAGTTTTTATTGATACAGAAGG'
                    'strand': 'F',
                    'codon_dist': 0
                }
        """
        target_dict = dict()
        target_dict['chrom'] = args['chrom']
        target_dict['from'] = int(args['target_bounds'][0])
        target_dict['to'] = int(args['target_bounds'][1])
        target_dict['size'] = args['seq'].__len__()
        target_dict['seq'] = args['seq'].upper()
        # Expected mismatch indexes in seq
        target_dict['mmc'] = [x - int(args['read_start']) for x in args['positions']]
        target_dict['codon_distance_up'] = args['codon_distance_up']
        target_dict['codon_distance_down'] = args['codon_distance_down']

        target_dict['strand'] = args['strand']

        return target_dict

    @staticmethod
    def build_indicator_dict(args):
        """
        Build an array of indicators that need to be modified
        :param args:
            {
                'samfile': SAM file containing reads,
                'chrom': chromosome of target,
                'positions': list of positions where mutations should be,
                'bases': list of bases corresponding to each position,
                'seq': sequence of your target,
                'target_bounds': coordinates of your target,
                'logLevel': level of logging
            }
        :return:
        [
            {'chrom': 'chr17', 'pos': 56774102, 'ref': 'G', 'required_alt': 'C'},
            {'chrom': 'chr17', 'pos': 56774132, 'ref': 'G', 'required_alt': 'A'}
        ]
        """
        indicator_dict = list()

        for x in range(args['positions'].__len__()):
            tmpdict = dict()
            tmpdict['chrom'] = args['chrom']
            tmpdict['pos'] = int(args['positions'][x])
            tmpdict['required_alt'] = args['bases'][x]
            indicator_dict.append(tmpdict)
            del tmpdict
        return indicator_dict

    # noinspection DuplicatedCode
    def process_sam(self):
        """
        Main function for parsing SAM files for SSM experiments

        :return: dictionary containing
        {
            'results': Total_CR, # A dictionary with additional dictionaries:
                        indicator dict: chrom, pos, required alt that was searched
                        target_dict: chrom, start, stop, target sequence
                        results_dict: each key is the index posiotion of the target sequence while each entry
                            contains the position, reference base, number of each nucleotides that met the criteria,
                            as well as number of reads that failed because they didn't match the required alt base(s)
            'sample_name':  sample_name,
            'total_count':  total_count, Total number of reads in the library
            'off_target': Reads that didn't start in the right posiiton
            'cigar_count':  Fails CIGAR filter
            'mm_count':     Has too many amino acid changes
            'invalid_count': Required DNA changes not present
            'frameshift_count':  How many frameshifts
            'nochange_count': How many with no changes
            'synonymous_count': How many synonymous?
            'missense_count': How many missense
            'stop_count': How many stop-gained mutations
            'mixed_count': How many Synonmymous-Missense hybrid reads
            'multi_syn': How many reads with multiple synonymous variants were found
            'usable_count':  synonymous_count + missense_count + stop_count

        }
        """
        sample_name = os.path.basename(self.args['samfile']).replace('.sam', '').replace('.gz', '')

        dna_out = TSVWriter(sample_name)
        metrics = MetricsWriter(sample_name)
        f = VCFWriter(sample_name, style='all')
        final_dict = dict()

        i = 0
        with open(self.args['samfile'], 'r') as r:
            # Put in dummy values as a placeholder
            read = 'ReadName\tFlag\tChrom\t0\tMAPQ\tCIGAR\tMateChrom\tMatePos\tTLEN\tSEQ\tQuality'.split('\t')
            Total_CR = CountReads(read, self.indicator_dict, self.target_dict)
            for row in r:
                i += 0
                # Parse the files, skipping over headers
                if row.startswith('@'):
                    continue
                # Now we're counting reads
                line = row.strip().split('\t')
                self.total_count += 1

                # Make sure the read starts exactly in the same place that you're expecting one to be
                if int(line[3]) != self.read_start or line[2] != self.target_dict['chrom']:
                    self.off_target += 1
                    continue

                # Skip if the CIGAR string is messy
                cig = filter_matches(split_cigar(line[5]), min_matches=self.cigar_min)
                if cig is False:
                    self.cigar_count += 1
                    continue

                # Create the model object for the SAM line
                CR = CountReads(line, self.indicator_dict, self.target_dict)

                # Ensure that all elements of the SAM file are present
                if CR.validateSAMLine():
                    # Limit the length of the read to the desired target sequence, and if necessary,
                    # skip the required number of bases until you hit the amino acid/frame of interest
                    # (Codon distance parameter) and [Downstream distance parameter]
                    # Remove variants that go beyond the coding sequence
                    # [Upstream distance parameter]
                    # target = CR.trim_to_target()
                    # Get index positions of DNA mismatches
                    mmc, _ = CR.get_mismatch_counts(line[9], aa=False)

                    # Make sure all required DNA changes are present
                    valid = CR.check_validity(mmc, self.target_dict, strict=self.strict)
                    if not valid:
                        # Not all the expected sequences were found
                        self.invalid_count += 1
                        mmc_len = ','.join(x for x in mmc)
                        continue

                    # Get amino acid-based results rather than DNA (if coding)
                    mmca, mmca_dict = CR.get_mismatch_counts(line[9], aa=True)

                    if mmca_dict['frameshift']:
                        self.frameshift_count += 1
                        continue

                    if mmca_dict['event_type'] == 'NoChanges':
                        self.nochange_count += 1
                        continue

                    if mmca_dict['event_type'] == 'TooMany':
                        self.mm_count += 1
                        continue

                    if mmca_dict['event_type'] == 'Mixed':
                        self.mixed_count += 1
                        continue

                    if mmca_dict['event_type'] == 'MultipleSynonymous':
                        self.multi_synonymous_count += 1
                        continue

                    # Everything below, I need an output for
                    if mmca_dict['event_type'] in ['UpstreamNoncoding', 'PartialCodingAndUpstreamNoncoding',
                                                   'PartialCodingUpStream', 'DownstreamNonCoding',
                                                   'DownstreamNonCodingAndPartialCodingDownstream',
                                                   'PartialCodingDownStream']:
                        self.noncoding_count += 1
                        final_dict = add_to_dict(final_dict,
                                                 CR.target_dict['strand'],
                                                 mmca_dict, dna_results_dict=None,
                                                 chrom=line[2],
                                                 read_start=self.read_start)
                        continue

                    if mmca_dict['event_type'] == 'Synonymous':
                        self.synonymous_count += 1
                        dna_results_dict = dnaf.position_wrapper(line[9], CR.target_dict, self.read_start, CR.other,
                                                                 self.args['strand'])
                        dna_out.write_result(dna_results_dict)
                        final_dict = add_to_dict(final_dict,
                                                 CR.target_dict['strand'],
                                                 mmca_dict, dna_results_dict=dna_results_dict,
                                                 chrom=line[2],
                                                 read_start=self.read_start)
                        continue

                    if mmca_dict['event_type'] == 'Missense':
                        self.missense_count += 1
                        dna_results_dict = dnaf.position_wrapper(line[9], CR.target_dict, self.read_start, CR.other,
                                                                 self.args['strand'])
                        dna_out.write_result(dna_results_dict)
                        final_dict = add_to_dict(final_dict,
                                                 CR.target_dict['strand'],
                                                 mmca_dict, dna_results_dict=dna_results_dict,
                                                 chrom=line[2],
                                                 read_start=self.read_start)
                        continue

                    if mmca_dict['event_type'] == 'StopGain':
                        self.stop_count += 1
                        ##TODO: if AA==2710 print out the read for day 14
                        dna_results_dict = dnaf.position_wrapper(line[9], CR.target_dict, self.read_start, CR.other,
                                                                 self.args['strand'])
                        AA_number=2660 + dna_results_dict['AApos']
                        if AA_number  == 2710:
                            print(line[9])
                            print(mmca_dict, mmca)
                        dna_out.write_result(dna_results_dict)
                        final_dict = add_to_dict(final_dict,
                                                 CR.target_dict['strand'],
                                                 mmca_dict, dna_results_dict=dna_results_dict,
                                                 chrom=line[2],
                                                 read_start=self.read_start)
                        continue

                    raise Exception(f"This should never be reached.\n"
                                    f"mmc: {mmc}, mmca: {mmca}, mmca_dict: {mmca_dict},\n"
                                    f"Line: {line}")

        self.usable_count = self.stop_count + self.missense_count + self.synonymous_count
        info = f'|{self.sample_name}|{self.total_count}|' \
               f'{self.off_target} ({self.off_target / self.total_count * 100:.2f}%)|' \
               f'{self.cigar_count} ({self.cigar_count / self.total_count * 100:.2f}%)|' \
               f'{self.invalid_count} ({self.invalid_count / self.total_count * 100:.2f}%)|' \
               f'{self.mm_count} ({self.mm_count / self.total_count * 100:.2f}%)|' \
               f'{self.frameshift_count} ({self.frameshift_count / self.total_count * 100:.2f}%)|' \
               f'{self.nochange_count} ({self.nochange_count / self.total_count * 100:.2f}%)|' \
               f'{self.synonymous_count} ({self.synonymous_count / self.total_count * 100:.2f}%)|' \
               f'{self.missense_count} ({self.missense_count / self.total_count * 100:.2f}%)|' \
               f'{self.stop_count} ({self.stop_count / self.total_count * 100:.2f}%)|' \
               f'{self.mixed_count} ({self.mixed_count / self.total_count * 100:.2f}%)|' \
               f'{self.multi_synonymous_count} ({self.multi_synonymous_count / self.total_count * 100:.2f}%)|' \
               f'{self.usable_count} ({self.usable_count / self.total_count * 100:.2f}%)| '
        self.logger.info(info)

        metrics.write(info)
        f.write(final_dict)

        f.close()
        dna_out.close()
        metrics.close()

        res = {
            'results': Total_CR,
            'sample_name': self.sample_name,
            'total_count': self.total_count,
            'cigar_count': self.cigar_count,
            'mm_count': self.mm_count,
            'invalid_count': self.invalid_count,
            'frameshift_count': self.frameshift_count,
            'nochange_count': self.nochange_count,
            'synonymous_count': self.synonymous_count,
            'missense_count': self.missense_count,
            'stop_count': self.stop_count,
            'usable_count': self.usable_count
        }
        return res


if __name__ == "__main__":
    args, logger = parse_args()
    args = args.__dict__
    SP = SamParser(args)
    res = SP.process_sam()
