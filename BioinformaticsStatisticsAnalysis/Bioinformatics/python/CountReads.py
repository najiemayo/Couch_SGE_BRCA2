import logging
from Bio.Seq import Seq
from dna_functions import position_wrapper


class CountReads:
    """
        self.read = 'M00738:252:000000000-JF4NM:1:1101:15954:14434\t0\tchr17\t56772411\t60\t120M31S\t*\t0\t0\t'\
                'GAACTTCTTGAGCAGGAGCATACCCAGGGCTTCATAATCACATTCTGTTCAGCACTAGATGATATTCTTGGGGGTGGAGTGCCCTTAATGAA'\
                'AACAACAGAAATTTGTGGTGTACCAGGTTAGTATTTTGTACTATCGTCAGGAAACCAAA\t]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]'\
                ']]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]#]]]]]]]]]]]]]]'\
                ']]]]]]]]]]]]]]]]]]]]]]]]\tNM:i:2\tMD:Z:41C70C7\tAS:i:110\tXS:i:0\tSA:Z:chr17,56772604,+,120S31M'\
                ',60,0;'.split('\t')
        self.indicator_dict = {'chrom': 'chr17', 'pos': 56772452, 'ref': 'C', 'required_alt': 'A'}
        self.target_dict = {'chrom': 'chr17', 'from': 56772452, 'to': 56772560, 'size': 108, 'seq': 'CTTCTGTTCAGCACTAGATGATATTCTTGGGGGTGGAGTGCCCTTAATGAAAACAACAGAAATTTGTGGTGCACCAGGTGTTGGAAAAACACAATTATGGTAAAATAAA'}
    """

    ############################################################
    #   SAM SPECIFICATION
    #   0. ReadName
    #   1. Flag
    #   2. Chrom
    #   3. Pos
    #   4. MAPQ
    #   5. CIGAR
    #   6. MateChrom
    #   7. MatePos
    #   8. TLEN (ignore)
    #   9. SEQ
    #   10. Quality
    ############################################################

    def __init__(self, read, indicator_dict, target_dict):
        # super().__init__()
        self.read = read
        self.indicator_dict = indicator_dict
        self.target_dict = target_dict
        self.results_dict = {}
        self.read_start = read[3]
        self.reference_aa, self.reference_dna, self.other = self.get_aa()
        self.aa_list = self.build_aa_ref_list()

    def get_aa(self, reference=None):
        strand = self.target_dict.get('strand')
        if reference is None:
            read_start = int(self.read_start)
            cds_full_start = int(self.target_dict['from']) - read_start
            cds_full_stop = int(self.target_dict['to']) - read_start
            reference = self.target_dict.get('seq')

            # Get all upstream sequences (Left side)
            """
            Examples:
                Assumption1:    upstream_seq = 'ATTTGTCCAG',            cds_up = 0
                    Results:    upstream_noncoding_seq = 'ATTTGTCCAG',  cds_up_seq = ''
                Assumption2:    upstream_seq = 'ATTTGTCCAG',            cds_up = 1
                    Results:    upstream_noncoding_seq = 'ATTTGTCCA',  cds_up_seq = 'G'
            """
            cds_up = self.target_dict['codon_distance_up']
            upstream_seq = reference[:cds_full_start]
            upstream_noncoding_seq = upstream_seq[:upstream_seq.__len__() - cds_up]  # All non-coding
            cds_up_seq = upstream_seq[upstream_seq.__len__() - cds_up:]  # Partial coding regions

            # Get all downstream sequences (Right side)
            """
            Examples:
                Assumption1:    downstream_seq = ''GTATGATGTAT'',         cds_down = 0
                    Results:    downstream_noncoding_seq = 'GTATGATGTAT', cds_down_seq = ''
                Assumption2:    downstream_seq = 'GTATGATGTAT',           cds_down = 1
                    Results:    downstream_noncoding_seq = 'GTATGATGTAT', cds_down_seq = 'G'
            """
            cds_down = self.target_dict['codon_distance_down']
            cds_down_seq = reference[cds_full_stop:cds_full_stop + cds_down]  # Partial coding regions
            downstream_noncoding_seq = reference[cds_full_stop + cds_down + 1:]  # All non-coding

            reference = reference[cds_full_start:cds_full_stop + 1]
            reference = Seq(reference)
            cds_down_seq = Seq(cds_down_seq)
            downstream_noncoding_seq = Seq(downstream_noncoding_seq)
            cds_up_seq = Seq(cds_up_seq)
            upstream_noncoding_seq = Seq(upstream_noncoding_seq)
            reference_seq_fstrand = reference
            if reference.__len__() % 3 != 0:
                raise Exception(f"reference: ({reference}) is not divisible by 3 ({reference.__len__() % 3})")
            if strand == 'R':
                reference_seq = reference.reverse_complement().transcribe().translate()
            else:
                reference_seq = reference.transcribe().translate()

            other_seq_dict = {'downstream_noncoding_seq': downstream_noncoding_seq,
                              'cds_down_seq': cds_down_seq,
                              'reference_seq_fstrand': reference_seq_fstrand,
                              'cds_up_seq': cds_up_seq,
                              'upstream_noncoding_seq': upstream_noncoding_seq
                              }
            return reference_seq, reference.upper(), other_seq_dict

    def validateSAMLine(self):
        """
        Ensures that all elements of the SAM file are present
        :return: bool
        """
        if self.read[2] == self.target_dict.get('chrom') and \
                self.read.__len__() >= 10:
            """"
            int(self.read[3]) + len(self.read[9]) > self.target_dict.get('from') and \
            int(self.read[3]) < self.target_dict.get('to') and 
            """
            return True
        else:
            return False

    def is_correct_base(self, seq_base):
        return self.indicator_dict.get('required_alt') == seq_base

    def trim_to_target(self):
        """
        Given a SAM line, return the string that starts with the editable sequence
        Only returns up to a string that is possible to be in the target area

        :return: trimmed read, upstream seq, downstream seq
        """
        upstream_trim = self.target_dict['from'] - int(self.read_start)
        downstream_trim = self.target_dict['to'] - int(self.read_start) + 1
        tmp_read = self.read[9]
        if upstream_trim == 0:
            return tmp_read[:downstream_trim]
        else:
            return tmp_read[upstream_trim:downstream_trim]

    @staticmethod
    def check_validity(mmc, target_dict, strict=True):
        """
        Checks whether all the required changes are present in the read
        :param mmc: index position of mismatches
        :param target_dict: target dictionary
        :param strict: require all inserted bases
        :return: bool
        """
        # There can be others, but we need to ensure that the exact changes are incorporated
        if strict:
            if all(str(x) in mmc for x in target_dict['mmc']):
                return True
            else:
                return False

        # Only used if strict is false
        if any(str(x) in mmc for x in target_dict['mmc']):
            return True
        else:
            return False

    def get_mismatch_counts(self, full_read, aa=False):
        """
        trimmed_read: read that I want to compare with reference
        full_read: unmodified read from SAM file
        aa: False if DNA comparison is to be made, otherwise it goes into amino acid comp

        returns:
            If aa=False, returns a Tuple ([mismatch indexes], {empty})
            If aa=True, returns a Tuple ([mismatch indexes], {last variant})
        """
        alternate_dnaseq = full_read
        reference_dnaseq = self.target_dict['seq']

        # These should always be the same length
        dna_min_val = min(alternate_dnaseq.__len__(), reference_dnaseq.__len__())
        dna_mismatches = [str(i) for i in range(dna_min_val) if alternate_dnaseq[i] != reference_dnaseq[i]]

        if aa is False:
            return dna_mismatches, {}

        # Trim reference and alternate to be of the same length
        strand = self.target_dict.get('strand')

        reference_dnaseq = self.other['reference_seq_fstrand'].__str__()
        fro = self.target_dict['from'] - int(self.read_start)
        to = self.target_dict['to'] - int(self.read_start) + 1
        alternate_dnaseq = Seq(full_read[fro:to]).__str__()

        # All analyses below are relative to amino acids
        if strand == 'F':
            # Alt is set to false to indicate there is a stop codon
            alt_AAseq = Seq(alternate_dnaseq).transcribe().translate(to_stop=False).__str__()
            reference_AAseq = Seq(reference_dnaseq).transcribe().translate(to_stop=True).__str__()
        else:
            alt_AAseq = Seq(full_read[fro:to]).reverse_complement().transcribe().translate()
            ref_AASeq = Seq(self.other['reference_seq_fstrand'].__str__()).reverse_complement().translate()

        mm_dict = dict()
        # Remove the expected DNA mutations
        for x in self.target_dict['mmc']:
            if str(x) in dna_mismatches:
                dna_mismatches.remove(str(x))

        if dna_mismatches.__len__() == 0:
            return [], {'frameshift': False,
                        'event_type': 'NoChanges',
                        'pos': 0,
                        'ref': '0',
                        'alt': '0'}
        tmp_dict = position_wrapper(full_read, self.target_dict, int(self.read_start), self.other, strand)
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        """
        Overwrite variables
        """
        reference_AAseq = tmp_dict['RefAA']
        alt_AAseq = tmp_dict['AltAA']
        reference_dnaseq = tmp_dict['refDNA']
        alternate_dnaseq = tmp_dict['altDNA']
        aa_mismatches = [str(i + tmp_dict['AApos']) for i in range(alt_AAseq.__len__()) if
                         alt_AAseq[i] != reference_AAseq[i]]
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # If they are equal, then there is no change. However, I want to know if it was synonymous or unusable
        if (alt_AAseq == reference_AAseq and alternate_dnaseq == reference_dnaseq) or alt_AAseq == '':
            dna_mismatches = [int(x) for x in dna_mismatches]
            min_pos = min(dna_mismatches)
            max_pos = max(dna_mismatches)

            ncd_len = self.other['downstream_noncoding_seq'].__str__().__len__()
            cdd_len = self.other['cds_down_seq'].__str__().__len__()
            c___len = self.other['reference_seq_fstrand'].__str__().__len__()
            cdu_len = self.other['cds_up_seq'].__str__().__len__()
            ncu_len = self.other['upstream_noncoding_seq'].__str__().__len__()

            ref = self.target_dict['seq'][min_pos:max_pos + 1]
            alt = full_read[min_pos:max_pos + 1]
            res_dict = {'frameshift': False,
                        'pos': min_pos,
                        'ref': ref,
                        'alt': alt
                        }
            if max_pos < ncu_len:
                # All changes are in upstream non-coding sequences
                res_dict['event_type'] = 'UpstreamNoncoding'
                return [], res_dict

            if max_pos <= ncu_len + cdu_len:
                if min_pos >= ncu_len:
                    res_dict['event_type'] = 'PartialCodingUpStream'
                else:
                    res_dict['event_type'] = 'PartialCodingAndUpstreamNoncoding'
                return [], res_dict

            if min_pos > ncu_len + cdu_len + c___len + cdd_len:
                # All changes are in downstream non-coding sequences
                res_dict['event_type'] = 'DownstreamNonCoding'
                return [], res_dict

            if max_pos >= ncu_len + cdu_len + c___len:
                if max_pos < ncu_len + cdu_len + c___len + cdd_len:
                    res_dict['event_type'] = 'PartialCodingDownStream'
                else:
                    res_dict['event_type'] = 'DownstreamNonCodingAndPartialCodingDownstream'
                return [], res_dict

            if ncu_len > min_pos and max_pos > ncu_len:
                res_dict['event_type'] = 'Mixed'
                return [], res_dict

            raise Exception("This stage should never have been reached")

        if alt_AAseq.find('*') >= 0:
            # Anything here must be a frameshift or a stop codon
            # If there is only 1 difference in amino acids, then its a stop codon

            if aa_mismatches.__len__() > 1:
                # Frameshift
                mm_dict['pos'] = aa_mismatches  # The first position where they differ
                mm_dict['ref'] = tmp_dict['RefAA']
                mm_dict['alt'] = tmp_dict['AltAA']
                mm_dict['event_type'] = 'Frameshift'
                mm_dict['frameshift'] = True
                return aa_mismatches, mm_dict

            elif aa_mismatches.__len__() == 0:
                # This is a multiple synonymous event that just happens to have a stop codon
                return aa_mismatches, {'frameshift': False,
                                       'event_type': 'MultipleSynonymous',
                                       'pos': tmp_dict['AApos'],
                                       'ref': reference_AAseq,
                                       'alt': reference_AAseq
                                       }

            else:
                # Stop codon
                mm_dict['pos'] = tmp_dict['AApos']
                mm_dict['ref'] = tmp_dict['RefAA']
                mm_dict['alt'] = tmp_dict['AltAA']
                mm_dict['event_type'] = 'StopGain'
                mm_dict['frameshift'] = False
                return aa_mismatches, mm_dict
        else:
            """
            Remaining Possibilities
                Missense, 
                Mixed,
                Synonymous, 
                Frameshift that is non-truncating, 
                Too many
                
            Non-truncating Frameshift and TooMany are indistinguishable
            """

            # If the amino acid sequences are the same, then they are synonymous
            if alt_AAseq == reference_AAseq:
                # There can only be synonymous changes
                if aa_mismatches.__len__() > 1:
                    # IN this case there are multiple synonymous variants
                    aas = [int(x) for x in aa_mismatches]
                    reference_AAseq = [reference_AAseq[x - tmp_dict['AApos']] for x in
                                       aas]  # [reference_AAseq[x] for x in aas]
                    return aa_mismatches, {'frameshift': False,
                                           'event_type': 'MultipleSynonymous',
                                           'pos': aas,
                                           'ref': reference_AAseq,
                                           'alt': reference_AAseq
                                           }

                if tmp_dict['RefAA'].__len__() > 1:
                    # This is truly mixed
                    mm_dict = {'frameshift': False,
                               'event_type': 'MultipleSynonymous',
                               'pos': tmp_dict['AApos'],
                               'ref': tmp_dict['RefAA'],  # reference_AAseq[int(aa_mismatches[0])],
                               'alt': tmp_dict['AltAA']  # reference_AAseq[int(aa_mismatches[0])]
                               }
                    return [tmp_dict['AApos'].__str__()], mm_dict

                # This is a synonymous change
                mm_dict = {'frameshift': False,
                           'event_type': 'Synonymous',
                           'pos': tmp_dict['AApos'],
                           'ref': reference_AAseq,  # reference_AAseq[int(aa_mismatches[0])],
                           'alt': reference_AAseq  # reference_AAseq[int(aa_mismatches[0])]
                           }
                return [tmp_dict['AApos'].__str__()], mm_dict

            """Remaining Possibilities
                            Missense, 
                            Mixed, 
                            Frameshift that is non-truncating, 
                            Too many
                        Non-truncating Frameshift and TooMany are indistinguishable
            """

            if aa_mismatches.__len__() > 1:
                # There are > 1 amino acids changed
                # Too many changes are made to use
                return aa_mismatches, {'frameshift': False,
                                       'event_type': 'TooMany',
                                       'pos': aa_mismatches,
                                       'ref': tmp_dict['RefAA'],
                                       'alt': tmp_dict['AltAA']}

            """Remaining Possibilities
                            Missense, 
                            Mixed
            """
            # separate DNA matches that change AA from those that don't
            # At this point we know that the aadiff == 1
            if aa_mismatches.__len__() == 1 and tmp_dict['RefAA'].__len__() == 1:
                # Then there is only 1 change, and it is a missense
                return aa_mismatches, {'frameshift': False,
                                       'event_type': 'Missense',
                                       'pos': tmp_dict['AApos'],
                                       'ref': tmp_dict['RefAA'],
                                       'alt': tmp_dict['AltAA']}
            """Remaining Possibilities 
                            Mixed            
            aa_mismatches is an index of amino acids that correlate to DNA changes
            aa_diff is a list of indexes where the amino acids are different
            
            So, we need to remove the amino acid change first, then treat the rest as synonymous
            
            {'pos_ref_alt': [2, 6, 'TGAA', 'AGAC'], 'RefAA': 'IE', 'AltAA': 'ID', 'cons': False, 'AApos': 0, 'refDNA': 'ATTGAA', 'altDNA': 'ATAGAC'}
            So in this case, the first mutatio is Syn, so we'll exclude that.  Should return
            
            """

            amino_acid_that_changed = aa_mismatches[0]
            try:
                aa_mismatches.remove(str(amino_acid_that_changed))
            except ValueError:
                # This is a rare case where there is an amino acid change at the fixed nucleotide position
                return aa_mismatches, {'frameshift': False,
                                       'event_type': 'Mixed',
                                       'pos': tmp_dict['AApos'] + int(amino_acid_that_changed),
                                       'ref': tmp_dict['RefAA'][int(amino_acid_that_changed)],
                                       'alt': tmp_dict['AltAA'][int(amino_acid_that_changed)]}

            # TODO: Add test case for reading through an exon into an intron
            # if int(amino_acid_that_changed) >= self.reference_aa.__len__():
            # In this case the exon has stopped and we are now in the intron at the end of the transcript.

            return aa_mismatches, {'frameshift': False,
                                   'event_type': 'Mixed',
                                   'pos': int(amino_acid_that_changed),
                                   'ref': tmp_dict['RefAA'][tmp_dict['AApos'] - int(amino_acid_that_changed)],
                                   'alt': tmp_dict['AltAA'][tmp_dict['AApos'] - int(amino_acid_that_changed)]}

    @staticmethod
    def get_mismatch_bases(read, positions):
        return [read[int(x)] for x in positions]

    def add_to_dict(self, read, aln_start):

        # Count every base in every read
        a = self.target_dict['seq'].__len__()
        b = read.__len__()
        for i in range(min(a, b)):
            seq_base = read[i]
            ref_base = self.target_dict['seq'][i]
            if str(i) not in self.results_dict.keys():
                self.results_dict[str(i)] = {'POS': aln_start + i,
                                             'REF': ref_base,
                                             'A': 0, 'G': 0, 'T': 0, 'C': 0,
                                             'N': 0}

            self.results_dict[str(i)][seq_base] += 1

    def add_aa_dict(self, aa_dict):
        """
        :param aa_dict: position, reference, and alternate amino acid

        Example:
               aa_dict = {'frameshift': False, 'event_type': 'Missense', 'pos': 8, 'ref': 'N', 'alt': 'H'}
        """
        if aa_dict['event_type'] == 'MultipleSynonymous':
            for a in range(aa_dict['ref'].__len__()):
                index = aa_dict['pos'] + a
                alt_aa = aa_dict['alt'][a]
                alt_aa = alt_aa.replace('*', "X")
                self.aa_list[index]['TOTAL'] += 1
                self.aa_list[index]['changes'][alt_aa] += 1
            return

        # if aa_dict['event_type'] == 'Mixed':
        #     "MIXED classes have a special field (other) that are synonymous. After which you treat them as WT"
        #     for i in aa_dict['other']:
        #         ref_aa = self.aa_list[int(i)]['REF']
        #         self.aa_list[int(i)]['changes'][ref_aa] += 1
        # Increment the total amino acid and the alternate
        index = aa_dict['pos']
        alt_aa = aa_dict['alt']
        alt_aa = alt_aa.replace('*', "X")
        try:
            self.aa_list[index]['TOTAL'] += 1
            self.aa_list[index]['changes'][alt_aa] += 1
        except IndexError:
            logging.warning(f"Index {index} exceeds the expected amino acid length.")

    def build_aa_ref_list(self):
        """ Initialize the reference dictionary with all possible amino acids"""
        a, d, other = self.get_aa()
        amino_acid_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                            'W', 'Y', 'X']
        res = []

        for i in range(a.__len__()):
            p1 = {'POS': i, 'REF': a[i], 'TOTAL': 0, 'changes': dict()}
            for j in amino_acid_codes:
                p1['changes'][j] = 0
            res.append(p1)
        return res
