import logging


class VCFWriter:

    def __init__(self, filename, style='all'):
        self.filename = filename + '.' + style + '.vcf'
        self.fo = self._create_file()
        self.style = style

    def close(self):
        self.fo.close()

    def _create_file(self):
        """
        Write header

        :return: file object
        """
        f = open(self.filename, 'w')
        f.write('##fileformat=VCFv4.3\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Placeholder ofr a dummy genotype value">\n')
        f.write('##INFO=<ID=EventType,Number=1,Type=String,Description="Type of event">\n')
        f.write('##INFO=<ID=EventCount,Number=1,Type=Integer,Description="Number of events">\n')
        f.write('##INFO=<ID=RefCodon,Number=1,Type=String,Description="Reference Codon">\n')
        f.write('##INFO=<ID=AltCodon,Number=1,Type=String,Description="Alternate codon">\n')
        f.write('##INFO=<ID=AApos,Number=1,Type=Integer,Description="Amino Acid Position events">\n')
        f.write('##INFO=<ID=RefAA,Number=1,Type=String,Description="Reference amino acid">\n')
        f.write('##INFO=<ID=AltAA,Number=1,Type=String,Description="Alternate Amino Acid">\n')
        f.write("##contig=<ID=chr1,length=248956422>\n")
        f.write("##contig=<ID=chr2,length=242193529>\n")
        f.write("##contig=<ID=chr3,length=198295559>\n")
        f.write("##contig=<ID=chr4,length=190214555>\n")
        f.write("##contig=<ID=chr5,length=181538259>\n")
        f.write("##contig=<ID=chr6,length=170805979>\n")
        f.write("##contig=<ID=chr7,length=159345973>\n")
        f.write("##contig=<ID=chr8,length=145138636>\n")
        f.write("##contig=<ID=chr9,length=138394717>\n")
        f.write("##contig=<ID=chr10,length=133797422>\n")
        f.write("##contig=<ID=chr11,length=135086622>\n")
        f.write("##contig=<ID=chr12,length=133275309>\n")
        f.write("##contig=<ID=chr13,length=114364328>\n")
        f.write("##contig=<ID=chr14,length=107043718>\n")
        f.write("##contig=<ID=chr15,length=101991189>\n")
        f.write("##contig=<ID=chr16,length=90338345>\n")
        f.write("##contig=<ID=chr17,length=83257441>\n")
        f.write("##contig=<ID=chr18,length=80373285>\n")
        f.write("##contig=<ID=chr19,length=58617616>\n")
        f.write("##contig=<ID=chr20,length=64444167>\n")
        f.write("##contig=<ID=chr21,length=46709983>\n")
        f.write("##contig=<ID=chr22,length=50818468>\n")
        f.write("##contig=<ID=chrX,length=156040895>\n")
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
        return f

    def write(self, final_dict):
        """
        Write a VCF output
        :param final_dict:
            {'chr13:24:GCT:CAC':
            {'chr': 'chr13', 'pos': 24, 'ref': 'GCT', 'alt': 'CAC', 'RefAA': 'A', 'AltAA': 'H',
            'RefCodon': 'GCT', 'AltCodon': 'CAC', 'AApos': 8, 'EventType': 'Missense', 'count': 1}

        :return:
        """

        for k, v in final_dict.items():
            chrom = v['chr']
            pos = v['pos']
            _id = '.'
            ref = v['ref']
            alt = v['alt']
            qual = '.'
            _filter = '.'
            info = 'RefAA=' + v['RefAA'] + ';' + \
                   'AltAA=' + v['AltAA'] + ';' + \
                   'AApos=' + str(v['AApos']) + ';' + \
                   'RefCodon=' + v['RefCodon'] + ';' + \
                   'AltCodon=' + v['AltCodon'] + ';' + \
                   'AltCodon=' + v['AltCodon'] + ';' + \
                   'EventCount=' + str(v['count']) + ';' + \
                   'EventType=' + v['EventType'] + ';'
            format = 'GT'
            sample = '0/1'
            line = '\t'.join([str(x) for x in [chrom, pos, _id, ref, alt, qual, _filter, info, format, sample]])
            self.fo.write(line + '\n')

    def write_noncoding(self, input_dict, target_dict):
        """

        :param input_dict: dictionary describing the event
        :param target_dict: dictionary containing the chromosome name, since it's not included in the input_dict
        :return: None

        Example:
            input_dict ={'frameshift': False, 'pos': 179, 'ref': 'AAATTAA', 'alt': 'TGGTAAT',
                'event_type': 'DownstreamNonCoding'}
            CR.target_dict = {'chrom': 'chr13', 'from': 32363357, 'to': 32363533, 'size': 189, 'seq': 'ATTGAACTTACAGATGG
                GTGGTATGCTGTTAAGGCCCAGTTAGATCCTCCCCTCTTAGCTGTCTTAAAGAATGGCAGACTGACAGTTGGTCAGAAGATTATTCTTCATGGAGCAGAACTGGTGG
                GCTCTCCTGATGCCTGTACACCTCTTGAAGCCCCAGAATCTCTTATGTTAAAGGTAAATTAATTT',
                'mmc': [35, 122], 'codon_distance_up': 0, 'codon_distance_down': 0, 'strand': 'F'}


        Give these inputs, output a VCF file
        """
        row = [
            str(target_dict['chrom']),
            str(input_dict['pos']),
            '.',
            str(input_dict['ref']),
            str(input_dict['alt']),
            '.',
            '.',
            'EventType=' + input_dict['event_type'] + ';',
            'GT',
            '0/1'
        ]
        row = '\t'.join(row)
        self.fo.write(row + '\n')


class TSVWriter:

    def __init__(self, filename):
        self.filename = filename + '.codingdna.tsv'
        self.fo = self._create_file()

    def _create_file(self):
        """
        Write header

        :return: file object
        """
        f = open(self.filename, 'w')
        f.write(
            'Start\tEnd\tRef_DNA\tAltDNA\tRefAA\tAltAA\tConsecutive\tAApos\tAltAA_length\trefDNACodons\taltDNACodons\n')
        return f

    def close(self):
        self.fo.close()

    def write_result(self, input_dict):
        """

        :param input_dict:  {
            `pos_ref_alt`:  'start_pos_idx'_'stop_pos_idx'_'DNA_ref'_'DNA_alt',
            `RefAA`: what the reference Amino acid(s) should be
            `AltAA`: what the alternate Amino acid(s) are
            `cons`: A boolean as to whether or not all the changes were consecutive
        }

        :return: None
        """
        self.fo.write('\t'.join([str(x) for x in input_dict['pos_ref_alt']]))
        self.fo.write('\t' + input_dict['RefAA'])
        self.fo.write('\t' + input_dict['AltAA'])
        self.fo.write('\t' + input_dict['cons'].__str__())
        self.fo.write('\t' + input_dict['AApos'].__str__())
        self.fo.write('\t' + input_dict['AltAA'].__len__().__str__())
        self.fo.write('\t' + input_dict['refDNA'].__str__())
        self.fo.write('\t' + input_dict['altDNA'].__str__())
        self.fo.write('\n')


class MetricsWriter:

    def __init__(self, filename):
        self.filename = filename + '.metrics.tsv'
        self.fo = self._create_file()

    def close(self):
        self.fo.close()

    def _create_file(self):
        """
        Write header

        :return: file object
        """
        w = open(self.filename, 'w')
        info = '|sample_name    |total_count|off_target       |cigar_count        |invalid_count    |mm_count     |frameshift_count |nochange_count   |synonymous_count |missense_count   |stop_count   |mixed_count  |multi_syn  |usable_count   |'
        logging.info(f"{info}")
        w.write(info + '\n')
        info = '|---------------|-----------|-----------------|-------------------|-----------------|-------------|-----------------|-----------------|-----------------|-----------------|-------------|-------------|-----------|---------------|'
        logging.info(f"{info}")
        w.write(info + '\n')
        return w

    def close(self):
        self.fo.close()

    def write(self, info):
        self.fo.write(info + '\n')
