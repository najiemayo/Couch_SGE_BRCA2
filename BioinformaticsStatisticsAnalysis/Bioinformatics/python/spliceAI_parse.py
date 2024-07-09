#!~/python/3.7.5/bin/python

"""
This script takes in a vcf file and parses out
spliceAI annotation in the info field

INFO ANNOTATION
DS_AG Delta score (acceptor gain)
DS_AL Delta score (acceptor loss)
DS_DG Delta score (donor gain)
DS_DL Delta score (donor loss)
DP_AG Delta position (acceptor gain)
DP_AL Delta position (acceptor loss)
DP_DG Delta position (donor gain)
DP_DL Delta position (donor loss)

#New Annotation keys
SpliceAI_DS_AG
SpliceAI_DS_AL
SpliceAI_DS_DG
SpliceAI_DS_DL
SpliceAI_DP_AG
SpliceAI_DP_AL
SpliceAI_DP_DG
SpliceAI_DP_DL

"""

import argparse

from cyvcf2 import VCF


def main():
    args = parse_args()
    splice_parse(args.vcf_file, args.out_file)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', dest='vcf_file', required=True,
                        help="input vcf file")
    parser.add_argument('-o', dest='out_file', default=None,
                        help="output vcf file")
    args = parser.parse_args()
    return args


def splice_parse(vcf_file, out_file):
    """
    using cyvcf2 parse spliceAI fields into ; separated info fields
    INPUT: VCF file with spliceAI annotations
    OUTPUT: VCF file to stdout
    """
    # Read vcf with cyvcf2
    vcf = VCF(vcf_file)

    # Add info fields to vcf header
    vcf.add_info_to_header({'ID': "SpliceAI_ALLELE", 'Description': 'SpliceAI Allele',
                            'Type': 'String', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_GENE", 'Description': 'SpliceAI Gene',
                            'Type': 'String', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DS_AG", 'Description': 'SpliceAI Delta score (acceptor gain)',
                            'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DS_AL", 'Description': 'SpliceAI Delta score (acceptor loss)',
                            'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DS_DG", 'Description': 'SpliceAI Delta score (donor gain)',
                            'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DS_DL", 'Description': 'SpliceAI Delta score (donor loss)',
                            'Type': 'Float', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DP_AG", 'Description': 'SpliceAI Delta position (acceptor gain)',
                            'Type': 'Integer', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DP_AL", 'Description': 'SpliceAI Delta position (acceptor loss)',
                            'Type': 'Integer', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DP_DG", 'Description': 'SpliceAI Delta position (donor gain)',
                            'Type': 'Integer', 'Number': '1'})
    vcf.add_info_to_header({'ID': "SpliceAI_DP_DL", 'Description': 'SpliceAI Delta position (donor loss)',
                            'Type': 'Integer', 'Number': '1'})

    # remove existing spliceAI info field in header
    header_list = vcf.raw_header.split('\n')
    new_header_list = [i if not i.startswith('##INFO=<ID=SpliceAI,Number=.,') else False for n, i in
                       enumerate(header_list)]

    # print out new header
    if out_file:
        o = open(out_file, 'w')
    else:
        o = open(vcf_file.replace(".vcf", ".parsed.vcf"), 'w')

    for i in new_header_list:
        if i:
            o.write(str(i) + '\n')

    # iterating through each vcf record and separting the sliceAI INFO field into separate info annotations
    for variant in vcf:
        info_gen = variant.INFO
        variant_info_list = [i[0] for i in info_gen if "SpliceAI" not in i]
        try:
            value = variant.INFO.get('SpliceAI').split("|")
            ALLEL = "SpliceAI_ALLELE" + "=" + value[0]
            GENE = "SpliceAI_GENE" + "=" + value[1]
            DS_AG = "SpliceAI_DS_AG" + "=" + value[2]
            DS_AL = "SpliceAI_DS_AL" + "=" + value[3]
            DS_DG = "SpliceAI_DS_DG" + "=" + value[4]
            DS_DL = "SpliceAI_DS_DL" + "=" + value[5]
            DP_AG = "SpliceAI_DP_AG" + "=" + value[6]
            DP_AL = "SpliceAI_DP_AL" + "=" + value[7]
            DP_DG = "SpliceAI_DP_DG" + "=" + value[8]
            DP_DL = "SpliceAI_DP_DL" + "=" + value[9]
            out = ";".join(str(i) for i in [ALLEL,
                                            GENE,
                                            DS_AG,
                                            DS_AL,
                                            DS_DG,
                                            DS_DL,
                                            DP_AG,
                                            DP_AL,
                                            DP_DG,
                                            DP_DL])
            info_fields = ";".join(
                str(i) for i in [x + "=" + str(variant.INFO.get(x)) for x in variant_info_list] + [out])
            vcf_out = str(variant).strip('\n').split('\t')[:7] + [info_fields] + str(variant).strip('\n').split('\t')[
                                                                                 8:]
            vcf_out = "\t".join(str(i) for i in vcf_out)
            o.write(vcf_out + '\n')

        except AttributeError:
            # If SpliceAI did not provide an annotation
            vcf_out = str(variant).strip('\n').split('\t')[:7] + ['.'] + str(variant).strip('\n').split('\t')[8:]
            vcf_out = "\t".join(str(i) for i in vcf_out)
            o.write(vcf_out + '\n')


if __name__ == '__main__':
    main()
