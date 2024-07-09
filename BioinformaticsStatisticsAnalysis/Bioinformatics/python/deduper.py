def add_to_dict(final_dict, strand, mmmc_dict, dna_results_dict=None, chrom=None, read_start=None):
    """
    Create a non-redundant set of annotations

    For non-coding variants, the dna_results_dict will be none

    final_dict

    Amino acid-based
    mmmc_dict = {'frameshift': False, 'event_type': 'Missense', 'pos': 8, 'ref': 'A', 'alt': 'H'}
    or
    DNA-based
    mmmc_dict = {'frameshift': False, 'pos': 183, 'ref': 'TAA', 'alt': 'GGG', 'event_type': 'DownstreamNonCoding'}

    Only present when there is a coding mutation
    dna_results_dict = {'pos_ref_alt': [89, 89, 'T', 'G'], 'RefAA': 'G', 'AltAA': 'G', 'cons': True, 'AApos': 29, 'refDNA': 'GGT', 'altDNA': 'GGG'}
    """

    # Create a dictionary key if it doesn't exist
    if dna_results_dict is not None:
        # This is for coding mutations
        dict_key = str(chrom) + ":" + \
                   str(dna_results_dict['pos_ref_alt'][0]) + ":" + \
                   str(dna_results_dict['pos_ref_alt'][2]) + ":" + \
                   str(dna_results_dict['pos_ref_alt'][3])
        chr = str(chrom)
        pos = read_start + int(dna_results_dict['pos_ref_alt'][0])
        ref = dna_results_dict['pos_ref_alt'][2]
        alt = dna_results_dict['pos_ref_alt'][3]
        RefAA = dna_results_dict['RefAA']
        AltAA = dna_results_dict['AltAA']
        RefCodon = dna_results_dict['refDNA']
        AltCodon = dna_results_dict['altDNA']
        AApos = dna_results_dict['AApos']
        EventType = mmmc_dict['event_type']
    else:
        # This is a non-coding mutation
        dict_key = str(chrom) + ":" + \
                   str(mmmc_dict['pos']) + ":" + \
                   str(mmmc_dict['ref']) + ":" + \
                   str(mmmc_dict['alt'])
        chr = str(chrom)
        pos = read_start + int(mmmc_dict['pos'])
        ref = mmmc_dict['ref']
        alt = mmmc_dict['alt']
        RefAA = '.'
        AltAA = '.'
        RefCodon = '.'
        AltCodon = '.'
        AApos = '.'
        EventType = mmmc_dict['event_type']

    if dict_key not in final_dict.keys():
        final_dict[dict_key] = dict()
        final_dict[dict_key]['count'] = 0

    new_dict = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt,
                'RefAA': RefAA, 'AltAA': AltAA,
                'RefCodon': RefCodon, 'AltCodon': AltCodon,
                'AApos': AApos, 'EventType': EventType, 'count': final_dict[dict_key]['count']}
    final_dict[dict_key] = new_dict
    final_dict[dict_key]['count'] += 1
    return final_dict
