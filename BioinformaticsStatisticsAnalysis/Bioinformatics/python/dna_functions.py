from Bio.Seq import Seq


def get_dna_changes(reference_seq, alt_seq):
    """
    Get a dictionary of changes and their locations

    :param reference_seq: Expected DNA sequence
    :param alt_seq: observed DNA sequence

    :returns: dictionary of changes [{'pos':135: 'ref': 'G', 'alt': 'A'}]

    """
    diff_list = list()
    dna_min_val = min(alt_seq.__len__(), reference_seq.__len__())
    reference_seq = reference_seq.upper()
    alt_seq = alt_seq.upper()
    dna_mismatches = [str(i) for i in range(dna_min_val) if alt_seq[i] != reference_seq[i]]
    for x in dna_mismatches:
        diff_dict = dict()
        diff_dict['pos'] = int(x)
        diff_dict['ref'] = reference_seq[int(x)]
        diff_dict['alt'] = alt_seq[int(x)]
        diff_list.append(diff_dict)
    return diff_list


def remove_expected_changes(dna_mismatches, target_mmc, codon_distance):
    """
    Removes expected DNA mismatches for mismatch list

    :param codon_distance:
    :param dna_mismatches: list of DNA mismatch index positions observed
    :param target_mmc: list of DNA mismatch index positions that one would expect to observe

    :return:
        list of DNA mismatch index positions observed (excluding the expected)
    """
    clean_mismatches = list()
    target_positions = [x + codon_distance for x in target_mmc]
    # Remove the expected DNA mutations
    for x in dna_mismatches:
        if x['pos'] not in target_positions:
            clean_mismatches.append(x)
    return clean_mismatches


def is_consecutive(mmc, max_distance=2):
    """
    Test to see if results are a consecutive set of changes

    :param mmc: list of dictionary elements of mismatches
    :param max_distance: how far apart should the variants be and still be called consecutive
                            (Since the library process uses 3 contiguous Ns to make mutations, we'll keep this at 3)

    :return:
        False if  [{'pos': 5, 'ref': 'A', 'alt': 'T'}, {'pos': 35, 'ref': 'A', 'alt': 'T'}]
        True if  [{'pos': 5, 'ref': 'A', 'alt': 'T'}, {'pos': 6, 'ref': 'A', 'alt': 'T'}]
    """
    deltas = list()
    mutated_positions = [x['pos'] for x in mmc]

    for x in range(mutated_positions.__len__() - 1):
        deltas.append(mutated_positions[x + 1] - mutated_positions[x])

    if sum(deltas) > max_distance:
        return False
    else:
        return True


def trim_contig(mmc, reference_seq, alt_seq):
    """
    Trim to the smallest contig that contains all the changes

    :param mmc: list of dictionary elements of mismatches
    :param reference_seq: Expected DNA sequence
    :param alt_seq: observed DNA sequence

    :returns: A sequence of [start_pos_idx, stop_pos_idx, ref, alt]
    """

    mutated_positions = [x['pos'] for x in mmc]
    # No change
    if mutated_positions.__len__() == 0:
        return [0, 0, None, None]

    # One Change
    if mutated_positions.__len__() == 1:
        return [mmc[0]['pos'], mmc[0]['pos'], mmc[0]['ref'], mmc[0]['alt']]

    min_index = min(mutated_positions)
    max_index = max(mutated_positions)

    trim_ref = reference_seq[min_index:max_index + 1]
    trim_alt = alt_seq[min_index:max_index + 1]
    return [min_index, max_index, trim_ref, trim_alt]


def get_affected_amino_acid(pos_ref_alt, seq, distance_to_change, other_dict, strand):
    """
    Given a set of inputs that are not relative to codons,
        1. Get their codon structure
        2. Get the amino acid changes

    :param distance_to_change:
    :param other_dict:
    :param strand:
    :param pos_ref_alt: a list of start, stop, refDNA, altDNA
    :param seq: The sequence of the target contig
    return refAA, altAA, refDNA, altDNA, AApos

    Example 1, Coding [34, 36, 'TT', 'AA']
        Coding starts at position 10
        If the first change happens before the coding sequence starts, then you can skip the codon stuff
        Nucleotides 34-36 belong to amino acid 8 [(34-10)//3 8 and (36-10)//3 8]
        Because (Start-CodonDistance) % 3 is 0, that means that all nucleotides at the beginning of the codon are accounted for
        Otherwise, on needs to get the remainder of bases from the target sequence (new_start = start - (Start-CodonDistance) % 3)
        Once at the beginning of the codon, add bases on the end from the sequence object until the number of nucleotides is divisible by 3

        For example, given [9, 11, 'GG','CC'] and the codon position is 10,
            the first G doesn't count, so the sequence should be trimmed to [10, 11, 'G', 'C']
    """

    start, stop, ref, alt = pos_ref_alt
    new_start = 0
    downstream_trim = 0
    # [36,36,'A','G'] -> ref A, alt G, seq, codon_index 0, bases to add to left side=0 AAT->GAT
    # [37,37,'A','G']->  ref A, alt G, seq, codon_index 0, bases to add to left side=1 AAT->AGT
    # [34,36,'CAA','CAG']->  ref A, alt G, seq, codon_index 0, bases to add to left side=1
    # [51, 53, 'CA', 'GT'] -> Ref CA, Alt GT, codon_index 0,bases to add to left side=0
    # [41, 41, 'A', 'G'] -> AAG->GAG, K -> E,
    # Where does the event start?
    # Remove any bases that are before the starting position
    # Where does the event end?
    # Remove any variants outside coding region
    unc_len = other_dict['upstream_noncoding_seq'].__str__().__len__()
    usc_len = other_dict['cds_up_seq'].__str__().__len__()
    dsc_len = other_dict['cds_down_seq'].__str__().__len__()
    dnc_len = other_dict['downstream_noncoding_seq'].__str__().__len__()
    us_len = unc_len + usc_len + other_dict['reference_seq_fstrand'].__str__().__len__()
    seq_len = other_dict['reference_seq_fstrand'].__str__().__len__()
    # syn1:[36, 36, 'A', 'G']
    # upp:pos_ref_alt [15, 36, 'ACCTAGAGGGAAAGCTTACCA', 'TCCTAGAGCGAAAGCTTACCT'] <- should be 35
    # up: pos_ref_alt [35, 35, 'A', 'T']
    # unc:            [15, 15, 'A', 'T']
    # 36: unc_len + usc_len
    # so, if stop < unc_len + usc_len, then all noncoding
    # The last base is 36 for no change, so we also allow for =

    if stop < unc_len + usc_len:
        # All events occur upstream or upstream + partial
        # No coding changes
        #return '', '', ref, alt, -1
        return '', '', '', '', -1
    if start >= unc_len + usc_len + seq_len:
        # All events are downstream
        # No coding changes
        #return '', '', ref, alt, -1
        return '', '', '', '', -1
    if start > us_len:
        # Starts in coding region. May/May not have downstream noncoding
        if stop > seq_len:
            # Trim the downstream bases
            ref = ref[:-(stop - seq_len) + 1]
            alt = alt[:-(stop - seq_len) + 1]
            new_start = start
    else:
        if start >= unc_len + usc_len and stop <= unc_len + usc_len + seq_len:
            # Starts and ends in coding, so no changes
            new_start = start
        elif start < unc_len and stop < us_len:
            # Starts upstream, ends in coding region
            # Can't use this
            return '', '', '', '', -1
        else:
            # Has downstream noncoding to trim off
            ref = ref[: -(stop - (unc_len + usc_len + seq_len))]
            alt = alt[: -(stop - (unc_len + usc_len + seq_len))]
            new_start = start
            downstream_trim = (stop - (unc_len + usc_len + seq_len))


    codon_index = (new_start - unc_len - usc_len) % 3

    refDNA, altDNA = get_codon(ref, alt, seq, codon_index, new_start)

    # Define the reference amino acid sequence
    if strand == 'F':
        refAA = Seq(refDNA).transcribe().translate(to_stop=False).__str__()
        altAA = Seq(altDNA).transcribe().translate(to_stop=False).__str__()
        AApos = (new_start - unc_len - usc_len) // 3
    else:
        refAA = Seq(refDNA).reverse_complement().transcribe().translate(to_stop=False).__str__()
        altAA = Seq(altDNA).reverse_complement().transcribe().translate(to_stop=False).__str__()
        AApos = (us_len - new_start - 1) // 3   # Should be 0-based
    return refAA, altAA, refDNA, altDNA, AApos


def get_codon(refDNA, altDNA, seq, codon_index, start):
    """
    Get the correct codon triplets being used.

    First, start by adding nucleotides to the left side, if necessary
    Second, add additional reads to the right until you have a fragment divisible by 3

    :param refDNA: ref allele, trimmed to only include coding bases
    :param altDNA: alt allele, trimmed to only include coding bases
    :param seq: reference sequence of the target region that is coding
    :param codon_index: 0, 1, or 2
    :param start: where the variant starts, relative to seq

    :return: refDNA-codons, altDNA-codons
    """
    # Add bases to the left side
    new_start = start - codon_index
    refDNA = seq[new_start:start] + refDNA
    altDNA = seq[new_start:start] + altDNA

    # Add bases to the right side until divisible by 3
    while refDNA.__len__() % 3 > 0:
        try:
            newBase = seq[new_start + refDNA.__len__()]
            altDNA += newBase
            refDNA += newBase
        except IndexError:
            refDNA = '?'
            altDNA = '?'

    while '?' in refDNA or refDNA.__len__() % 3 > 0:
        # When you go too far into the read you can run out of sequence.
        # In this case, we just drop back to the previous codon
        refDNA = refDNA[:-1]
        altDNA = altDNA[:-1]
    return refDNA, altDNA


def position_wrapper(read, target_dict, read_start, other_dict, strand):
    """
    Performs all of the functions in this document

    :param read_start:
    :param other_dict:
    :param strand:
    :param read: full read from SAM file
    :param target_dict: must contain:
            `seq`: the target sequence from reference genome
            `mmc`: list of DNA mismatch index positions that one would expect to observe
            `codon_distance`: how many bases from start of the read to the first full codon position

    :return:
    {
        `pos_ref_alt`:  'start_pos_idx'_'stop_pos_idx'_'DNA_ref'_'DNA_alt',
        `RefAA`: what the reference Amino acid(s) should be
        `AltAA`: what the alternate Amino acid(s) are
        `cons`: A boolean as to whether or not all the changes were consecutive
    }
    """
    mmc = get_dna_changes(target_dict['seq'], read)
    mmc = remove_expected_changes(mmc, target_dict['mmc'], 0)
    cons = is_consecutive(mmc)
    pos_ref_alt = trim_contig(mmc, target_dict['seq'], read)
    RefAA, AltAA, refDNA, altDNA, AApos = get_affected_amino_acid(pos_ref_alt, target_dict['seq'],
                                                                  target_dict['from'] - read_start, other_dict, strand)
    return dict(pos_ref_alt=pos_ref_alt, RefAA=RefAA, AltAA=AltAA, cons=cons, AApos=AApos, refDNA=refDNA, altDNA=altDNA)
