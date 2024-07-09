import re

line = ['M00738:252:000000000-JF4NM:1:1101:15954:14434', '2048', 'chr17', '56772604', '60', '120H31M', '*', '0', '0',
        'TAGTATTTTGTACTATCGTCAGGAAACCAAA', ']]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]', 'NM:i:0', 'MD:Z:31', 'AS:i:31', 'XS:i:19',
        'SA:Z:chr17,56772411,+,120M31S,60,2;']


def split_cigar(cigar):
    keys = list(filter(None, re.split('[0-9]+', cigar)))
    values = list(filter(None, re.split('\D+', cigar)))
    return dict(zip(keys, values))


def trim_soft(keys, values):
    if keys[0] == 'M':
        return 0, int(values[0])
    elif keys[0] == 'S' and keys[1] == 'M':
        return 0, int(values[1])
    elif keys[0] == 'H' and keys[1] == 'M':
        return 0, int(values[1])
    else:
        return None, None


def filter_matches(cigar_dict, min_matches=150):
    if 'M' in cigar_dict.keys():
        if int(cigar_dict['M']) >= min_matches:
            return True
    return False
