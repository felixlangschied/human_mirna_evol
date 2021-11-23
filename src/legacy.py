

def read_matmir_gff(path):
    """
    TODO: correct index calculation before using this function
    :param path:
    :return:
    """
    out_dict = {}
    with open(path, 'r') as fh:
        for line in fh:
            if not line.startswith('#'):
                linedata = line.strip().split('\t')
                if linedata[2] == 'pre_miRNA':
                    mirna = linedata[-1].split(';')[0].replace('ID=', '').replace('_pre', '')
                    pre_start, pre_end = map(int, linedata[3:5])
                    matline = next(fh).split(';')[0].strip()
                    if matline.endswith('*'):
                        matline = next(fh) .split(';')[0].strip()
                    mat_start, mat_end = map(int, matline.split('\t')[3:5])

                    start_index = mat_start - pre_start
                    end_index = pre_end - pre_start
                    print(start_index, end_index)
                    out_dict[mirna] = (start_index, end_index)
    return out_dict


def find_mature(preseq, mature, trim=2):
    seq_degap = preseq.replace('-', '')
    # find position of gaps in alignment
    gap_pos = [index for index in range(len(preseq)) if preseq.startswith('-', index)]
    # use partial mature sequence (some ncortho hits are 1-2 nucleotides truncated)
    trim_matseq = mature[trim:-trim]

    # search for mature sequence
    mat_found = [m for m in re.finditer(re.escape(trim_matseq), seq_degap)]
    if len(mat_found) == 0:
        print('Could not find mature sequence in the pre-seq')
        return False, False
    elif len(mat_found) > 1:
        print('Found multiple hits of mature sequence in pre-seq. Trimming was too unspecific')
        return False, False

    match_start, match_stop = mat_found[0].span()
    if match_stop - match_start < 0.7 * len(mature):
        print('Hit too short')
        return False, False

    # increase index by number of gaps before mature miRNA
    inserted = [gap for gap in gap_pos if gap <= match_start]
    num_inserted = len(inserted)
    match_start += num_inserted
    match_stop += num_inserted

    # untrim mature sequence
    match_start -= trim
    match_stop += trim

    # increase end of mature miRNA index by number of gaps inside the mature miRNA
    matchrange = range(match_start, match_stop)
    mat_gaps = [gap for gap in gap_pos if gap in matchrange]
    match_stop += len(mat_gaps)


    return match_start, match_stop