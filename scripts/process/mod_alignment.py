"""
Take alignments of vertebrate miRNA orthologs and identify the columns that align
with the mature miRNA in human (reference species) identified by their annotation in MirGeneDB.

Notes:
    * In the case of large insertions in the region of the mature miRNA (many gaps),
    we will delete a large number of columns during the jackknifing

"""
from Bio import AlignIO
import glob
import os
import re
import random as rd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

################################################################################################
align_dir = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\raw\refseq_alignments'
gff_path = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\raw\hsa.gff'
mat_path = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\raw\hsa.fas'

nomat_outdir = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\processed\no_mat_alns'
jack_outdir = r'C:\Users\felix\PycharmProjects\human_mirna_evol\data\processed\jackknifed_alns'
################################################################################################


def read_matmir_fasta(path):
    out_dict = {}
    with open(path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                mirna = line.strip().replace('>', '').split('_')[0]
            else:
                matseq = line.strip().replace('U', 'T')
                out_dict[mirna] = matseq
    return out_dict


def make_mat_regex(seq):
    regex_list = ['-*']
    for nucl in seq:
        regex_list.append(nucl)
        regex_list.append('-*')
    regex_list.append('-*')
    regex = r''.join(regex_list)
    return regex


def regex_mature(preseq, mature):
    trim = 0
    mat_regex = make_mat_regex(mature)
    res = re.search(mat_regex, preseq)
    if res is None:
        trim = 1
        trim_mat_regex = make_mat_regex(mature[trim:-trim])
        res = re.search(trim_mat_regex, preseq)
        if res is None:
            trim = 2
            trim_mat_regex = make_mat_regex(mature[trim:-trim])
            res = re.search(trim_mat_regex, preseq)
            if res is None:
                print('No Hit found')
                return '', ''
    return res.group().strip('-'), trim


def get_mature_pos(alignment, matseq):
    for row in alignment:
        if row.name == 'Homo_sapiens':
            human_seq = str(row.seq)
            mat_hit, padding = regex_mature(human_seq, matseq)
            if mat_hit:
                mat_start = human_seq.index(mat_hit) - padding
                mat_end = mat_start + len(mat_hit) + padding
                mat_length = mat_end - mat_start
                return mat_start, mat_end, mat_length
            else:
                return '', '', ''


def alignment_jackknife(alignment, matlength):
    rand_ind = np.sort(rd.sample(range(alignment.get_alignment_length()), matlength))

    aln_mat = np.array([list(rec) for rec in alignment], 'S1')
    sampled_aln_mat = np.delete(aln_mat, rand_ind, axis=1)

    sample_list = []
    for rowc, record in enumerate(alignment):
        rowseq = sampled_aln_mat[rowc, :].tobytes().decode('utf-8')
        row = SeqRecord(Seq(rowseq), id=record.id)
        sample_list.append(row)
    sample_aln = MultipleSeqAlignment(sample_list)
    return sample_aln


################################################################################################

mat_dict = read_matmir_fasta(mat_path)

aln_files = glob.glob(f'{align_dir}/*')
mirnas_skipped = 0
for file in aln_files:
    mirna = file.split(os.sep)[-1].replace('.aln', '')
    mature = mat_dict[mirna]
    full_aln = AlignIO.read(file, 'fasta')

    start, end, length = get_mature_pos(full_aln, mature)
    if not start:  # skip alignments where mature not found in ncOrtho hit (2 in total)
        continue
    no_mat_aln = full_aln[:, :start] + full_aln[:, end:]
    nomat_out = f'{nomat_outdir}/{mirna}.aln'
    AlignIO.write(no_mat_aln, nomat_out, "fasta")
    print('# Finished writing alignments of miRNAs without the mature sequence')

    jack_aln = alignment_jackknife(full_aln, length)
    jack_out = f'{jack_outdir}/{mirna}.aln'
    AlignIO.write(jack_aln, jack_out, 'fasta')
    print('# Finished writing jackknifed alignments')
