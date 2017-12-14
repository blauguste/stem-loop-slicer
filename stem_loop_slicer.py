from Bio import AlignIO
from Bio.AlignIO import StockholmIO
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import generic_rna
from itertools import takewhile
import numpy as np
import sys

def stem_loop_sep(infile):
    infile_bn = infile.split('_')[0]
    loop_out = infile_bn + '_loop.fas'
    stem_out = infile_bn + '_stem.fas'
    with open(infile, 'r') as stk_in:
        alignment = AlignIO.read(stk_in, 'stockholm')
        SS_cons = alignment.annotations['GC']['SS_cons'][0]
        print(alignment)
        print(alignment.annotations['GC']['SS_cons'][0])
        # Loop through the secondary structure to store loop and stem positions, respectively
        loop_pos = ([pos for pos, char in enumerate(SS_cons) if char == '.'])
        stem_pos = ([pos for pos, char in enumerate(SS_cons) if char == '('])
        stem_pos_b = ([pos for pos, char in enumerate(SS_cons) if char == ')'])
        full_stem_pos = stem_pos + stem_pos_b
        loop_records = []
        stem_records = []
        for rec in alignment:
            loop_only = Seq.Seq(''.join(list(rec.seq[i] for i in loop_pos)))
            new_loop_rec = SeqRecord.SeqRecord(loop_only, id=rec.id + '_loop')
            loop_records.append(new_loop_rec)
            stem_only = Seq.Seq(''.join(list(rec.seq[i] for i in full_stem_pos)))
            new_stem_rec = SeqRecord.SeqRecord(stem_only, id=rec.id + '_stem')
            print(len(new_loop_rec)+ len(new_stem_rec))
            stem_records.append(new_stem_rec)
        with open(loop_out, 'w') as loop_out:
            SeqIO.write(loop_records, loop_out, 'fasta')
        with open(stem_out, 'w') as stem_out:
            SeqIO.write(stem_records, stem_out, 'fasta')

if __name__ == '__main__':
    if len(sys.argv) == 2:
         stem_loop_sep(sys.argv[1])
    else:
         print("Usage: stem_loop_slicer.py seqs_in.stk")
         sys.exit(0)

        #alignment_loops_only = alignment[:, :9]
        #print(alignment_loops_only)
        #align_array = np.array([list(rec) for rec in alignment], np.character)
        #print(align_array)
