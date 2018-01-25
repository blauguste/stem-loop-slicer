from Bio import AlignIO
from Bio.AlignIO import StockholmIO
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import generic_rna
import sys

def stem_loop_sep(infile):

    infile_bn = infile.split('_')[0]
    loop_out = infile_bn + '_loop.fas'
    stem_out = infile_bn + '_stem.fas'
    trim_out = infile_bn + '_trim.fas'

    with open(infile, 'r') as stk_in:
        alignment = AlignIO.read(stk_in, 'stockholm')
        SS_cons = alignment.annotations['GC']['SS_cons'][0]
        
        # Remove gapped sites
        lead_trim = 0
        trail_trim = 0

        for i in range(alignment.get_alignment_length()):
                if '-' in alignment[:, i]:
                    lead_trim = i + 1
                else: break

        for i in range(1, alignment.get_alignment_length()):
                if '-' in alignment[:, -i]:
                    trail_trim = -i
                else: break

        if trail_trim != 0:
            trim_align = alignment[:, lead_trim:trail_trim]
            trim_struc = SS_cons[lead_trim:trail_trim]
        else:
            trim_align = alignment[:, lead_trim:]
            trim_struc = SS_cons[lead_trim:]

        gap_ct = 0

        for i in range(trim_align.get_alignment_length()-1):
            if '-' in trim_align[:, i]:
                gap_ct += 1

        print('Original alignment length:', alignment.get_alignment_length())
        print('Leading bases trimmed:', lead_trim)
        print('Trailing bases trimmed:', trail_trim*-1)
        print('New alignment length:', trim_align.get_alignment_length())
        print('Number of gapped sites:', gap_ct)

        # Loop through the secondary structure to store loop and stem positions, respectively
        loop_pos = ([pos for pos, char in enumerate(trim_struc) if char == '.'])
        stem_pos = ([pos for pos, char in enumerate(trim_struc) if char == '(' or char == ')'])
        loop_records = []
        stem_records = []

        for rec in trim_align:
            loop_only = Seq.Seq(''.join(list(rec.seq[i] for i in loop_pos)))
            new_loop_rec = SeqRecord.SeqRecord(loop_only, id=rec.id + '_loop')
            loop_records.append(new_loop_rec)
            stem_only = Seq.Seq(''.join(list(rec.seq[i] for i in stem_pos)))
            new_stem_rec = SeqRecord.SeqRecord(stem_only, id=rec.id + '_stem')
            stem_records.append(new_stem_rec)

        with open(loop_out, 'w') as loop_out:
            SeqIO.write(loop_records, loop_out, 'fasta')
        with open(stem_out, 'w') as stem_out:
            SeqIO.write(stem_records, stem_out, 'fasta')
        with open(trim_out, 'w') as trim_out:
            AlignIO.write(trim_align, trim_out, 'fasta')

if __name__ == '__main__':
    if len(sys.argv) == 2:
         stem_loop_sep(sys.argv[1])
    else:
         print("Usage: stem_loop_slicer.py seqs_in.stk")
         sys.exit(0)
