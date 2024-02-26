#!/usr/bin/env python

##A set of tools, that can be used to convert sequences to all different sequencing
##format types. Also convert DNA to RNA and both to amino acids

##import libraries used in script
import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
#import subprocess
#import glob
import pandas as pd
import sys

##set up arguments
parser = argparse.ArgumentParser(description=(
	'Sequence conversion script. Specify formats/types in lower case. Formats accepted: fasta, clustal, emboss, nexus, phylip, phylip-sequential and phylip-relaxed'
						)
				)
parser.add_argument('-d', '--working_directory', dest='directory', help='directory input file(s) are kept in', required=True)
parser.add_argument('-i', '--input', dest='in_filename', help='name of input file, including extension', required=True)
parser.add_argument('-t', '--input_format', dest='in_form', help= "alignment format of input file ('fasta', 'phylip', 'nexus')", required =True)
parser.add_argument('-ot', '--output_format', dest='out_form', help="alignment format of input file ('fasta', 'phylip', 'nexus')", required=True)
parser.add_argument('-od', '--output_directory', dest='out_directory', help='specify directory for output file if different')
parser.add_argument('-1', '--op', dest='op', help='specify operation you wish to be carried out (seqeunce conversion = seq_con; sequence selection = seq_sel; list sequences = seq_lst; change to RNA = seq_tranl; change to prot = seq_tranc; sequence lengths = seq_len (if you do not give a out_filename it will print to screen...)', required=True)
parser.add_argument('-s', '--isolate_id', dest='id', help= 'id of fasta sequence you wish to extract', required = 'seq_sel' in sys.argv and '-f' not in sys.argv) #conditional argument based on prior arguments
parser.add_argument('-f', '--isolate_id_list', dest='id_file', help= "single column text file with sequences ids you wish to select. Ensure header is 'ids'", required = 'seq_sel' in sys.argv and '-s' not in sys.argv) #conditional argument based on prior arguments
parser.add_argument('-o', '--output', dest='out_filename', help='name of output file, including extension', required = 'seq_con' in sys.argv or 'seq_tranc' in sys.argv or 'seq_tranl' in sys.argv) #conditional argument based on prior arguments
parser.add_argument('-st', '--sequence_type', dest='seq_type', help='sequence type of the input file: dna, rna', required = 'seq_tranc' in sys.argv) #conditional argument based on prior arguments


args = parser.parse_args()

##set working directory
os.chdir(args.directory)

##opperation selection
if args.op == 'seq_con':

	##load fasta sequence (AlignIO.parse returns a MSA as an interator, to get
	#access to seperate sequences it has to be saved into a list)
	with open(args.out_filename, 'w') as output_handle:
		#convert and write out sequence file
		alignment = AlignIO.read(input_handle, args.in_form)
		seq_out = AlignIO.write(alignment, output_handle, args.out_form)
		#some alignment stats output automatically
		print('Alignment length is %i' % alignment.get_alignment_length())
		print('%i sequences in alignment' % len(alignment))

	#if args.fastq:

elif args.op == 'seq_sel':

	#load optional arguments for thispartivular operation
	if args.id:
		wanted = str(args.id)
		fasta_dict = SeqIO.index(args.in_filename, args.in_form) #Create a dictionary out of the fasta so that is can be indexed for a particular sequence

		end = False
		with open(wanted+'.'+args.out_form, "w") as f:
			curr_seq = fasta_dict[wanted]
			SeqIO.write(curr_seq, f, args.out_form)

	if args.id_file:
		id_df = pd.read_table(args.id_file)
		wanted = id_df['ids']
		fasta_dict = SeqIO.index(args.in_filename, args.in_form) #Create a dictionary out of the fasta so that is can be indexed for a particular sequence

		end = False
		for seq in wanted:
			with open(seq+'.'+args.out_form, "w") as f:
				curr_seq = fasta_dict[seq]
				SeqIO.write(curr_seq, f, args.out_form)

elif args.op == 'seq_lst':

	fasta_dict = SeqIO.index(args.in_filename, args.in_form)
	print(list(fasta_dict.keys()))

elif args.op == 'seq_len':

	if args.out_filename:
		with open(args.out_filename, 'w') as nt_fa:
			#access seperate sequences with for loop)
			for index, dna_record in enumerate(SeqIO.parse(args.in_filename, args.in_form)):
				nt_fa.write(
					"ID = %s, length %i\n"
					% (dna_record.id, len(dna_record.seq))
					)

	else:
		#access seperate sequences with for loop)
		for index, dna_record in enumerate(SeqIO.parse(args.in_filename, args.in_form)):
			print(
				"ID = %s, length %i\n"
				% (dna_record.id, len(dna_record.seq))
				)

elif args.op == 'seq_tranl':

	##load fasta sequence
	with open(args.out_filename, 'w') as aa_fa:
		#access seperate sequences with for loop)
		for dna_record in SeqIO.parse(args.in_filename, args.in_form):
			#variable of the particular sequence iterate
			dna_seqs = [dna_record.seq, dna_record.seq.reverse_complement()]
			#generate all translation frames
			aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in dna_seqs)
			# select the longest one
			max_aa = max(aa_seqs, key=len)
			# write new record
			aa_record = SeqRecord(max_aa, id=dna_record.id, description="translated sequence")
			SeqIO.write(aa_record, aa_fa, args.out_form)

elif args.op == 'seq_tranc':
	#conditional action if sequence type is DNA
	if args.seq_type == 'dna':
		##load fasta sequence
		with open(args.out_filename, 'w') as transc_fa:
			#access seperate sequences with for loop)
			for dna_record in SeqIO.parse(args.in_filename, args.in_form):
				#variable of the particular sequence iterate
				dna_seqs = dna_record.seq
				#generate all translation frames
				transc_seqs = dna_seqs.transcribe()
				# write new record
				transc_record = SeqRecord(transc_seqs, id=dna_record.id, description="transcribed sequence")
				SeqIO.write(transc_record, transc_fa, args.out_form)

	#conditional action if sequence type is RNA
	if args.seq_type == 'rna':
		##load fasta sequence
		with open(args.out_filename, 'w') as transc_fa:
			#access seperate sequences with for loop)
			for rna_record in SeqIO.parse(args.in_filename, args.in_form):
				#variable of the particular sequence iterate
				rna_seqs = rna_record.seq
				#generate all translation frames
				transc_seqs = rna_seqs.back_transcribe()
				# write new record
				transc_record = SeqRecord(transc_seqs, id=rna_record.id, description="transcribed sequence")
				SeqIO.write(transc_record, transc_fa, args.out_form)
