# seq_tools

Set of commands to help manage/manipulate multiple/aligned sequence files. 

The current list of operations includes:

1: List all sequence in input file

2: Select sequences to output into separate fasta files

3: Convert an aligned file of any format listed to any other format listed

4: Translate sequences in input file 

5: Transcribe sequences in input file

6: Get length of all sequences in input file

------------------------------------------------------------------------------------------------

Example command:
```
seq_tools.py -d ./ -i seq_test.fasta -t fasta -1 seq_tranc -o tranc_test.fasta -ot fasta -st dna
```

Script help file:
```
-h, --help            show this help message and exit
-d DIRECTORY, --working_directory DIRECTORY
                        directory input file(s) are kept in
-i IN_FILENAME, --input IN_FILENAME
                        name of input file, including extension
-t IN_FORM, --input_format IN_FORM
                        alignment format of input file ('fasta', 'phylip', 'nexus')
-ot OUT_FORM, --output_format OUT_FORM
                        alignment format of input file ('fasta', 'phylip', 'nexus')
-1 OP, --op OP        specify operation you wish to be carried out (seqeunce conversion = seq_con; sequence selection =
                        seq_sel; list sequences = seq_lst; change to amino acids (aa) = seq_tranl; change to mRNA =
                        seq_tranc; sequence lengths = seq_len (if you do not give a out_filename it will print to screen...)
-s ID, --isolate_id ID
                        id of fasta sequence you wish to extract
-f ID_FILE, --isolate_id_list ID_FILE
                        single column text file with sequences ids you wish to select. Ensure header is titled 'ids'
-o OUT_FILENAME, --output OUT_FILENAME
                        name of output file, including extension
-st SEQ_TYPE, --sequence_type SEQ_TYPE
                        sequence type of the input file: dna, rna
```
