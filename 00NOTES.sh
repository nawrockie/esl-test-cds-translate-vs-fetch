# EPN, Fri Jun 19 15:02:29 2015
# pwd: /panfs/pan1/dnaorg/programs/15_0331_esl_test_cds_translate_vs_fetch
# 
# # NOTE: the name of this directory references the last day of March but this
# # file wasn't written until June 19, which was when this script was added
# # to the /panfs/pan1/dnaorg/programs/ directory and came under git 
# # revision control. 
# 
# Example of running esl-test-cds-translate-vs-fetch.pl to test a translated
# CDS sequence against the fetched protein sequence that is linked to it. 
#
# This is a Perl script that uses the Bio-Easel library which uses
# Inline to interact and call functions from the Easel sequence
# analysis library from Sean Eddy's group.
#
# git repository: https://github.com/nawrockie/esl-test-cds-translate-vs-fetch.git
#
#######################
# More information
#######################
#
# See /home/nawrocke/notebook/15_0331_esl_test_cds_translate_vs_fetch
# for notes on development and testing of this script.
# 
#######################
# Prerequisites
#######################
# 
# Directories that include the BioEasel perl modules must be part of your
#  $PERL5LIB environment variable.
# 
# To modify your PERL5LIB environment variable appropriately:
# For bash shell users
source /panfs/pan1/dnaorg/programs/setup-bio-easel.bash.sh
# For C shell or C shell compatible users
source /panfs/pan1/dnaorg/programs/setup-bio-easel.csh.sh

#######################
# Usage and options
#######################
#esl-test-cds-against-aa.pl [OPTIONS] <input fasta file output from esl-fetch-cds.pl>
#	OPTIONS:
#		-v         : be verbose; output translated and fetched protein sequences
#		-p         : print all sequences, even those that pass all tests
#		-subset    : only perform first six tests, not all nine
#		-incompare : input file was created by dnaorg_compare_genomes.pl -protid -codonstart
#		-skipinc   : skip examination of incomplete CDS'
#
############################
# Example command and output
############################
# 
perl esl-test-cds-translate-vs-fetch.pl sample.cds.fa 
#
# Output:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Complete CDS:
#protein-accession      incomplete?  start     T1   T2   T3   T4   T5   T6   T7   T8    pass/fail
AAH45158.1                       no      1      0    0    0    0    0    1    0    0    fail      # position 241 mismatch A ne P (translated ne fetched)
#
#protein-accession      incomplete?  start     T1   T2   T3   T4   T5   T6   T7   T8    pass/fail
# num-fails-complete             no    N/A      0    0    0    0    0    1    0    0    N/A
# num-fails-incomplete          yes    N/A      0    0    0    0    0    0    0    0    N/A
# num-fails-all                 N/A    N/A      0    0    0    0    0    1    0    0    N/A
#
# Summary:
#
# category    num-pass  num-fail  fract-fail
# complete           1         1      0.5000
# incomplete         2         0      0.0000
# all                3         1      0.2500
#
# Test 1 (T1): Length of CDS (DNA) is a multiple of 3 (or stop is incomplete)
# Test 2 (T2): Fetched protein starts with an M (or is annot. as incomplete on 5' end)
# Test 3 (T3): CDS ends with a stop codon (or is annot. as incomplete on 3' end)
# Test 4 (T4): CDS has 0 Ns that make AA assignment ambiguous
# Test 5 (T5): Translated CDS and fetched protein are identical length
#              (ignoring final codon, if it's a stop or incomplete)
# Test 6 (T6): Translated CDS and fetched protein are identical sequence
#              (ignoring hangovers if lengths differ and final codon (if stop) and codons with ambiguous nts)
# Test 7 (T7): No internal stops in translated CDS
# Test 8 (T8): No internal stops in fetched protein
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The output explains the tests run for each CDS sequence and which
# failed. Each line that does not begin with a '#' lists information
# for a different sequence. By default, only sequences which have at
# least one failure are printed, use -p to print all sequences, even
# those that pass all tests.
# 
###################################################
# Example command and output with -incompare option
###################################################
# 
# This script was extended to work with sequence files output from
# dnaorg_compare_genomes.pl (see notes in
# /panfs/pan1/dnaorg/programs/15_0528_dnaorg_compare_genomes) when
# using the options -protid and -codonstart with
# dnaorg_compare_genomes.pl.  Those fasta files have a special naming
# convention format that requires special processing to check the CDS.
# 
# Here is an example of running the script on those files:
#
perl esl-test-cds-translate-vs-fetch.pl -incompare sample.cds.compare.in
#
# Output:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Complete CDS:
#protein-accession      incomplete?  start     T1   T2   T3   T4   T5   T6   T7   T8    pass/fail
AIJ19251.1                       no      1      0    0    1    0    1    0    0    0    fail     
#
#protein-accession      incomplete?  start     T1   T2   T3   T4   T5   T6   T7   T8    pass/fail
# num-fails-complete             no    N/A      0    0    1    0    1    0    0    0    N/A
# num-fails-incomplete          yes    N/A      0    0    0    0    0    0    0    0    N/A
# num-fails-all                 N/A    N/A      0    0    1    0    1    0    0    0    N/A
#
# Summary:
#
# category    num-pass  num-fail  fract-fail
# complete         176         1      0.0056
# incomplete         2         0      0.0000
# all              178         1      0.0056
#
# Test 1 (T1): Length of CDS (DNA) is a multiple of 3 (or stop is incomplete)
# Test 2 (T2): Fetched protein starts with an M (or is annot. as incomplete on 5' end)
# Test 3 (T3): CDS ends with a stop codon (or is annot. as incomplete on 3' end)
# Test 4 (T4): CDS has 0 Ns that make AA assignment ambiguous
# Test 5 (T5): Translated CDS and fetched protein are identical length
#              (ignoring final codon, if it's a stop or incomplete)
# Test 6 (T6): Translated CDS and fetched protein are identical sequence
#              (ignoring hangovers if lengths differ and final codon (if stop) and codons with ambiguous nts)
# Test 7 (T7): No internal stops in translated CDS
# Test 8 (T8): No internal stops in fetched protein
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#############
# Input files
#############
# 
# This script takes as input fasta files created by esl-fetch-cds.pl.
#
# Expected input:
# If esl-fetch-cds.pl was run in default mode with a particular format
# of input file. Those sequences will be
# named something like the ones in the sample input file in this
# dir 'sample.cds.in'. For example,
#
# >AAH45158.1:codon_start1:BC045158.1:4:870:+:
# 
# Tokens are delimited by ':'. The first is simply from the input
# file, as explained below. The second indicates the position of the
# first nucleotide to be translated as the single digit after the
# string 'codon_start'. Usually this is 1, meaning the first nt is
# translated. If it is 2, it means the first nucleotide is *not*
# translated and the second nucleotide is the first one to be
# translated. The only other valid value is 3. Tokens 3 to N come
# in sets of 4. The first in each set gives a nucleotide accession
# that encodes a segment (i.e. exon) of the CDS, followed by the
# start position, then the stop, then the strand. So in this 
# example the protein AAH45158.1 is encoded by BC045158.1 from
# positions 4 to 870 on the positive strand.
#
#########
#
# Special input from dnaorg_compare_genomes.pl:
#
# In June, 2015 I extended the script to work with a special
# type of fasta file output from dnaorg_compare_genomes.pl.
#
# If you run dnaorg_compare_genomes.pl with cmd line options
# '-s -protid -codonstart', you will get fasta files with
# sequences named like the ones in the sample input
# file 'sample.cds.compare.fa'. An example is:
# >L13994 L13994:1376:1840:+: product:preS protein_id:gb|AAC31754.1| codon_start:-
# 
# For this type of input fasta file, the description (not the
# sequence name) contains the relevant information.
# In this sequence name, the first white-space delimited token
# after the name ('L13994:1376:1840:+:') is in the same format
# as the 3rd to N sets of 4 tokens in the default case. 
# First ':' delimited token is a nucleotide accession encoding
# a segment (e.g. exon), second token is start posn, third
# token is stop posn, fourth token is strand. There can be more
# than one set of 4 tokens here, but there will always be a 
# multiple of 4 total ':'-delimited tokens.
# 
# The protein accession is in the token beginning with 'protein_id:'.
# In this case it is 'AAC31754.1'.
# 
# The codon_start is likewise in the token beginning with 'codon_start:'.
# If it is not '1', '2', or '3', it will be treated as if it is '1'
# as this is the default value for codon_start in NCBI annotation.
#
##############################################
# Last updated: EPN, Tue Jun 23 14:20:43 2015
##############################################
