# takes as input a list of fasta files and outputs a qsub file that will run esl-test-cds-translate-vs-fetch.pl for each
$usage  = "perl jiffy-list2qsub.pl\n";
$usage .= "\t<list of files to run esl-test-cds-translate-vs-fetch.pl on>\n";
$usage .= "\t<0 or more options to add to esl-test-cds-translate-vs-fetch.pl options>\n";

$edir = "/panfs/pan1/dnaorg/programs";

if(scalar(@ARGV) < 1) { die $usage; }

$listfile = shift @ARGV;
$options = "";
while(scalar(@ARGV) > 0) { 
  if($options ne "") { $options .= " "; }
  $options .= shift @ARGV;
}

$ctr = 0;
open(IN, $listfile) || die "ERROR unable to open $listfile for reading";
while($file = <IN>) { 
  chomp $file;
  $root = $file; 
  if($root !~ m/\.fa$/) { die "ERROR file $file does not end in .fa"; }
  # NOTE THAT WE DO NOT REMOVE THE DIR PATH, SO ERR AND OUT FILES WILL BE PLACED IN SAME DIR AS FASTA FILE
  $root =~ s/\.fa$//;
  $jobname = $root;
  $errfile = $root . ".cds-test.err";
  $outfile = $root . ".cds-test";
  printf("qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n \"perl $edir/esl-test-cds-translate-vs-fetch.pl $options $file > $file.name_accn.report\"\n");
}
