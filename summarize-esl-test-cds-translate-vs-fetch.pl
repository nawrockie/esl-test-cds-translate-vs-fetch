#!/usr/bin/env perl
# 
# summarize-esl-test-cds-translate-vs-fetch.pl: 
#
#  Given a file that lists multiple output files from 
#  esl-test-cds-translate-vs-fetch.pl, summarize their output.
#
# EPN, Wed Jun 24 08:55:07 2015

use strict;
use Getopt::Long;

my $in_listfile = "";    # name of input file to 
my $outroot     = undef; # changed with -oroot <s>

&GetOptions( "oroot=s" => \$outroot);

my $usage;
$usage  = "summarize-esl-test-cds-translate-vs-fetch.pl [OPTIONS] <file with list of multiple esl-test-cds-against input fasta file output from esl-fetch-cds.pl>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-oroot <s>: output sequences that fail each test to <s>.T<n>.fail-list and seqs that fail any test to <s>.any.fail-list\n";

if(scalar(@ARGV) != 1) { die $usage; }
($in_listfile) = @ARGV;

if(! -s $in_listfile) { die "ERROR, $in_listfile does not exist"; }

my @extra_lines_A = ();

my %pass_all_ct_H = ();
my %fail_all_ct_H = ();
my %num_Ns_all_ct_H = ();
my %num_oth_all_ct_H = ();
my %t1_all_ct_H = ();
my %t2_all_ct_H = ();
my %t3_all_ct_H = ();
my %t4_all_ct_H = ();
my %t5_all_ct_H = ();
my %t6_all_ct_H = ();
my %t7_all_ct_H = ();
my %t8_all_ct_H = ();

my %pass_incomp_ct_H = ();
my %fail_incomp_ct_H = ();
my %num_Ns_incomp_ct_H = ();
my %num_oth_incomp_ct_H = ();
my %t1_incomp_ct_H = ();
my %t2_incomp_ct_H = ();
my %t3_incomp_ct_H = ();
my %t4_incomp_ct_H = ();
my %t5_incomp_ct_H = ();
my %t6_incomp_ct_H = ();
my %t7_incomp_ct_H = ();
my %t8_incomp_ct_H = ();

# output lines for each test failure or ambig chars
my @num_Ns_A = ();
my @num_oth_A = ();
my @t1_A = ();
my @t2_A = ();
my @t3_A = ();
my @t4_A = ();
my @t5_A = ();
my @t6_A = ();
my @t7_A = ();
my @t8_A = ();

my @file_A = (); # list of file names, as read from $in_listfile
# open list file and parse each file listed in it

# open output files, if nec:
my $N_list_file   = undef;
my $ab_list_file  = undef;
my $t1_fail_file  = undef;
my $t2_fail_file  = undef;
my $t3_fail_file  = undef;
my $t4_fail_file  = undef;
my $t5_fail_file  = undef;
my $t6_fail_file  = undef;
my $t7_fail_file  = undef;
my $t8_fail_file  = undef;
my $any_fail_file = undef;

if(defined $outroot) { 
  $N_list_file   = $outroot . ".Ns.list";
  $ab_list_file  = $outroot . ".ab.list";
  $t1_fail_file  = $outroot . ".T1.fail-list";
  $t2_fail_file  = $outroot . ".T2.fail-list";
  $t3_fail_file  = $outroot . ".T3.fail-list";
  $t4_fail_file  = $outroot . ".T4.fail-list";
  $t5_fail_file  = $outroot . ".T5.fail-list";
  $t6_fail_file  = $outroot . ".T6.fail-list";
  $t7_fail_file  = $outroot . ".T7.fail-list";
  $t8_fail_file  = $outroot . ".T8.fail-list";
  $any_fail_file = $outroot . ".any.fail-list";
  open(OUTN,  ">" . $N_list_file)  || die "ERROR unable to open $N_list_file for writing"; 
  open(OUTAB, ">" . $ab_list_file) || die "ERROR unable to open $ab_list_file for writing"; 
  open(OUTT1, ">" . $t1_fail_file) || die "ERROR unable to open $t1_fail_file for writing"; 
  open(OUTT2, ">" . $t2_fail_file) || die "ERROR unable to open $t2_fail_file for writing"; 
  open(OUTT3, ">" . $t3_fail_file) || die "ERROR unable to open $t3_fail_file for writing"; 
  open(OUTT4, ">" . $t4_fail_file) || die "ERROR unable to open $t4_fail_file for writing"; 
  open(OUTT5, ">" . $t5_fail_file) || die "ERROR unable to open $t5_fail_file for writing"; 
  open(OUTT6, ">" . $t6_fail_file) || die "ERROR unable to open $t6_fail_file for writing"; 
  open(OUTT7, ">" . $t7_fail_file) || die "ERROR unable to open $t7_fail_file for writing"; 
  open(OUTT8, ">" . $t8_fail_file) || die "ERROR unable to open $t8_fail_file for writing"; 
}

my $header_line1 = sprintf("%-50s  %9s  %9s  %9s  %9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s\n", 
                           "#", "num-pass", "num-fail", "num-w-Ns", "num-w-oth", "    T1", "    T2", "    T3", "    T4", "    T5", "    T6", "    T7", "    T8");
my $header_line2 = sprintf("%-50s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s  %4s %4s\n",
                           "# filename", "tot", "inc", "tot", "inc", "tot", "inc", "tot", "inc", 
                           "tot", "inc", "tot", "inc", "tot", "inc", "tot", "inc", "tot", "inc", "tot", "inc", "tot", "inc", "tot", "inc");

print $header_line1; 
print $header_line2; 

open(LIST, $in_listfile) || die "ERROR unable to open $in_listfile for reading";
my $file_ctr = 0;
while(my $line = <LIST>) { 
  chomp $line;
  $file_ctr++;
  my $file = $line;
  if(! -s $file) { die "ERROR, unable to open $file listed in $in_listfile"; }
  # Lines we parse (except first one):
  #protein-accession      incomplete?  start  tr-len   num-Ns  num-oth     T1   T2   T3   T4   T5   T6   T7   T8    pass/fail
  # AAO63223.1                       no      1    1509        0        2      0    0    0    0    0    0    0    0    pass     
  # num-fails-complete             no      -       -        -        -      0    0    0    0    0    0    0    0    N/A
  # num-fails-incomplete          yes      -       -        -        -      0    0    0    0    0    0    0    0    N/A
  # num-fails-all                   -      -       -        -        -      0    0    0    0    0    0    0    0    N/A

  # category    num-pass  num-fail  fract-fail
  # complete         187         0      0.0000
  # incomplete        20         0      0.0000
  # all              207         0      0.0000
  open(IN, $file) || die "ERROR unable to open $file for reading"; 
  push(@file_A, $file);
  $num_Ns_all_ct_H{$file}     = 0;
  $num_oth_all_ct_H{$file}    = 0;
  $num_Ns_incomp_ct_H{$file}  = 0;
  $num_oth_incomp_ct_H{$file} = 0;
  while($line = <IN>) { 
    if($line !~ /^\#/) { 
      #AAO63223.1                       no      1    1509        0        2      0    0    0    0    0    0    0    0    pass     
      if($line =~ /^(\S+)\s+(\S+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\S\s+(\S+)\s+\d\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/)
      { 
        my ($paccn, $ntaccn, $incomplete, $num_Ns, $num_oth, $t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);
        my $ntaccn_str  = $ntaccn;
        $ntaccn_str  =~ s/\,/\\n/g;
        $ntaccn_str .= "\n";
        if($num_Ns > 0) { 
          $num_Ns_all_ct_H{$file} += $num_Ns;
          if(defined $N_list_file) { print OUTN $ntaccn_str; }
          if($incomplete =~ m/yes/) { $num_Ns_incomp_ct_H{$file} += $num_Ns; }
        }
        if($num_oth > 0) { 
          $num_oth_all_ct_H{$file} += $num_oth;
          if(defined $ab_list_file) { print OUTAB $ntaccn_str; }
          if($incomplete =~ m/yes/) { $num_oth_incomp_ct_H{$file} += $num_oth; }
        }
        if($num_Ns  > 0) { push(@num_Ns_A, $line); }
        if($num_oth > 0) { push(@num_oth_A, $line); }
        if($t1 != 0) { push(@t1_A, $line); }
        if($t2 != 0) { push(@t2_A, $line); }
        if($t3 != 0) { push(@t3_A, $line); }
        if($t4 != 0) { push(@t4_A, $line); }
        if($t5 != 0) { push(@t5_A, $line); }
        if($t6 != 0) { push(@t6_A, $line); }
        if($t7 != 0) { push(@t7_A, $line); }
        if($t8 != 0) { push(@t8_A, $line); }
      }
      else { 
        die "ERROR unable to parse line in $file: $line";
      }
    } # end of 'if($line !~ /^\#/)'
    if($line =~ /^\# num-fails-incomplete\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) { 
      # num-fails-incomplete          yes      -       -        -        -      0    0    0    0    0    0    0    0    N/A
      $t1_incomp_ct_H{$file} = $1;
      $t2_incomp_ct_H{$file} = $2;
      $t3_incomp_ct_H{$file} = $3;
      $t4_incomp_ct_H{$file} = $4;
      $t5_incomp_ct_H{$file} = $5;
      $t6_incomp_ct_H{$file} = $6;
      $t7_incomp_ct_H{$file} = $7;
      $t8_incomp_ct_H{$file} = $8;
    }
    if($line =~ /^\# num-fails-all\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) { 
      # num-fails-all                   -      -       -        -        -      0    0    0    0    0    0    0    0    N/A
      $t1_all_ct_H{$file} = $1;
      $t2_all_ct_H{$file} = $2;
      $t3_all_ct_H{$file} = $3;
      $t4_all_ct_H{$file} = $4;
      $t5_all_ct_H{$file} = $5;
      $t6_all_ct_H{$file} = $6;
      $t7_all_ct_H{$file} = $7;
      $t8_all_ct_H{$file} = $8;
    }
    if($line =~ m/^\# incomplete\s+(\d+)\s+(\d+)/) { 
      # incomplete        20         0      0.0000
      $pass_incomp_ct_H{$file} = $1;
      $fail_incomp_ct_H{$file} = $2;
    }
    if($line =~ m/^\# all\s+(\d+)\s+(\d+)/) { 
      # all              207         0      0.0000
      $pass_all_ct_H{$file} = $1;
      $fail_all_ct_H{$file} = $2;
    }
    if($file_ctr == 1) { 
      if($line =~ m/^\# Test/) { push @extra_lines_A, $line; }
    }
  }
  close(IN);
  if(! exists $t1_all_ct_H{$file})      { die "ERROR didn't read num-fails-all line in $file"; }
  if(! exists $t1_incomp_ct_H{$file})   { die "ERROR didn't read num-incomplete-all line in $file"; }
  if(! exists $pass_incomp_ct_H{$file}) { die "ERROR didn't read incomplete line in $file"; }
  if(! exists $pass_all_ct_H{$file})    { die "ERROR didn't all line in $file"; }
# output:
#             num-pass  num-fail                        T1     T2 
#  filename   tot inc   tot  inc  num-w-Ns num-w-oth tot inc tot inc ...

  printf("%-50s  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d  %4d %4d\n",
         $file, 
         $pass_all_ct_H{$file},    $pass_incomp_ct_H{$file},
         $fail_all_ct_H{$file},    $fail_incomp_ct_H{$file},
         $num_Ns_all_ct_H{$file},  $num_Ns_incomp_ct_H{$file},
         $num_oth_all_ct_H{$file}, $num_oth_incomp_ct_H{$file},
         $t1_all_ct_H{$file},      $t1_incomp_ct_H{$file},      
         $t2_all_ct_H{$file},      $t2_incomp_ct_H{$file},      
         $t3_all_ct_H{$file},      $t3_incomp_ct_H{$file},      
         $t4_all_ct_H{$file},      $t4_incomp_ct_H{$file},      
         $t5_all_ct_H{$file},      $t5_incomp_ct_H{$file},      
         $t6_all_ct_H{$file},      $t6_incomp_ct_H{$file},      
         $t7_all_ct_H{$file},      $t7_incomp_ct_H{$file},      
         $t8_all_ct_H{$file},      $t8_incomp_ct_H{$file});      
}
print("#\n");
foreach my $line (@extra_lines_A) { 
  print $line;
}
close(LIST);
