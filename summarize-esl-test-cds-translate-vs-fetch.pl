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
my $do_subset   = 0;  # changed to '1' with -subset, run only first $nsubset tests instead of all $nall
# these need to be manually kept in sync with variables of the same name in esl-test-cds-translate-vs-fetch.pl
my $nall        = 11; # number of total tests
my $nsubset     = 9;  # number of tests run if -subset is used

&GetOptions( "oroot=s" => \$outroot, 
             "subset"  => \$do_subset);

my $usage;
$usage  = "summarize-esl-test-cds-translate-vs-fetch.pl [OPTIONS] <file with list of multiple esl-test-cds-translate-vs-fetch.pl output files>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-oroot <s>: output sequences that fail each test to <s>.T<n>.fail-list and seqs that fail any test to <s>.any.fail-list\n";
$usage .= "\t\t-subset:    esl-test-cds-translate-vs-fetch was run with the -subset option\n";

if(scalar(@ARGV) != 1) { die $usage; }
($in_listfile) = @ARGV;

if(! -s $in_listfile) { die "ERROR, $in_listfile does not exist"; }

my $last_test = ($do_subset) ? $nsubset : $nall;

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
my %t9_all_ct_H = ();
my %t10_all_ct_H = ();
my %t11_all_ct_H = ();

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
my %t9_incomp_ct_H = ();
my %t10_incomp_ct_H = ();
my %t11_incomp_ct_H = ();

# counts of totals across all files
my $pass_all_ct = 0;
my $fail_all_ct = 0;
my $pass_incomp_ct = 0;
my $fail_incomp_ct = 0;
my $nlist_n  = 0;
my $nlist_ab = 0;
my $nfail_t1 = 0;
my $nfail_t2 = 0;
my $nfail_t3 = 0;
my $nfail_t4 = 0;
my $nfail_t5 = 0;
my $nfail_t6 = 0;
my $nfail_t7 = 0;
my $nfail_t8 = 0;
my $nfail_t9 = 0;
my $nfail_t10 = 0;
my $nfail_t11 = 0;
my $nfail_any = 0;

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
my $t9_fail_file  = undef;
my $t10_fail_file  = undef;
my $t11_fail_file  = undef;
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
  $t9_fail_file  = $outroot . ".T9.fail-list";
  $t10_fail_file = $outroot . ".T10.fail-list";
  $t11_fail_file = $outroot . ".T11.fail-list";
  $any_fail_file = $outroot . ".any.fail-list";
  open(OUTN,   ">" . $N_list_file)   || die "ERROR unable to open $N_list_file for writing"; 
  open(OUTAB,  ">" . $ab_list_file)  || die "ERROR unable to open $ab_list_file for writing"; 
  open(OUTT1,  ">" . $t1_fail_file)  || die "ERROR unable to open $t1_fail_file for writing"; 
  open(OUTT2,  ">" . $t2_fail_file)  || die "ERROR unable to open $t2_fail_file for writing"; 
  open(OUTT3,  ">" . $t3_fail_file)  || die "ERROR unable to open $t3_fail_file for writing"; 
  open(OUTT4,  ">" . $t4_fail_file)  || die "ERROR unable to open $t4_fail_file for writing"; 
  open(OUTT5,  ">" . $t5_fail_file)  || die "ERROR unable to open $t5_fail_file for writing"; 
  open(OUTT6,  ">" . $t6_fail_file)  || die "ERROR unable to open $t6_fail_file for writing"; 
  open(OUTT7,  ">" . $t7_fail_file)  || die "ERROR unable to open $t7_fail_file for writing"; 
  open(OUTT8,  ">" . $t8_fail_file)  || die "ERROR unable to open $t8_fail_file for writing"; 
  open(OUTT9,  ">" . $t9_fail_file)  || die "ERROR unable to open $t9_fail_file for writing"; 
  if(! $do_subset) { 
    open(OUTT10,  ">" . $t9_fail_file)  || die "ERROR unable to open $t9_fail_file for writing"; 
    open(OUTT11,  ">" . $t9_fail_file)  || die "ERROR unable to open $t9_fail_file for writing"; 
  }
  open(OUTANY, ">" . $any_fail_file) || die "ERROR unable to open $any_fail_file for writing"; 
}

my $filename_w = 63;
my $header_line1 = sprintf("%-*s  %13s  %13s\n", 
                           $filename_w, "#", "num-pass", "num-fail");
my $header_line2 = "";
if($do_subset) { 
  $header_line2 = sprintf("%-*s  %6s %6s  %6s %6s  %9s  %9s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s\n", 
                          $filename_w, "# filename", "tot", "inc", "tot", "inc", "num-w-Ns", "num-w-oth", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9");
}
else { 
  $header_line2 = sprintf("%-*s  %6s %6s  %6s %6s  %9s  %9s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s\n", 
                          $filename_w, "# filename", "tot", "inc", "tot", "inc", "num-w-Ns", "num-w-oth", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11");
}

print $header_line1; 
print $header_line2; 

open(LIST, $in_listfile) || die "ERROR unable to open $in_listfile for reading";
my $file_ctr = 0;
while(my $line = <LIST>) { 
  chomp $line;
  $file_ctr++;
  my $file = $line;
  my $file2print = $file;
  $file2print =~ s/^.+\///;
  if(! -s $file) { die "ERROR, unable to open $file listed in $in_listfile"; }
  # Lines we parse (except first one):
  ##protein-accession  nt-accession  mincoord  maxcoord  tr-len  nexon  str  incomplete?  start   num-Ns  num-oth     T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11   pass/fail
  #AAO63223.1          AY173956            37      1545    1509      1    +           no      1        0        2      0    0    0    0    0    0    0    0    0    0    0   pass     
  #AAO63240.1          AY173958           138      1616    1479      1    +           no      1        0        1      0    0    0    0    0    0    0    0    0    0    0   pass     
  #AAV90743.1          AY819715           157      1656    1500      1    +           no      1        0        6      0    0    0    0    0    0    0    0    0    0    0   pass     
  #
  #OR
  # 
  # if ($do_subset): 
  ##protein-accession  nt-accession  mincoord  maxcoord  tr-len  nexon  str  incomplete?  start   num-Ns  num-oth     T1   T2   T3   T4   T5   T6   T7   T8   T9    pass/fail
  #AAO63223.1          AY173956            37      1545    1509      1    +           no      1        0        2      0    0    0    0    0    0    0    0    0    pass     
  #AAO63240.1          AY173958           138      1616    1479      1    +           no      1        0        1      0    0    0    0    0    0    0    0    0    pass     
  #AAV90743.1          AY819715           157      1656    1500      1    +           no      1        0        6      0    0    0    0    0    0    0    0    0    pass     
  #
  # AND
  # 
  # category    num-pass  num-fail  fract-fail
  # complete         184         0      0.0000
  # incomplete        20         0      0.0000
  # all              204         0      0.0000
  open(IN, $file) || die "ERROR unable to open $file for reading"; 
  push(@file_A, $file);
  $num_Ns_all_ct_H{$file}     = 0;
  $num_oth_all_ct_H{$file}    = 0;
  $num_Ns_incomp_ct_H{$file}  = 0;
  $num_oth_incomp_ct_H{$file} = 0;
  while($line = <IN>) { 
    chomp $line;
    if($line !~ /^\#/) { 
      #AAO63223.1                       no      1    1509        0        2      0    0    0    0    0    0    0    0    pass     
      my ($paccn, $ntaccn, $incomplete, $num_Ns, $num_oth, $t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8, $t9, $t10, $t11);
      my $have_parseable_line = 0;
      if($line =~ /^(\S+)\s+(\S+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\S\s+(\S+)\s+\d\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S\S\S\S/)
      { 
        if(! $do_subset) { die "ERROR, read line with only $nsubset pass/fail values. Did you mean to use -subset because you used -subset with esl-fetch-cds-translate-vs-fetch.pl?"; }
        ($paccn, $ntaccn, $incomplete, $num_Ns, $num_oth, $t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8, $t9) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12);
        $have_parseable_line = 1;
      }
      if($line =~ /^(\S+)\s+(\S+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\S\s+(\S+)\s+\d\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S\S\S\S/) 
      { 
        if($do_subset) { die "ERROR, read line with $nall pass/fail values but you used the -subset option. It does not appear that -subset was also used with esl-fetch-cds-translate-vs-fetch.pl?"; }
        ($paccn, $ntaccn, $incomplete, $num_Ns, $num_oth, $t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8, $t9, $t10, $t11) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14);
        $have_parseable_line = 1;
      }
      if($have_parseable_line) { 
        my $outline = ""; # what we'll print to the output files
        my @ntaccn_A = split(",", $ntaccn);
        foreach my $indi_ntaccn (@ntaccn_A) { 
          $outline .= "$paccn $indi_ntaccn\n";
        }
        if($num_Ns > 0) { 
          $num_Ns_all_ct_H{$file}++;
          if(defined $N_list_file) { print OUTN $outline; }
          if($incomplete =~ m/yes/) { $num_Ns_incomp_ct_H{$file} += $num_Ns; }
          $nlist_n++;
        }
        if($num_oth > 0) { 
          $num_oth_all_ct_H{$file}++;
          if(defined $ab_list_file) { print OUTAB $outline; }
          if($incomplete =~ m/yes/) { $num_oth_incomp_ct_H{$file} += $num_oth; }
          $nlist_ab++;
        }
        my $fail_flag = 0;
        if($t1 != 0)  { $fail_flag = 1; $nfail_t1++;  if(defined $t1_fail_file)  { print OUTT1 $outline; } }
        if($t2 != 0)  { $fail_flag = 1; $nfail_t2++;  if(defined $t2_fail_file)  { print OUTT2 $outline; } }
        if($t3 != 0)  { $fail_flag = 1; $nfail_t3++;  if(defined $t3_fail_file)  { print OUTT3 $outline; } }
        if($t4 != 0)  { $fail_flag = 1; $nfail_t4++;  if(defined $t4_fail_file)  { print OUTT4 $outline; } }
        if($t5 != 0)  { $fail_flag = 1; $nfail_t5++;  if(defined $t5_fail_file)  { print OUTT5 $outline; } }
        if($t6 != 0)  { $fail_flag = 1; $nfail_t6++;  if(defined $t6_fail_file)  { print OUTT6 $outline; } }
        if($t7 != 0)  { $fail_flag = 1; $nfail_t7++;  if(defined $t7_fail_file)  { print OUTT7 $outline; } }
        if($t8 != 0)  { $fail_flag = 1; $nfail_t8++;  if(defined $t8_fail_file)  { print OUTT8 $outline; } }
        if($t9 != 0)  { $fail_flag = 1; $nfail_t9++;  if(defined $t9_fail_file)  { print OUTT9 $outline; } }
        if(! $do_subset) { 
          if($t10 != 0) { $fail_flag = 1; $nfail_t10++; if(defined $t10_fail_file) { print OUTT10 $outline; } }
          if($t11 != 0) { $fail_flag = 1; $nfail_t11++; if(defined $t11_fail_file) { print OUTT11 $outline; } }
        }
        if($fail_flag) { $nfail_any++; if(defined $any_fail_file) { print OUTANY $outline; } }
      }
      else { 
        die "ERROR unable to parse line in $file: $line";
      }
    } # end of 'if($line !~ /^\#/)'
    if($line =~ /^\# num-fails-incomplete\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$/) { 
      if(! $do_subset) { die "ERROR, read summary line with only $nsubset pass/fail values. Did you mean to use -subset because you used -subset with esl-fetch-cds-translate-vs-fetch.pl?"; }
      # num-fails-incomplete      0    0    0    0    0    0    0    0    0
      $t1_incomp_ct_H{$file} = $1;
      $t2_incomp_ct_H{$file} = $2;
      $t3_incomp_ct_H{$file} = $3;
      $t4_incomp_ct_H{$file} = $4;
      $t5_incomp_ct_H{$file} = $5;
      $t6_incomp_ct_H{$file} = $6;
      $t7_incomp_ct_H{$file} = $7;
      $t8_incomp_ct_H{$file} = $8;
      $t9_incomp_ct_H{$file} = $9;
    }
    if($line =~ /^\# num-fails-incomplete\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$/) { 
      if($do_subset) { die "ERROR, read summary line with $nall pass/fail values but you used the -subset option. It does not appear that -subset was also used with esl-fetch-cds-translate-vs-fetch.pl?"; }
      # num-fails-incomplete      0    0    0    0    0    0    0    0    0    0    0 
      $t1_incomp_ct_H{$file}  = $1;
      $t2_incomp_ct_H{$file}  = $2;
      $t3_incomp_ct_H{$file}  = $3;
      $t4_incomp_ct_H{$file}  = $4;
      $t5_incomp_ct_H{$file}  = $5;
      $t6_incomp_ct_H{$file}  = $6;
      $t7_incomp_ct_H{$file}  = $7;
      $t8_incomp_ct_H{$file}  = $8;
      $t9_incomp_ct_H{$file}  = $9;
      $t10_incomp_ct_H{$file} = $10;
      $t11_incomp_ct_H{$file} = $11;
    }
    if($line =~ /^\# num-fails-all\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$/) { 
      if(! $do_subset) { die "ERROR, read summary line with only $nsubset pass/fail values. Did you mean to use -subset because you used -subset with esl-fetch-cds-translate-vs-fetch.pl?"; }
      # num-fails-all             0    0    0    0    0    0    0    0    0
      $t1_all_ct_H{$file} = $1;
      $t2_all_ct_H{$file} = $2;
      $t3_all_ct_H{$file} = $3;
      $t4_all_ct_H{$file} = $4;
      $t5_all_ct_H{$file} = $5;
      $t6_all_ct_H{$file} = $6;
      $t7_all_ct_H{$file} = $7;
      $t8_all_ct_H{$file} = $8;
      $t9_all_ct_H{$file} = $9;
    }
    if($line =~ /^\# num-fails-all\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$/) { 
      if($do_subset) { die "ERROR, read summary line with $nall pass/fail values but you used the -subset option. It does not appear that -subset was also used with esl-fetch-cds-translate-vs-fetch.pl?"; }
      # num-fails-all             0    0    0    0    0    0    0    0    0
      $t1_all_ct_H{$file}  = $1;
      $t2_all_ct_H{$file}  = $2;
      $t3_all_ct_H{$file}  = $3;
      $t4_all_ct_H{$file}  = $4;
      $t5_all_ct_H{$file}  = $5;
      $t6_all_ct_H{$file}  = $6;
      $t7_all_ct_H{$file}  = $7;
      $t8_all_ct_H{$file}  = $8;
      $t9_all_ct_H{$file}  = $9;
      $t10_all_ct_H{$file} = $10;
      $t11_all_ct_H{$file} = $11;
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
      if($line =~ m/^\# Test/) { push @extra_lines_A, $line . "\n"; }
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

  printf("%-*s  %6d %6d  %6d %6d  %9d  %9d  %5d  %5d  %5d  %5d  %5d  %5d  %5d  %5d  %5d%s\n",
         $filename_w,
         $file2print, 
         $pass_all_ct_H{$file},    $pass_incomp_ct_H{$file},
         $fail_all_ct_H{$file},    $fail_incomp_ct_H{$file},
         $num_Ns_all_ct_H{$file},  
         $num_oth_all_ct_H{$file}, 
         $t1_all_ct_H{$file},      
         $t2_all_ct_H{$file},      
         $t3_all_ct_H{$file},      
         $t4_all_ct_H{$file},      
         $t5_all_ct_H{$file},      
         $t6_all_ct_H{$file},      
         $t7_all_ct_H{$file},      
         $t8_all_ct_H{$file},
         $t9_all_ct_H{$file}, 
         ($do_subset) ? "" : sprintf("  %5d  %5d", $t10_all_ct_H{$file}, $t11_all_ct_H{$file}));

  $pass_all_ct    += $pass_all_ct_H{$file};
  $pass_incomp_ct += $pass_incomp_ct_H{$file};
  $fail_all_ct    += $fail_all_ct_H{$file};
  $fail_incomp_ct += $fail_incomp_ct_H{$file};
}

printf("%-*s  %6d %6d  %6d %6d  %9d  %9d  %5d  %5d  %5d  %5d  %5d  %5d  %5d  %5d  %5d%s\n",
       $filename_w,
       "# total", 
       $pass_all_ct,  $pass_incomp_ct, 
       $fail_all_ct,  $fail_incomp_ct,
       $nlist_n, $nlist_ab,
       $nfail_t1, 
       $nfail_t2, 
       $nfail_t3, 
       $nfail_t4, 
       $nfail_t5, 
       $nfail_t6, 
       $nfail_t7, 
       $nfail_t8, 
       $nfail_t9, 
       ($do_subset) ? "" : sprintf("  %5d  %5d", $nfail_t10, $nfail_t11));

print("#\n");
print("# Explanation of column headings:\n");
printf("# 'num-pass':          'tot': number of gene sequences that pass all tests (T1 through T%s); 'inc': subset that are 'incomplete' CDS\n", ($do_subset) ? "9" : "11");
printf("# 'num-fail':          'tot': number of gene sequences that fail >=1 test  (T1 through T%s); 'inc': subset that are 'incomplete' CDS\n", ($do_subset) ? "9" : "11");
print("# 'num-w-Ns':          total number of gene sequences that have >= 1 'N' ambiguous character\n");
print("# 'num-w-oth':         total number of gene sequences that have >= 1 non-'N' ambiguous character\n");
printf("# 'T1'..'T<n>'..'T%s':  total number of gene sequences that fail test <n>\n", ($do_subset) ? "9" : "11");
print("#\n");

foreach my $line (@extra_lines_A) { 
  print $line;
}
if(defined $outroot) { 
  print("#\n");
  close(OUTN);   printf("# Saved list of %4d accessions with > 0 Ns to $N_list_file\n", $nlist_n);
  close(OUTAB);  printf("# Saved list of %4d accessions with > 0 non-N ambiguous characters to $ab_list_file\n", $nlist_ab);
  close(OUTT1);  printf("# Saved list of %4d accessions that failed test 1 to $t1_fail_file\n", $nfail_t1);
  close(OUTT2);  printf("# Saved list of %4d accessions that failed test 2 to $t2_fail_file\n", $nfail_t2);
  close(OUTT3);  printf("# Saved list of %4d accessions that failed test 3 to $t3_fail_file\n", $nfail_t3);
  close(OUTT4);  printf("# Saved list of %4d accessions that failed test 4 to $t4_fail_file\n", $nfail_t4);
  close(OUTT5);  printf("# Saved list of %4d accessions that failed test 5 to $t5_fail_file\n", $nfail_t5);
  close(OUTT6);  printf("# Saved list of %4d accessions that failed test 6 to $t6_fail_file\n", $nfail_t6);
  close(OUTT7);  printf("# Saved list of %4d accessions that failed test 7 to $t7_fail_file\n", $nfail_t7);
  close(OUTT8);  printf("# Saved list of %4d accessions that failed test 8 to $t8_fail_file\n", $nfail_t8);
  close(OUTT9);  printf("# Saved list of %4d accessions that failed test 9 to $t9_fail_file\n", $nfail_t9);
  if(! $do_subset) { 
    close(OUTT10);  printf("# Saved list of %4d accessions that failed test 10 to $t10_fail_file\n", $nfail_t10);
    close(OUTT11);  printf("# Saved list of %4d accessions that failed test 11 to $t11_fail_file\n", $nfail_t11);
  }
  close(OUTANY); printf("# Saved list of %4d accessions that failed any of the %d tests to $any_fail_file\n", ($do_subset) ? "9" : "11", $nfail_any);
}
close(LIST);
