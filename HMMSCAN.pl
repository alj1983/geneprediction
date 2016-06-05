#!/usr/bin/perl -w

use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



if ( $ARGV[0] eq '' ) {
    die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. One contains the contigs of your target species the other contains the sequences of the best blastp hits\n"; }

if ( $ARGV[1] eq '' ) {
    die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. One contains the contigs of your target species the other contains the sequences of the best blastp hits\n"; }

if ( $ARGV[2] eq '' ) {
    die "Provide the Pfam path (including the file) as third argument\n",}


open my $fh, '<', $ARGV[0] or die "error opening $ARGV[0]: $!";
my $data1 = do { local $/; <$fh> };

open my $fh2, '<', $ARGV[1] or die "error opening $ARGV[1]: $!";
my $data2 = do { local $/; <$fh2> };

open (MYOUTFILE, ">alldata.fasta");
my $alldata=$data1.$data2;
print MYOUTFILE $alldata;
close MYOUTFILE;

# running hmmscan from the command line
my $output = `hmmscan --domtblout hmmscanoutput.tab $ARGV[2] alldata.fasta`;
open (MYOUTFILE2, ">hmmscan.out");
print MYOUTFILE2 "queryname\tlen\tdomname\tbegin\tend\tdomaccession\tdomlength\tquerylength\tEval\tscore\tiEval\tdomscore\n";


open(IN, 'hmmscanoutput.tab') or die "could not open file hmmscanoutput.tab\n";
# Since one domain can be reported more than once, I will select only the unique ones
my @uniquelines;
while (<IN>) {
    unless (/^\#/) {    # avoid all lines beginning
	# with the '#' character
	my @columns  = split(/ +/);
	my $domname  = $columns[0];
	my $queryname = $columns[3];
#	$queryname1 =~ /.*\|.*$/g;
	#$queryname    =~ s/.*\|//;
	my $len      = $columns[5];
	my $begin    = $columns[19]; # This is the domain envelope that is also reported in SMART.
	my $end      = $columns[20];
 	my $domaccession      = $columns[1];
 	my $domlength      = $columns[2];
 	my $querylength      = $columns[5];
 	my $Eval      = $columns[6];
 	my $score      = $columns[7];
 	my $iEval     = $columns[12];
 	my $domscore    = $columns[13];
	
	if ( $iEval < 1e-10 ) {
	    push @uniquelines, "$queryname\t$len\t$domname\t$begin\t$end\t$domaccession\t$domlength\t$querylength\t$Eval\t$score\t$iEval\t$domscore\n";
	}
    }
}
close IN;

my @unique = uniq @uniquelines;
foreach ( @unique ) {
    print MYOUTFILE2 $_;
}

close MYOUTFILE2;
