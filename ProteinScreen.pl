#!/usr/bin/perl -w

use strict;
use warnings;

if ( $ARGV[0] eq '' ) {
    die "A fasta file with query sequences has to be provided as first input\n";}

if ( $ARGV[1] eq '' ) {
    die "The filepath of the Pfam-A database has to be provided as second input\n";}

if ( $ARGV[2] eq '' ) {
    die "The fasta file of the local transcriptome database needs to be provided as third argument to the script\n"; }

if ( $ARGV[3] eq '' ) {
    die "The protein database (nr) of taxa related to your target species has to be provided as fourth argument to the script\n"; }

if ( $ARGV[4] eq '' ) {
    die "fasta file with EST sequence(s) of your target species has to be provided as fifth argument to the
script\n"; }

my $retVal0 = `perl LocalDatabaseSearch.pl $ARGV[0] $ARGV[2]`;
my $retVal01 = `perl Parse_LocalDatabaseSearch.pl tblastn.LocalDatabase`;
my $retVal02 = `perl BestHitsToFasta.pl Cfin.compinfo $ARGV[2]`;
my $retVal03 = `perl translation.pl Comphits.fasta Cfin.compinfo`;
my $retVal1 = `perl DeduplicateFasta.pl TranslatedCompHits.fasta`;
my $retVal3 = `perl Parse_DeCypher2.pl Cfin.compinfo deduplicated.fasta`;
my $retVal4 = `perl PeptideExtraction.pl UniqueCompHits.fasta`;
my $retVal5 = `perl ReciprocalArthropodaBlast.pl LongestPolypeptide.fasta $ARGV[3]`;
my $retVal6 = `perl Parse_ArthropodaBlastp.pl  Blastp.outfiles $ARGV[3]`;
my $retVal7 = `perl HMMSCAN.pl LongestPolypeptide.fasta besthits.fasta $ARGV[1]`;
my $retVal8 = `perl MAFFT.pl LongestPolypeptide.fasta besthits.fasta blastp.besthits`;
my $retVal9 = `perl Parse_mafft.pl mafftalignments.outfilenames`;
my $retVal10 = `perl Vetting.pl LongestPolypeptide.fasta $ARGV[4]`;
my $retVal11 = `perl Parse_Vetting.pl tblastn.ESTout`;
my $retVal12 = `R CMD BATCH Pfam.r`;
my $retVal13 = `perl CreateReport.pl`;
my $retVal14 = `pdflatex --interaction=nonstopmode ResultReport.tex`;
my $retVal15 = `pdflatex --interaction=nonstopmode ResultReport.tex`;
