#!/usr/bin/perl -w

# Script to translate a fasta file of DNA sequences to a fasta file of protein sequences with at a given reading frame
# Adjusted from Samuelsson T. ,  Genomics and Bioinformatics - An Introduction to Programming Tools for Life Scientists
use strict;
use warnings;
use Bio::SeqIO;



# The genetic code is represented as a hash. Each codon (hash key)
# is associated with an amino acid (hash value).
# A stop codon is shown as '*'
my %code = (
     'UUU' => 'F', 'UUC' => 'F', 'UUA' => 'L', 'UUG' => 'L',
     'CUU' => 'L', 'CUC' => 'L', 'CUA' => 'L', 'CUG' => 'L',
     'AUU' => 'I', 'AUC' => 'I', 'AUA' => 'I', 'AUG' => 'M',
     'GUU' => 'V', 'GUC' => 'V', 'GUA' => 'V', 'GUG' => 'V',
     'UCU' => 'S', 'UCC' => 'S', 'UCA' => 'S', 'UCG' => 'S',
     'CCU' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
     'ACU' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
     'GCU' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
     'UAU' => 'Y', 'UAC' => 'Y', 'UAA' => '*', 'UAG' => '*',
     'CAU' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
     'AAU' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
     'GAU' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
     'UGU' => 'C', 'UGC' => 'C', 'UGA' => '*', 'UGG' => 'W',
     'CGU' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
     'AGU' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
     'GGU' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G' 
);

if ($ARGV[0] eq '') { die "A fasta input file is needed"};
if ($ARGV[1] eq '') { die "The file Cfin.compinfo is needed as second input file"};


open(COMPINFO, $ARGV[1]) or die "Can not open file $ARGV[1]";

# Get the comp ids and in which target frame they shall be translated
my @comps;
my @frames;
while(<COMPINFO>) {
    chomp;
    my $actualline=$_;
    if ($actualline =~ /\S*\t\S*\t\d*\t\S*\t\S*\t\S*/){
	my @columns=split(/\t/);
	push @comps,$columns[0];
	push @frames,$columns[5];
    }
}
close COMPINFO;

my $compnumber=$#comps;

# Split these lines and get only compID and targetframe

my $fastafile   = $ARGV[0];
my $seqio  = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");
open(MYOUTFILE, ">TranslatedCompHits.fasta"); #open for write

# Get only the unique sequence * frame combinations
my %unique;      
my @uniquecomps;
my @uniqueframes;
for (my $a = 0 ; $a < $compnumber+1 ; $a++){
    my $hashkey ="$comps[$a]"."$frames[$a]";
    unless(exists($unique{$hashkey})) {
	push @uniquecomps,$comps[$a];
	push @uniqueframes,$frames[$a];
	$unique{$hashkey} +=1;	  
    }
}
my $uniquecompnumber=$#uniquecomps;

while(my $seqs = $seqio->next_seq) {
  my $id  = $seqs->display_id;
  my $seq = $seqs->seq;
  my $translation;
  for (my $b = 0 ; $b < $uniquecompnumber+1 ; $b++){
      #print "$uniquecomps[$b]\n";
      if ($id =~ /$uniquecomps[$b]/){
	  #print "ok\n";
	  print MYOUTFILE ">$id\n";    
	  # translastion for positive readframes
	  if ($uniqueframes[$b]>0){
	      my $rnaseq = $seq;
	      $rnaseq =~ tr/T/U/;       # the RNA sequence
	      my $startpoint = (abs $uniqueframes[$b]) - 1;
	      for ( my $i = $startpoint ; $i < length($rnaseq) - 2 ; $i = $i + 3 ) {
		  my $codon = substr( $rnaseq, $i, 3 );
		  my $amino_acid = $code{$codon};
		  $translation=$translation.$amino_acid;
	      }
	      print MYOUTFILE "$translation\n";    
	  }
	  
	  
	  # translation for negative read frames
	  if ($uniqueframes[$b]<0){
	      my $revseq=reverse($seq);
	      $revseq =~ tr/natgcbdkrvhmyxwsNATGCBDKRVHMYXWS/ntacgvhmybdkrxswNTACGVHMYBDKRXSW/;
	      my $rnaseq = $revseq;
	      $rnaseq =~ tr/T/U/;       # the RNA sequence
	      my $startpoint = (abs $uniqueframes[$b]) - 1;
	      for ( my $i = $startpoint ; $i < length($rnaseq) - 2 ; $i = $i + 3 ) {
		  my $codon = substr( $rnaseq, $i, 3 );
		  my $amino_acid = $code{$codon};
		  $translation=$translation.$amino_acid;
		  
	      }
	      print MYOUTFILE "$translation\n";    
	  }
	  
	  
	  
	  
      }
  
  }
  
}

close MYOUTFILE;



