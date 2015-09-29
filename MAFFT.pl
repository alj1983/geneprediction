#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;

if ( $ARGV[0] eq '' ) {
    die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. The first contains the contigs of your target species and the second contains the sequences of the best blastp hits\n"; }

if ( $ARGV[1] eq '' ) {
    die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. The first contains the contigs of your target species and the second contains the sequences of the best blastp hits\n"; }

if ( $ARGV[2] eq '' ) {
    die "Provide the file ending with .besthits as third argument\n",}


# First loading in the information which fasta sequences form a pair 
open(IN, $ARGV[2]) or die "could not open file $ARGV[2]\n";

my @queries;
my @besthits;
while (<IN>) {
    chomp;
    if ($.!=1){# skip the header line
	my @row = split(/\t/);
	push @queries,$row[0];
	push @besthits,$row[1];
    }
}
close IN;

#### Extract the protein sequences for the query-besthit pairs and
#### submit them for alignment to mafft

open(OUT, ">mafftalignments.outfilenames"); # this will list the outputfilenames in one file.

# Creating a temporary output file that stores the sequence pairs in fasta format
my $t;
my $tlength=@queries;
for ($t=0;$t<$tlength;$t++){
    open (TEMPORARY, ">temporary_fastapair.fasta");
    print "$t\n";
## Getting the query sequence
    my $query_in = Bio::SeqIO->new(
	-file   => "<$ARGV[0]",
	-format => 'Fasta',
	);
    
    my $actualquery;
    while (my $seq = $query_in->next_seq) {
	my $fastaid=$seq->display_id;
	if ($fastaid =~ /.*$queries[$t].*/){
	    $actualquery = $seq->seq();
	}
    }
    
    
    print TEMPORARY ">$queries[$t]\n";
    print TEMPORARY "$actualquery\n\n";
    
## getting the besthit sequence
    my $besthit_in = Bio::SeqIO->new(
	-file   => "<$ARGV[1]",
	-format => 'Fasta',
	);
    
    my $actualbesthit;
    while (my $seq = $besthit_in->next_seq) {
	my $fastaid=$seq->display_id;
	if ($fastaid =~ /.*$besthits[$t].*/){
	    $actualbesthit = $seq->seq();
	}
    }
    print TEMPORARY ">$besthits[$t]\n";
    print TEMPORARY "$actualbesthit\n\n";

    close TEMPORARY;

    my $fastapairname1="$queries[$t]\_\_$besthits[$t]\_\_MAFFTalignmentClustalW\.txt";
    print OUT "$fastapairname1\n";
    my $maffrun1 = `mafft --clustalout temporary_fastapair.fasta  > $fastapairname1`;

    my $fastapairname2="$queries[$t]\_\_$besthits[$t]\_\_MAFFTalignment\.fasta";
    my $maffrun2 = `mafft temporary_fastapair.fasta  > $fastapairname2`;
    
    


}
close OUT;
