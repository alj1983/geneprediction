#!/usr/bin/perl -w

package ProteinScreen;
use match::smart; ## Replaces the smartmatch operator ~~
use strict;
use warnings;
use Bio::SeqIO;
use LWP::Simple;
use List::MoreUtils qw(uniq); # This allows to find all unique values in an array


### Function definitions

## local_database_search
sub local_database_search {

}      


## parse_local_database_search
# Script to extract Best hit in the local transcriptome/genome with
# information such as evalue and reading frame from tblastn output and
# then selecting only those with an e-value <1e-05
sub parse_local_database_search {


}



## best_hits_to_fasta

sub best_hits_to_fasta {

}

## translation
sub translation {

}

## deduplicate_fasta
sub deduplicate_fasta {

}

## parse_decypher2
sub parse_decypher2 {

}

## peptide_extraction
sub peptide_extraction {

}

## reciprocal_blast
sub reciprocal_blast {

}

## parse_blast
sub parse_blastp {


}

## hmmscan
sub hmmscan {


}

## mafft
sub mafft {


}

## parse_mafft
sub parse_mafft {


}

## vetting
sub vetting {
 
}

## parse_vetting
sub parse_vetting {
 

}


## create_resport
sub create_report {

}


### Function calls
sub ProteinScreen {
# argument 0: A fasta file with query sequences
# argument 1: The filepath of the Pfam-A database
# argument 2: The fasta file of the local transcriptome database
# argument 3: The protein database (nr) of taxa related to your target species
# argument 4: fasta file with EST sequence(s) of your target species
local_database_search($_[0], $_[2]);
parse_local_database_search("tblastn.LocalDatabase");
best_hits_to_fasta("Cfin.compinfo");
translation("Comphits.fasta", "Cfin.compinfo");
deduplicate_fasta("TranslatedCompHits.fasta");
parse_decypher2("Cfin.compinfo", "deduplicated.fasta");
peptide_extraction("UniqueCompHits.fasta");
reciprocal_blast("LongestPolypeptide.fasta", $_[3]);
parse_blastp("Blastp.outfiles");
hmmscan("LongestPolypeptide.fasta", "besthits.fasta", $_[1]);
mafft("LongestPolypeptide.fasta", "besthits.fasta", "blastp.besthits");
parse_mafft("mafftalignments.outfilenames");
vetting("LongestPolypeptide.fasta", $_[4]);
parse_vetting("tblastn.ESTout");
#my $retVal12 = `R CMD BATCH Pfam.r`;
create_report();
#my $retVal14 = `pdflatex --interaction=nonstopmode ResultReport.tex`;
#my $retVal15 = `pdflatex --interaction=nonstopmode ResultReport.tex`;
}

#ProteinScreen::ProteinScreen("/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/Fishes_reproduction_embl.fasta", "/home/alj/Dropbox.personal/Dropbox/Pfam/Pfam-A.hmm", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nonredundant_dbv1.fasta", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nr_bonyfishes", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/GuppyESTs.fasta");

1;

# Use the scripts in the guppy folder instead of the ones I have now. I think that will solve a lot of problems: 
# https://www.dropbox.com/work/Projects/2014Guppy/201408GuppyReproductionGenePrediction
# I downloaded them already and they are in the same folder as the actual script

# Replace $ARGV with $_


# Need to fix the uninitialized value warning?!

#XXX How to include Pfam.r execution??


# Running R script from within perl: http://www.perlmonks.org/?node_id=1009021
# XX Create reports in markdown and then convert it to pdf and html!

#XX Create usage of model in separate script file and use subroutines with NAM::function(args);
## XX Make an argument for an optional output folder
# XX Export only the ProteinScreen function and make the others inaccessible

# http://www.perlmonks.org/?node_id=304000
