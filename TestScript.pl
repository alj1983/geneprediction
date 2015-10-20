#!/usr/bin/perl -w
use strict;
use ProteinScreen;




ProteinScreen::ProteinScreen("/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/Fishes_reproduction_embl.fasta", "/home/alj/Dropbox.personal/Dropbox/Pfam/Pfam-A.hmm", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nonredundant_dbv1.fasta", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nr_bonyfishes", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/GuppyESTs.fasta");

#ProteinScreen::parse_local_database_search("tblastn.LocalDatabase");

#ProteinScreen::best_hits_to_fasta("Cfin.compinfo");
#best_hits_to_fasta("Cfin.compinfo");
