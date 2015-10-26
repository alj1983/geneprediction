#!/usr/bin/perl -w
use strict;
use ProteinScreen;




#ProteinScreen::ProteinScreen("/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/Fishes_reproduction_embl.fasta", "/home/alj/Dropbox.personal/Dropbox/Pfam/Pfam-A.hmm", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nonredundant_dbv1.fasta", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nr_bonyfishes", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/GuppyESTs.fasta");

#ProteinScreen::parse_local_database_search("tblastn.LocalDatabase");

#ProteinScreen::best_hits_to_fasta("Cfin.compinfo","/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nonredundant_dbv1.fasta");


#ProteinScreen::translation("Comphits.fasta", "Cfin.compinfo");
#Use of uninitialized value $amino_acid in concatenation (.) or string at ProteinScreen.pm line 681, <GEN0> line 351.

#ProteinScreen::deduplicate_fasta("TranslatedCompHits.fasta");

#ProteinScreen::parse_decypher2("Cfin.compinfo", "deduplicated.fasta");

#ProteinScreen::peptide_extraction("UniqueCompHits.fasta");
#ProteinScreen::reciprocal_blast("LongestPolypeptide.fasta", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nr_bonyfishes");
ProteinScreen::parse_blastp("Blastp.outfiles", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nr_bonyfishes");
