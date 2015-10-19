#!/usr/bin/perl -w

# This perl script will create a report file in latex

use strict;
use warnings;

open (MYOUTFILE,">ResultReport.tex");
#if ( $ARGV[0] eq '' ) {
#    die "A fasta file from DeCypher has to be provided as first input\n";}

#if ( $ARGV[1] eq '' ) {
#    die "A text file from DeCypher hast to be provided as second input\n";}

#if ( $ARGV[2] eq '' ) {
#    die "The filepath of the Pfam-A database has to be provided as third input\n";}
print MYOUTFILE "\\documentclass\[a4paper\]\{article\}\n";
print MYOUTFILE "\\usepackage\{graphicx\}\n";
print MYOUTFILE "\\usepackage\{hyperref\}\n";
print MYOUTFILE "\\hypersetup\{\n";
print MYOUTFILE "colorlinks   = true,\n";
print MYOUTFILE "urlcolor     = blue,\n";
print MYOUTFILE "linkcolor    = blue,\n";
print MYOUTFILE "citecolor   = red\n";
print MYOUTFILE "\}\n";
print MYOUTFILE "\\usepackage\[top=2cm, bottom=2cm, left=2.5cm, right=2.5cm\]\{geometry\}\n";
print MYOUTFILE "\\usepackage\{longtable\}\n";

print MYOUTFILE "\\begin\{document\}\n";
print MYOUTFILE "\\title\{Report\}\n";
print MYOUTFILE "\\author\{Protein prediction\}\n";
print MYOUTFILE "\\maketitle\n";
print MYOUTFILE "\\section\{Queries\}\n";
print MYOUTFILE "Accession numbers and descriptions of the queries used\:\\\\\\\\\n";
open(QUERYINFO,"Cfin.queryinfo");

while (<QUERYINFO>){
    chomp;
    my $line=$_;
    if ($line =~ /^Query\:\s(.*)\t(.*)$/){
	my $query=$1;
	my $querydescription=$2;
	$querydescription =~ s/\_/\\_/g;
	$querydescription =~ s/\=/\\=/g;
	my $searchquery=$query;
	$query =~ s/\_/\\_/g;
	print MYOUTFILE "\\href\{http://www.ncbi.nlm.nih.gov/protein/$searchquery\}{$query} - $querydescription\\\\\n";
    }
}
close QUERYINFO;

print MYOUTFILE "\\section\{Unique contig hits in the transcriptome of your target species}\n";
print MYOUTFILE "The following sections will summarize the results for each of these contigs\:\\\\\\\\\n";
open(COMPINFO,"Cfin.compinfo");
my @comphits;
my @comphitslatex;
while (<COMPINFO>){
    chomp;
    my $line=$_;
    if ($line =~ /^hit\:\s(.*)$/){
	my $hit=$1;
	push @comphits,$hit;
	$hit =~ s/\_/\\_/g;
	push @comphitslatex,$hit;
	print MYOUTFILE "$hit\\\\\n";
    }
}
close COMPINFO;

    print MYOUTFILE "\\\\The fasta file \\href\{UniqueComphits.fasta\}\{UniqueComphits.fasta\} contains the polypeptide sequences of these contigs.\\\\\n";        
    print MYOUTFILE "The fasta file \\href\{LongestPolypeptide.fasta\}\{LongestPolypeptide.fasta\} contains the polypeptide sequences of the longest open reading frames for these contigs.\\\\\n";        

foreach (@comphitslatex){
    my $actualcomphit=$_;
    my $comphit= $actualcomphit;
    $comphit =~ s/\\_/\_/g;
    print MYOUTFILE "\\clearpage\n";        
    print MYOUTFILE "\\section\{$actualcomphit\}\n";

############### Include a table on the information to which querys
############### this comphit showed up as a hit and in which quality
############### print MYOUTFILE "\\begin\{center\}\n";
    print MYOUTFILE "\\textbf\{Contig $actualcomphit showed up as a hit to the following queries (sorted by e-values). The targetframe refers to the  contig of your target species.\}\\\\\\\\\n";
    print MYOUTFILE "\\begin\{longtable\}\{c c c r r\}\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "Query \& Rank \& E\-value \& Queryframe \& Targetframe\\\\\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\endfirsthead\n";
    print MYOUTFILE "\\multicolumn\{5\}\{c\}\n";
    print MYOUTFILE "\{\\tablename\ -- \\textit\{Continued from previous page\}\}\\\\\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "Query \& Rank \& E\-value \& Queryframe \& Targetframe\\\\\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\endhead\n";
    print MYOUTFILE "\\multicolumn\{5\}\{r\}\n";
    print MYOUTFILE "\{\\textit\{Continued on next page\}\} \\\\\n";
    print MYOUTFILE "\\endfoot\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\endlastfoot\n";

    open(COMPINFO2,"Cfin.compinfo");
    my @queries;
    my @ranks;
    my @evalues;
    my @queryframes;
    my @targetframes;
    while (<COMPINFO2>){
	chomp;
	my $line=$_;
	if ($line =~ /^$comphit\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)$/){
	    push @queries,$1;
	    push @ranks,$2;
	    push @evalues,$3;
	    push @queryframes,$4;
	    push @targetframes,$5;
	}
    }
    close COMPINFO2;
    # rankig
    my @list_order = sort { $evalues[$a] <=> $evalues[$b] } 0 .. $#evalues;
    @queries=@queries[@list_order];
    @ranks=@ranks[@list_order];
    @evalues=@evalues[@list_order];
    @queryframes=@queryframes[@list_order];
    @targetframes=@targetframes[@list_order];
    my $e;
    my $elength=@evalues;
	for ($e=0;$e<$elength;$e++){
	    my $q=$queries[$e];
	    $q =~ s/\_/\\_/g;
	    print MYOUTFILE "\\href\{http://www.ncbi.nlm.nih.gov/protein/$queries[$e]\}{$q} \& $ranks[$e] \& $evalues[$e] \& $queryframes[$e] \& $targetframes[$e] \\\\\n";
    }
    print MYOUTFILE "\\end\{longtable\}\\\\\n";
#    print MYOUTFILE "\\end\{center\}\n";    



################# Include information on best hit in Arthropod
################# proteins on NCBI
    
    print MYOUTFILE "\\\\ \\\\ \\textbf\{Best blastp hit in the protein database:\}\\\\\\\\\n";
    my $besthits;
    open(BESTHIT,"blastp.besthits");
    while (<BESTHIT>){
	chomp;
	my $line=$_;
	if ($line =~ /^$comphit\t(.*)\t(.*)\t(.*)\t(.*)/g){
	    my $besthit=$1;
	    $besthits=$besthit;
	    my $desc=$2;	    
	    my $s=$3;
	    my $ev=$4;
	    if ($besthit =~ /\_/){
		$besthit =~ s/\_/\\_/;
	    }
	    if ($desc =~ /\_/){
		$desc =~ s/\_/\\_/;
	    }
	    if ($besthits){
		print MYOUTFILE "Accession: \\href\{http://www.ncbi.nlm.nih.gov/protein/$besthits\}\{$besthit}\\\\\n";
		print MYOUTFILE "Description: $desc\\\\\n";
		print MYOUTFILE "Score: $s\\\\\n";
		print MYOUTFILE "E-value: $ev\\\\\n";
	    }
	    else {
		print MYOUTFILE "No best hit found !\\\\\n";
	    }
	}
    }
    close BESTHIT;
    



############## Include information on the identity and similarity from the mafft alignment
    if ($besthits){    
	open(MAFFTOUT,"parsemafft.out");
	while (<MAFFTOUT>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^$comphit\t.*\t(.*)\t(.*)/g){
		my $identity=sprintf "%.2f", $1;
		my $similarity=sprintf "%.2f", $2;
		print MYOUTFILE "\\% Identity (based on alignment): $identity\\\\\n";
		print MYOUTFILE "\\% Similarity (based on alignment): $similarity\\\\\n";
	    }
	}
	close MAFFTOUT;
	
########### Make a link to the clustalw alignment file
	
	my $printcomphit=$comphit;
	$printcomphit =~ s/\_/\\_/g;
	my $printbesthit=$besthits;
	$printbesthit =~ s/\_/\\_/g;
	print MYOUTFILE "Link to alignment file:\\\\\n";
	print MYOUTFILE "\\href\{$comphit\_\_$besthits\_\_MAFFTalignmentClustalW\.txt\}\{$printcomphit\\_\\_$printbesthit\\_\\_MAFFTalignmentClustalW\.txt\}\\\\\n";
	

############# Include information on the protein domains in the Cfin contig and the best hit.
	print MYOUTFILE "\\\\\\textbf\{Domain predictions. Start and End refer to the non-aligned protein-sequences. Slength: length of the query sequence; Dlength: Length of the domain.\}\\\\\\\\\n";    
	
	print MYOUTFILE "\\begin\{longtable\}\{l r l l r r r l\}\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "Sequence \& Slength \& Domain \& Accession \& Dlength \& Start \& End \& i\-E\-value \\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endfirsthead\n";
	print MYOUTFILE "\\multicolumn\{8\}\{c\}\n";
	print MYOUTFILE "\{\\tablename\ -- \\textit\{Continued from previous page\}\}\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "Sequence \& Slength \& Domain \& Accession \& Dlength \& Start \& End \& i\-E\-value \\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endhead\n";
	print MYOUTFILE "\\multicolumn\{8\}\{r\}\n";
	print MYOUTFILE "\{\\textit\{Continued on next page\}\} \\\\\n";
	print MYOUTFILE "\\endfoot\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endlastfoot\n";

	
	open(PROTEINDOMAINS,"hmmscan.out");
	while (<PROTEINDOMAINS>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /$comphit/g){
		$line =~ /.*\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/g;
		my $seqlength=$1;
		my $domain=$2;
		my $accession=$5;
		my $domlength=$6;
		my $start=$3;
		my $end=$4;
		my $evalue=$10;
		$domain =~ s/\_/\\_/g;
		print MYOUTFILE "$actualcomphit \& $seqlength \& $domain \& \\href\{http://pfam.xfam.org/family/$accession\}\{$accession\} \& $domlength \& $start \& $end \& $evalue \\\\\n";
	    }
	    if ($line =~ /$besthits/g){
		$line =~ /.*\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/g;
		my $seqlength=$1;
		my $domain=$2;
		my $accession=$5;
		my $domlength=$6;
		my $start=$3;
		my $end=$4;
		my $evalue=$10;
		$domain =~ s/\_/\\_/g;
		print MYOUTFILE "$printbesthit \& $seqlength \& $domain \& \\href\{http://pfam.xfam.org/family/$accession\}\{$accession\} \& $domlength \& $start \& $end \& $evalue \\\\\n";
}
	}
	close PROTEINDOMAINS;
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\end\{longtable\}\\\\\n";


############### Include the figure of alignment and protein domain
	print MYOUTFILE "\\begin\{figure\}\[h\]\n";
	print MYOUTFILE "\\centering\n";
	print MYOUTFILE "\\includegraphics\[width\=\.6\\textwidth\]\{\{$comphit\_Pfam\}.png\}\n";
	print MYOUTFILE "\\caption\{Protein domains predicted along the protein sequence of your target species derived from $actualcomphit and its aligned best blastp hit $printbesthit.\}\n";
	print MYOUTFILE "\\end\{figure\}\n";
    }
    
    
############## Include information on the EST hits in the Cfin database
    print MYOUTFILE "\\\\\\textbf\{EST hits\:\}\\\\\\\\\n";    
    
    print MYOUTFILE "\\begin\{longtable\}\{l r r r l r r\}\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "EST \& Length \& Score \& Identity \& E\-value \& Alignmentstart \& Frame\\\\\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\endfirsthead\n";
    print MYOUTFILE "\\multicolumn\{7\}\{c\}\n";
    print MYOUTFILE "\{\\tablename\ -- \\textit\{Continued from previous page\}\}\\\\\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "EST \& Length \& Score \& Identity \& E\-value \& Alignmentstart \& Frame\\\\\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\endhead\n";
    print MYOUTFILE "\\multicolumn\{7\}\{r\}\n";
    print MYOUTFILE "\{\\textit\{Continued on next page\}\} \\\\\n";
    print MYOUTFILE "\\endfoot\n";
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\endlastfoot\n";
    
    
    
    open(VETTINGOUT,"Vetting.out");
    while (<VETTINGOUT>){
	chomp;
	my $line=$_;
	if ($line =~ /^$comphit\t.*\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t.*\t(.*)/g){
	    my $est=$1;
	    my $searchest=$est;
	    my $length=$2;
	    my $score=$3;
	    my $identity=$4;
	    my $evalue=$5;
	    my $start=$6;
	    my $readingframe=$7;
	    if ($est =~ /\_/){
		$est =~ s/\_/\\_/;
	    }
	    $identity =~ s/%/\\%/;
	    print MYOUTFILE "\\href\{http://www.ncbi.nlm.nih.gov/nucest/$searchest\}\{$est\} \& $length \& $score \& $identity \& $evalue \& $start \& $readingframe \\\\\n";
	}
    }
    close VETTINGOUT;
    print MYOUTFILE "\\hline\n";
    print MYOUTFILE "\\end\{longtable\}\\\\\n";
    

}



print MYOUTFILE "\\end\{document\}\n";


#\section{Conclusion}
#Write your conclusion here.




close MYOUTFILE;
