NAME

    ProteinScreen - identifying target proteins in transcriptomes

SYNOPSIS

     use strict;
     use ProteinScreen;
    
     my $queries = 'queries.fasta';
     my $pfam = 'Pfam-A.hmm';
     my $transcriptome = 'transcriptome.fasta';
     my $proteindatabase = 'nr_database';
     my $ests = 'ests.fasta';
    
     ProteinScreen::ProteinScreen($queries, $pfam, $transcriptome, $proteindatabase, $ests);

RATIONAL

    This module was designed with the purpose of automatising the process
    of identifying contigs that contain sequences similar to certain target
    proteins of interest. For example, this module was used to identify
    contigs encoding for heat shock proteins in the transcriptome of the
    copepod Calanus finmarchicus.

    XX Did I really use it on the genome? It was further used to identify
    genetic regions that contain proteins related to functions like growth
    and reproduction in the genome of guppy.

DESCRIPTION

    This module extracts contigs from a transcriptome that are in sequence
    and structure closely related to target proteins, such as for example
    heat shock proteins.

 Program and package requirements

    This program works only in a UNIX environment (Linux and Mac
    computers). The user needs to install

       =item hmmer 
      
       L<http://hmmer.org/download.html>
      
       =item BLAST+ applications 
      
       including makeblastdb, tblastn, blastdbcmd, blastp; 
       L<http://www.ncbi.nlm.nih.gov/books/NBK279671/>
      
       =item The Pfam-A.hmm database
      
       Available from
       L<ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/>.  This has
       to be formatted for HMMER searches with the UNIX command: C<hmmpress
       Pfam-A.hmm>
      
       =item Perl packages, including BioPython packages
       These are best installed with C<cpanm>
      
        =over
      
        =item
        List::MoreUtils
      
        =item
        LWP::Simple
      
        =item
        Bio::SeqIO
      
        =item
        List::MoreUtils
      
        =item
        Statistics::R
      
        =item
        Text::Markdown
      
        =back
      
       =item R
       Available from L<https://cran.r-project.org/>
      
       =item R packages
      
        =over
      
        =item
        scales
      
        =item
        seqinr
      
        =back
       
       =item  MAFFT
       Available from L<http://mafft.cbrc.jp/alignment/software/>

 Preparations

  Transcriptome

    The transcriptome must be prepared as a database with the following
    UNIX command:

     makeblastdb -in fastafilename.fasta -title fastafilename -out -parse-seqids -dbtype nucl

    Here, fastafilename needs to be changed to the name of your own
    fastafile.

  Protein database

    A database with protein sequences from species closely related to your
    target species must be prepared as follows:

    Download the nr.gz database from ftp://ftp.ncbi.nih.gov/blast/db/FASTA/

    unpack it with the Unix code

       gunzip nr.gz

    format it as blast database on the command line with

       makeblastdb -dbtype prot -in nr -parse_seqids

    Search the Entrez Protein database
    (http://www.ncbi.nlm.nih.gov/protein) for the wider taxon of your
    target species (e.g. bony fishes or arthropoda) with the query: "bony
    fishes"[ORGN] or "arthropoda"[ORGN]

    Select 'Send to File', choose the format 'GI list' and save it in a
    file sequence.gi.txt.

    Now, run the Unix command

        blastdb_aliastool -gilist sequence.gi.txt -db nt -out nt_bonyfishes -title nt_bonyfishes

      Change bonyfishes to the wider taxon name of your own target species.
      For example, if your tanscriptome belongs to a decapod species, you
      can choose arthropoda as wider taxon.

  ESTs

    Download from the NCBI EST database (http://www.ncbi.nlm.nih.gov/est)
    all entries for your target species as fasta file and format it as
    blast database with the command

     makeblastdb -in fastafilename.fasta -dbtype nucl -parse_seqids 

    Here, fastafilename needs to be replaced with the name of your own
    fasta file containing EST sequences.

 Running the main function

    The ProteinScreen module provides the function
    ProteinScreen::ProteinScreen. This function runs a pipeline to extract
    from a transcriptome those contigs that are closely related to
    previously selected target proteins, e.g. heat shock proteins. The user
    needs to provide five arguments:

    $queries

      Name of a fasta file with peptide sequences of target proteins from
      related species, including the file ending '.fasta'.

    $pfam

      Link to the protein database pfam, downloaded from
      http://pfam.xfam.org/ as Pfam-A.hmm file.

    $transcriptome

      Fasta file name of the transcriptome, including the ending '.fasta'.

    $proteindatabase

      Name of the local non-redundant protein database of the wider taxon
      of your target species.

    $ests

      Name of the fasta file containing EST sequences of your target
      species, including the ending '.fasta'.

  Output

    The wrapper function ProteinScreen::ProteinScreen summarizes all
    results in the file ResultReport.html. This file provides a link to
    'Queries', which lists all protein queries that were used, and a link
    to 'Best hits', which lists all the contigs in your target species'
    transcriptome that are closely related to one or several of your
    protein queries. 'Best hits' lists for each contig best hits in the
    protein database, as well as links to predicted domains, used protein
    queries, and related Expressed Sequence Tags (ESTs).

  Pipeline details

    The wrapper function ProteinScreen::ProteinScreen executes
    automatically the following functions in the following order:

    local_database_search

      Compares the fasta file of protein sequences (first argument) against
      a transriptome database (second argument) using the tblastn function
      of NCBI's blast+ package.

    parse_local_database_search

      The single argument of this function is a list of tblastn report
      files resulting from the previous function. This function extracts
      from each of the report files only the best hits (E-value < 1e-05) of
      protein queries in the transcriptome databse.

    best_hits_to_fasta

      Extracts sequences of the unique best contigs hits (first argument,
      output of the previous function) from the transcriptome database
      (second argument) and saves them in a single fasta file.

    translation

      A fasta file of nucleotide sequences (first argument) listed as best
      contig hits (second argument) is translated to protein sequences,
      which are saved in a new fasta file. The reading frame is based on
      the protein blast hits that were detected with the
      local_database_search function.

    deduplicate_fasta

      The single argument to this function is the fasta file with protein
      sequences that was put out by the previous function. This function
      removes potential duplicates of transcriptome contigs that are
      related to query proteins. Duplicates occur when one contig showed up
      as best blast hit to different protein queries.

    parse_decypher2

      XX Continue here - test the code

      XX I think that this function is not necessary at all. The output
      (UniqueComphits.fasta) is the same as deduplicated.fasta?! XX Try
      running the script when removing this function

    peptide_extraction

      Extracts from the fasta file (the function's single argument) of
      translated contigs, related to protein queries, the longest peptide
      sequence that starts with a methionine and ends with a stop codon.
      The longest peptides of each contig are saved to a new fasta file

    reciprocal_blast

      Peptides extracted with the previous function (first argument) are
      blasted against a protein database (second argument) of taxa related
      to your target species (whose transcriptome you are exploring). The
      function uses the blastp function of NCBI's BLAST+ package.

    parse_blastp

      Removes those blast hits from the output files of the previous
      function (first argument) that refer to 'puative', 'predicted',
      'provisional', 'partial', 'unknown', or 'hypothetic' proteins. The
      sequences of the remaining hits are extracted from a protein database
      of taxa related to your target species (second argument) and saved to
      a fasta file.

    hmmscan

      This function runs the program hmmscan to identify protein domains in
      the translated contigs of your target species (fasta file as first
      argument) and in the corresponding best blast hits from the previous
      function (fasta file as second argument). Hmmscan searches the
      protein sequences against Pfam, a large database of protein families
      (third argument).

    mafft

      The three arguments to this function are: 1) a fasta file with
      translated contigs of your target species that wer used as queries ,
      2) a fasta file with protein sequences of best hits in related
      species, and 3) a file informing on which sequences in these two
      fasta files form a pair. These pairs of protein sequences are then
      aligned against each other with the program mafft. The alignments are
      saved in both text (ClustalW alignment) and fasta files.

    parse_mafft

      This function goes through a list of ClustalW alignment files
      (argument of the function) and calculates the percent identity and
      percent similarity of the aligned sequences. The percent identity is
      calculated as the number of amino acids that match exactly (indicated
      by '*' in the alignment file) divided by the total number of aligned
      amino acids. The percent similarity is calculated as the number of
      amino acids that match exactly (indicated by '*' ), or show similar
      properties (indicated by '.' and ':' in the alignment file) divided
      by the total number of aligned amino acids.

    vetting

      uses the program tblastn to blast translated contigs of your target
      species (first argument that was saved in a fasta file by the
      function <peptide_extraction>) against a fasta file of EST sequences
      of your target species (second argument of the function). The file
      named 'tblastn.ESTout' which lists all the produced alignment files.
      The alignment files themselves end with ...tblastnresult.out

    parse_vetting

      Saves a table called Vetting.out in which it stores important
      information fextracted from the tblastn.ESTout files produced by the
      preceding function vetting.

    pfam

      Uses an R script to plot protein domains along the alignment of the
      translated contigs of your target species and their best hits.

    create_report

      This function creates six html files that are all accessible via
      links in the single ResultReport.html file. This makes all output
      information of all the previous functions are accessible to the user.

LICENSE

    This is released under the Artistic License. See perlartistic.

AUTHOR

    Alexander Jueterbock - http://marinetics.org/

