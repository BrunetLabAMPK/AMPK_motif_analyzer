AMPK Motif Analyzer
===================

Uses a position weight matrix to calculate a distance score for input motifs based on how closely they match the AMPK motif. Position weight matrix was generated by taking the log10 of the (frequency of each amino acid from N terminal 5 to C terminal 4 in the AMPK motif of 109 known AMPK substrates divided by a background frequency calculated by taking the average amino acid occurrence at each matching position in 10,000 matched background datasets containing human phosphoserines and threonines). The higher the score, the closer the motif to the AMPK motif.

Lowest motif score allowed for high confidence: 1.016513 (the mean minus one standard deviation of the 109 known AMPK substrates when scored with this algorithm).

Last update: Decebmer 16, 2014

Installation
------------

Download files from GITHUB on local computer. You just need Perl to be installed on your computer, and you are good to go.

Usage
------------

Run Perl from terminal like this

    perl Rank_Motifs.pl --frequency <filename or path> --input <filename or path> --output <filename or path>

The command line options are:
* --frequency  : The AMPK motif file provided with this package (AMPK_motif_109_standardized_log10.txt)
* --input      : The phosphorylation motifs to be scored. Motif to be scored MUST be in the first column with any attached data important to the user in additional columns.
* --output     : Name of the output file with scored motifs

AN example input and output file is provided with this package. Try to run

   perl Rank_Motifs.pl --frequency AMPK_motif_109_standardized_log10.txt --input Example_input.txt --output Example_input_scored.txt
   
