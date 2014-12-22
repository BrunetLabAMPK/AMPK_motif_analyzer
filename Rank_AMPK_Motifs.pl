#!/usr/bin/perl
#
#     Rank_AMPK_Motifs.pl  (v1)			Last Modified: 2014-12-17		© Bethany Schaffer
#
# ================ DESCRIPTION =====================
#
# This script scores any given phosphorylation motifs based on how well they match the AMPK phosphorylation 
# motif. The higher the score, the better. Each motif of interest is scored by summing the log10 of the 
# standardized AMPK motif frequencies (see "frequency file" below).
#
# ================= HOW TO RUN =====================

# Required files:
#
# (1) frequency file:  AMPK_motif_109_standardized_log10.txt contains the log10 of (frequencies of amino 
#                      acids in the AMPK motifs of 109 AMPK phosphorylation sites from N5 to C4 divided 
#                      by the average frequency of that amino acid at that location in 10,000 matched 
#                      background datasets of human phosphorylation sites
# (2) Input file:      Any tab delimited text file of the user's choosing that contains the phosphorylation 
#                      motif(s) of interest in the first column. Rest of the columns will be ignored.
# Output:
# A file, with the motif score added in the final column. User must provide the name of this file on terminal.
#
# To run type the following command in the termninal:
# perl Rank_AMPK_Motifs.pl --frequency AMPK_motif_109_standardized_log10.txt --input Example_input.txt --output Example_input_scored.txt
#
# --frequency  is the AMPK motif file provided with the script
# --input      is the file having motifs to be scored
# --output     is the output file with scored motifs
#
# ================= CONTACT =====================
#
# Direct any questions, suggestions or bugs to Bethany Schaffer (beschaff@stanford.edu)
#
# ================= Reference =====================
#
# Schaffer et al. 
# =================================================

use strict;
use warnings;
use Getopt::Long;

# Variables to hold file names
my $freq = "";
my $input = "";
my $output = "";

# Get options based on command line input
GetOptions ("frequency=s" => \$freq,    # file name is string
	          "input=s"   => \$input,      # file name is string
              "output=s"  => \$output)   # file name is string
or die("Error in command line arguments\n");

# Check if all the input are fine

if ($freq eq ""){die "provide data frequency file\n";}
if ($input eq ""){die "provide data input file\n";}
if ($output eq ""){die "provide data output file\n";}

# Open the file containing the log10 of the standardized AMPK frequencies
open AMPK, $freq || die;

my %subsbacklog10 = (); # Hash to store the AMPK motif information

# Add each frequency to %subsbacklog10
while (<AMPK>)
{
	chomp;
	my @ampk = split/\t/;	
	$subsbacklog10{$ampk[0]}{$ampk[1]} = $ampk[2]; 	# Add information to hash such that the position (ex neg5) and amino acid (ex G) are used to call the assocaited frequency
}
close AMPK;

open INPUT, $input || die "couldn't find the input file : $!"; # Open the user's input file
open (OUTPUT, '>', $output) or die "couldn't create new file : $!"; # Create output file

my %score_hash = (); # Hash to store the motifs, associated information, and scores


# For each line of input in the user's file, execute the following loop to score motifs
while (<INPUT>)
{
	chomp;
	my @data = split/\t/;	
	
	my $phosphomotif = $data[0]; 	# Define the phosphorylation motif, which MUST be in the first column
			
	# Create a variable to hold the AA at each site
	my $neg5 = substr $phosphomotif, 0, 1;
	my $neg4 = substr $phosphomotif, 1, 1;
	my $neg3 = substr $phosphomotif, 2, 1;
	my $neg2 = substr $phosphomotif, 3, 1;
	my $neg1 = substr $phosphomotif, 4, 1;
	my $site = substr $phosphomotif, 5, 1;
	my $pos1 = substr $phosphomotif, 6, 1;
	my $pos2 = substr $phosphomotif, 7, 1;
	my $pos3 = substr $phosphomotif, 8, 1;
	my $pos4 = substr $phosphomotif, 9, 1;
	
	# Only score the phosphoserine and threonines (phosphotryosines will be ignored)
	if(($site eq "S") || ($site eq "T"))
	{
		# Only score motifs that contain a basic residue within 5 N terminal positions of the phosphorylated residue
		if(($neg5 eq "R") || ($neg4 eq "R") || ($neg3 eq "R") || ($neg2 eq "R") || ($neg1 eq "R") || ($neg5 eq "K") || ($neg4 eq "K") || ($neg3 eq "K") || ($neg2 eq "K") || ($neg1 eq "K") || ($neg5 eq "H") || ($neg4 eq "H") || ($neg3 eq "H") || ($neg2 eq "H") || ($neg1 eq "H"))
		{
		
			my %motif = (); # Another hash to store each amino acid in the user's motif of interest
			
			# Assign amino acids to hash based on their location in the user's motif of interest
			$motif{neg5} = $neg5;
			$motif{neg4} = $neg4;
			$motif{neg3} = $neg3;
			$motif{neg2} = $neg2;
			$motif{neg1} = $neg1;
			$motif{pos1} = $pos1;
			$motif{pos2} = $pos2;
			$motif{pos3} = $pos3;
			$motif{pos4} = $pos4;
			
			# Set the initial distance score to 0
			my $distance = 0;
			
			# Now score the motifs
			# Start by looking at each site in the motif from N5 to C4; the phosphoserine and threonines are not incorporated into the score
			foreach my $sites (keys %motif)
			{
				# Open the AMPK motif hash
				foreach my $location (keys %subsbacklog10)
				{
					# For the key in the AMPK motif ($location) that matches the location in the input motif ($sites)
					if($location eq $sites)
					{
						# Go through the amino acids in the AMPK frequency hash to find the one that matches the input motif
						foreach my $AA (keys %{$subsbacklog10{$location}})
						{
							# Enter the loop only when the amino acid at the location of interest matches for the AMPK motif and the input motif
							if($AA eq $motif{$sites})
							{
								# Calculate the distance score based on the AMPK motif
								if($distance == 0) 	# If the distance is 0, the new distance score is equal to the log10 frequency in the AMPK motif for the amino acid of interest
								{
									$distance = $subsbacklog10{$location}{$AA};
								}
								else # If the distance score has already been calculated for some residues, then add the new frequency to the preexisting score
								{
									$distance = $distance + $subsbacklog10{$location}{$AA};
								}					
							}
						}
					}
				}
			}
			# Add the distance score to the line of data from the user's input file
			push (@data, $distance);	
			push (@{$score_hash{$distance}}, \@data); # Add the data line to a hash where the distance score is the key, and the value is an array of arrays -- to account for possible identical scores if two motifs are the same
			#$score_hash{$distance} = \@data;
			#print OUTPUT join("\t", @{$score_hash{$distance}}), "\n";
		}	
	}
}
close INPUT;


# Sort the hash containing the user's input data + the distance score such that the output file is rank ordered i.e. highest/best scoring motif first
foreach (sort {$b <=> $a} keys %score_hash)
{
	# Print out each dataline associated with a given distance score (most scores will be unique to one array)
	foreach my $line (@{$score_hash{$_}})
	{
		print OUTPUT join("\t", @$line), "\n";

	}
}
close OUTPUT;

