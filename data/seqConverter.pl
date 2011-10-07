#!/usr/bin/perl -w
#
# seqConverter.pl v1.2
# Last modified August 12, 2010 21:09
# (c) Olaf R.P. Bininda-Emonds
#
# Input:
#   Sequence data in any of fasta, GenBank, nexus, (classic or extended) phylip,
#	or Se-Al formats
#
# Output:
#   Sequence data in any of fasta, nexus, (classic or extended) phylip, and/or
#	Se-Al formats.
#
# Usage: seqConverter.pl -d<filename> -o<f|n|pc|pe|s> [-a] [-c<number>] [-g]
#		 [-G] [-H] [-i<f|g|n|p|s>] [-j] [-l<number>] [-n] [-r<a|i>] [-s] [-t]
#		 [-u] [-v] [-h]
#	options: -a = print out accession numbers instead of taxon labels for nexus
#				  and phylip output
#			 -c<number> = global genetic code (to be used for translation unless
#						  otherwise specified) (default = standard code; end
#						  with * to override local settings)
#			 -d<filename> = file containing raw sequence information; * = batch
#							convert all specified file types in working
#							directory (suffixes must be .fasta / .fas, .gb,
#							.nexus / .nex, .phylip / .ph, or .seal as
#							appropriate)
#			 -g = do not convert ~ gap characters (e.g., from BioEdit) to -
#			 -G = convert flanking gaps to Ns
#			 -H = convert sequence data to haplotype data
#            -i<f|g|n|p|s> = format of sequence file (fasta (f), GenBank (g),
#							 nexus (n), phylip (p), or Se-Al (s)); default =
#							 autodetect
#            -j = produce jackknifed data sets each with a single sequence
#				  deleted
#            -l<number> = interleave nexus- and/or phylip-formatted output
#						  (default = sequential) with specified sequence length
#						  (between 10 and 100 inclusive; default = 80);
#						  NOTE: also sets default interleave length for fasta
#						        output
#			 -n = convert ambiguous nucleotides to Ns
#            -o<f|n|pc|pe|s> = output results in any or all of fasta (f), nexus
#							   (n), classic or extended phylip (pc or pe),
#							   and/or Se-Al (s) formats (NOTE: flag can be
#							   invoked multiple times)
#			 -r<a|i> = order sequences in final output alphabetically by name
#					   (a; default) or in input order from file (i)
#			 -s = split data set into individual partitions according to nexus
#				  charset statements
#            -t = translate sequences into amino acids
#            -u = interactive user-input mode
#            -h = print this message and quit
#            -v = verbose output

use strict;
use POSIX;

# Set user-input defaults and associated parameters
	# Data set variables
		my $inputType = "";	# Options are "fasta", "GenBank", "nexus", "phylip", and "Se-Al"
		my ($readFile, @seqFiles);
			my $dataSource;
 		my $seqType = "nucleotide";
		my $globalGenCode = 0;
			my $globalOverride = 0;

		my (@accNum, %nameLabel, %sequence, %geneticCode, %accPresent);
		my (@charsetList, %charsetStart, %charsetEnd);
		my (%deletedSeq, %finalSeq);
		my $seqCount;
		my $ntax;
		my $nameClean = 1;
		my $skipDuplicates = 1;

	# User input variables
 		my $seqOrder = "alphabetical";	# Options are "alphabetical" (default) and "input"
 		my $seqSplit = 0;
 		my $jackknife = 0;
 		my $gapChange = 1;
 		my $ambigChange = 0;
 		my $translateSeqs = 0;
		my $haploTyping = 0;
			my $haploFile;
		my $flankGap = 0;

	# Translation variables and genetic codes
		my %transTable = ('1' => 'standard',
						  '2' => 'vertebrate mitochondrial',
						  '3' => 'yeast mitochondrial',
						  '4' => 'mold, protozoan and colenterate mitochondrial and  mycoplasam/spiroplasma',
						  '5' => 'invertebrate mitochondrial',
						  '6' => 'ciliate, dasycladacean and hexamita nuclear',
						  '9' => 'echinoderm mitochondrial',
						  '10' => 'euplotid nuclear',
						  '11' => 'bacterial and plant plastid',
						  '12' => 'alternative yeast nuclear',
						  '13' => 'ascidian mitochondrial',
						  '14' => 'alternative flatworm mitochondrial',
						  '15' => 'Blepharisma nuclear',
						  '16' => 'chlorophycean mitochondrial',
						  '21' => 'trematode mitochondrial',
						  '22' => 'Scenedesmus obliquus mitochondrial',
						  '23' => 'Thraustochytrium mitochondrial');
		my %gb2seal = ('1' => '0',
					   '2' => '1',
					   '3' => '2',
					   '4' => '3',
					   '4' => '4',
					   '5' => '5',
					   '6' => '6',
					   '9' => '7',
					   '10' => '8',
					   '11' => '9',
					   '12' => '10',
					   '13' => '11',
					   '14' => '12',
					   '15' => '13',
					   '16' => '1',
					   '21' => '1',
					   '22' => '1',
					   '23' => '1');
		my %seal2gb = ('0' => '1',
					   '1' => '2',
					   '2' => '3',
					   '3' => '4',
					   '4' => '4',
					   '5' => '5',
					   '6' => '6',
					   '7' => '9',
					   '8' => '10',
					   '9' => '11',
					   '10' => '12',
					   '11' => '13',
					   '12' => '14',
					   '13' => '15');
		my %genCodes;
			foreach my $code (qw(0 1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
				{
				$genCodes{$code} = 1;
				}
		my %DNAtoAA;
		my @ambigList = qw(A C G T M R W S Y K V H D B N);
		my %ambigCode = ('A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T',
						'AC' => 'M', 'AG' => 'R', 'AT' => 'W', 'CG' => 'S', 'CT' => 'Y', 'GT' => 'K',
						'ACG' => 'V', 'ACT' => 'H', 'AGT' => 'D', 'CGT' => 'B', 'ACGT' => 'N');
		my (%constitNTlist);
			while (my ($nt, $code) = each %ambigCode)	# Where $nt = key and $code = value
				{
				push @{$constitNTlist{$code}}, $_ foreach (split("",$nt));
				}

	# Output variables
		my $maxLength;
		my $fastaPrint = 0;
			my $fastaOut;
		my $nexusPrint = 0;
			my $nexusOut;
		my ($phylipTradPrint, $phylipExtPrint) = (0, 0);
			my $phylipOut;
		my $sealPrint = 0;
			my $sealOut;
		my $outFormat = "sequential";
			my $interleaveLength = 80;
			my $fastaLength = 80;
		my $accPrint = 0;

	# Miscellaneous variables
		my $verbose = 0;
		my $debug = 0;
		my $perlScript = "seqConverter.pl";
		my $version = "1.2";

# Read in user input
	if (not @ARGV or join(' ', @ARGV) =~ /\s-u/ or $ARGV[0] =~ /^-u/)	# Enter interactive user-input mode
		{
		print "Entering interactive user-input mode. Type \"q\" at any prompt to exit program.\n";

		# Get print format
			my $defaultAcc = ($accPrint) ? "y" : "n";
				undef $accPrint;
			until (defined $accPrint)
				{
				print "\tPrint out accession numbers only to nexus- and/or phylip output (y|n) [$defaultAcc]: ";
				$accPrint = <stdin>;
					chomp ($accPrint);
					exit(0) if ($accPrint eq "q");
				if (substr($accPrint, 0, 1) =~ /^n/i or $accPrint eq "")
					{
					$accPrint = 0;
					}
				elsif (substr($accPrint, 0, 1) =~ /^y/i)
					{
					$accPrint = 1;
					}
				else
					{
					print "\t\tInvalid input ($accPrint)\n";
					undef $accPrint;
					}
				}

		# Get datafile
			until (defined $readFile)
				{
				print "\tEnter name of data file (* = batch convert): ";
				$readFile = <stdin>;
					chomp ($readFile);
					exit(0) if ($readFile eq "q");
				unless (-e $readFile or $readFile eq "*")
					{
					print "\t\tFile '$readFile' does not exist\n";
					undef $readFile;
					}
				}
			push @seqFiles, $readFile;
		
		# Get format of datafile
			my $defaultInput = "autodetect";
				undef $inputType;
			until (defined $inputType)
				{
				print "\tEnter format of file $readFile (fasta|GenBank|nexus|phylip|Se-Al) [$defaultInput]: ";
				$inputType = <stdin>;
					chomp ($inputType);
					exit(0) if ($inputType =~ /^q/i);
					if (substr($inputType, 0, 1) =~ /^a/i or $inputType eq "")
						{
						$inputType = "autodetect";
						}
					elsif (substr($inputType, 0, 1) =~ /^f/i)
						{
						$inputType = "fasta";
						}
					elsif (substr($inputType, 0, 1) =~ /^g/i)
						{
						$inputType = "GenBank";
						}
					elsif (substr($inputType, 0, 1) =~ /^n/i)
						{
						$inputType = "nexus";
						}
					elsif (substr($inputType, 0, 1) =~ /^p/i)
						{
						$inputType = "phylip";
						}
					elsif (substr($inputType, 0, 1) =~ /^s/i)
						{
						$inputType = "Se-Al";
						}
					else
						{
						print "\t\tInvalid input ($inputType)\n";
						undef $inputType;
						}
				}
			$inputType = "" if ($inputType eq "autodetect");

		# Get whether or not to clean sequence labels
			my $defaultClean = ($nameClean) ? "y" : "n";
				undef $nameClean;
			until (defined $nameClean)
				{
				print "\tClean sequence labels by changing non-alphanumeric characters to underscores (y|n)? [$defaultClean]: ";
				$nameClean = <stdin>;
					chomp ($nameClean);
					exit(0) if ($nameClean =~ /^q/i);
					if (substr($nameClean, 0, 1) =~ /^y/i or $nameClean eq "")
						{
						$nameClean = 1;
						}
					elsif (substr($nameClean, 0, 1) =~ /^n/i)
						{
						$nameClean = 0;
						}
					else
						{
						print "\t\tInvalid input ($seqSplit)\n";
						undef $nameClean;
						}
				}
		
		# Get whether or not to split sequences
			my $defaultSplit = ($seqSplit) ? "y" : "n";
				undef $seqSplit;
			until (defined $seqSplit)
				{
				print "\tSplit into individual partitions following charset statements (y|n)? [$defaultSplit]: ";
				$seqSplit = <stdin>;
					chomp ($seqSplit);
					exit(0) if ($seqSplit =~ /^q/i);
					if (substr($seqSplit, 0, 1) =~ /^n/i or $seqSplit eq "")
						{
						$seqSplit = 0;
						}
					elsif (substr($seqSplit, 0, 1) =~ /^y/i)
						{
						$seqSplit = 1;
						}
					else
						{
						print "\t\tInvalid input ($seqSplit)\n";
						undef $seqSplit;
						}
				}
		
		# Get whether or not to jackknife sequences
			my $defaultJack = ($jackknife) ? "y" : "n";
				undef $jackknife;
			until (defined $jackknife)
				{
				print "\tProduce replicate, jackknifed data sets (y|n)? [$defaultJack]: ";
				$jackknife = <stdin>;
					chomp ($jackknife);
					exit(0) if ($jackknife =~ /^q/i);
					if (substr($jackknife, 0, 1) =~ /^n/i or $jackknife eq "")
						{
						$jackknife = 0;
						}
					elsif (substr($jackknife, 0, 1) =~ /^y/i)
						{
						$jackknife = 1;
						}
					else
						{
						print "\t\tInvalid input ($jackknife)\n";
						undef $jackknife;
						}
				}

		# Get whether or not to convert gaps
			my $defaultGaps = ($gapChange) ? "y" : "n";
				undef $gapChange;
			until (defined $gapChange)
				{
				print "\tConvert ~ gap characters to - (y|n)? [$defaultGaps]: ";
				$gapChange = <stdin>;
					chomp ($gapChange);
					exit(0) if ($gapChange =~ /^q/i);
					if (substr($gapChange, 0, 1) =~ /^n/i)
						{
						$gapChange = 0;
						}
					elsif (substr($gapChange, 0, 1) =~ /^y/i or $gapChange eq "")
						{
						$gapChange = 1;
						}
					else
						{
						print "\t\tInvalid input ($gapChange)\n";
						undef $gapChange;
						}
				}
		
		# Get whether or not to convert flanking gaps
			my $defaultFlank = ($flankGap) ? "y" : "n";
				undef $flankGap;
			until (defined $flankGap)
				{
				print "\tConvert flanking gap characters to Ns (y|n)? [$defaultFlank]: ";
				$flankGap = <stdin>;
					chomp ($flankGap);
					exit(0) if ($flankGap =~ /^q/i);
					if (substr($flankGap, 0, 1) =~ /^n/i)
						{
						$flankGap = 0;
						}
					elsif (substr($flankGap, 0, 1) =~ /^y/i or $flankGap eq "")
						{
						$flankGap = 1;
						}
					else
						{
						print "\t\tInvalid input ($flankGap)\n";
						undef $flankGap;
						}
				}
		
		# Get whether or not to convert ambiguous nucleotides
			my $defaultAmbig = ($ambigChange) ? "y" : "n";
				undef $ambigChange;
			until (defined $ambigChange)
				{
				print "\tConvert ambiguous nucleotides to Ns (y|n)? [$defaultAmbig]: ";
				$ambigChange = <stdin>;
					chomp ($ambigChange);
					exit(0) if ($ambigChange =~ /^q/i);
					if (substr($ambigChange, 0, 1) =~ /^n/i or $ambigChange eq "")
						{
						$ambigChange = 0;
						}
					elsif (substr($ambigChange, 0, 1) =~ /^y/i)
						{
						$ambigChange = 1;
						}
					else
						{
						print "\t\tInvalid input ($ambigChange)\n";
						undef $ambigChange;
						}
				}
		
		# Get output order of sequences
			my $defaultOrder = $seqOrder;
				undef $seqOrder;
			until (defined $seqOrder)
				{
				print "\tEnter output order for sequences (alphabetical|clustal|input file) [$defaultOrder]: ";
				$seqOrder = <stdin>;
					chomp ($seqOrder);
					exit(0) if ($seqOrder =~ /^q/i);
					if (substr($seqOrder, 0, 1) =~ /^i/i)
						{
						$seqOrder = "input";
						}
					elsif (substr($seqOrder, 0, 1) =~ /^a/i or $seqOrder eq "")
						{
						$seqOrder = "alphabetical";
						}
					else
						{
						print "\t\tInvalid input ($seqOrder)\n";
						undef $seqOrder;
						}
				}

		# Get whether or not to convert to haplotypes
			my $defaultHaplo = ($haploTyping) ? "y" : "n";
				undef $haploTyping;
			until (defined $haploTyping)
				{
				print "\tConvert sequence data to haplotypes (y|n)? [$defaultHaplo]: ";
				$haploTyping = <stdin>;
					chomp ($haploTyping);
					exit(0) if ($haploTyping =~ /^q/i);
					if (substr($haploTyping, 0, 1) =~ /^n/i or $haploTyping eq "")
						{
						$haploTyping = 0;
						}
					elsif (substr($haploTyping, 0, 1) =~ /^y/i)
						{
						$haploTyping = 1;
						}
					else
						{
						print "\t\tInvalid input ($haploTyping)\n";
						undef $haploTyping;
						}
				}

		# Get whether or not to convert to amino acids
			my $defaultTranslate = ($translateSeqs) ? "y" : "n";
				undef $translateSeqs;
			until (defined $translateSeqs)
				{
				print "\tConvert sequence data to amino acids (y|n)? [$defaultTranslate]: ";
				$translateSeqs = <stdin>;
					chomp ($translateSeqs);
					exit(0) if ($translateSeqs =~ /^q/i);
					if (substr($translateSeqs, 0, 1) =~ /^n/i or $translateSeqs eq "")
						{
						$translateSeqs = 0;
						}
					elsif (substr($translateSeqs, 0, 1) =~ /^y/i)
						{
						$translateSeqs = 1;
						}
					else
						{
						print "\t\tInvalid input ($translateSeqs)\n";
						undef $translateSeqs;
						}
				}

		# Get genetic code for translation
			my $defaultCode = $globalGenCode;
				$globalGenCode = 99;
			print "\tGenetic codes available for protein translation:\n";
				foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
					{
					print "\t\t$code: $transTable{$code}";
					print "\n";
					}
			until (defined $genCodes{$globalGenCode})
				{
				$globalOverride = 0;
				print "\tEnter global genetic code to be applied (0 = no translation; follow with an * to override local code) [$defaultCode]: ";
				$globalGenCode = <stdin>;
					chomp ($globalGenCode);
					exit(0) if ($globalGenCode =~ /^q/i);
					$globalGenCode = $defaultCode if ($globalGenCode eq "");
					if ($globalGenCode =~ /\*$/)
						{
						$globalOverride = 1;
						$globalGenCode =~ s/\*$//;
						}
					print "\t\tInvalid input ($globalGenCode)\n" unless (defined $genCodes{$globalGenCode} or $globalGenCode == 0);
				}
				
		# Get output formats
			my $defaultFasta = ($fastaPrint) ? "y" : "n";
				undef $fastaPrint;
			until (defined $fastaPrint)
				{
				print "\tOutput results in fasta format (y|n) [$defaultFasta]: ";
				$fastaPrint = <stdin>;
					chomp ($fastaPrint);
					exit(0) if ($fastaPrint =~ /^q/i);
					if (substr($fastaPrint, 0, 1) =~ /^y/i)
						{
						$fastaPrint = 1;
						}
					elsif (substr($fastaPrint, 0, 1) =~ /^n/i or $fastaPrint eq "")
						{
						$fastaPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($fastaPrint)\n";
						undef $fastaPrint;
						}
				}

			my $defaultNexus = ($nexusPrint) ? "y" : "n";
				undef $nexusPrint;
			until (defined $nexusPrint)
				{
				print "\tOutput results in nexus format (y|n) [$defaultNexus]: ";
				$nexusPrint = <stdin>;
					chomp ($nexusPrint);
					exit(0) if ($nexusPrint =~ /^q/i);
					if (substr($nexusPrint, 0, 1) =~ /^y/i)
						{
						$nexusPrint = 1;
						}
					elsif (substr($nexusPrint, 0, 1) =~ /^n/i or $nexusPrint eq "")
						{
						$nexusPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($nexusPrint)\n";
						undef $nexusPrint;
						}
				}

			my $defaultPhylip = ($phylipTradPrint) ? "y" : "n";
				undef $phylipTradPrint;
			until (defined $phylipTradPrint or $phylipExtPrint)
				{
				print "\tOutput results in traditional phylip format (y|n) [$defaultPhylip]: ";
				$phylipTradPrint = <stdin>;
					chomp ($phylipTradPrint);
					exit(0) if ($phylipTradPrint =~ /^q/i);
					if (substr($phylipTradPrint, 0, 1) =~ /^y/i)
						{
						$phylipTradPrint = 1;
						}
					elsif (substr($phylipTradPrint, 0, 1) =~ /^n/i or $phylipTradPrint eq "")
						{
						$phylipTradPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($phylipTradPrint)\n";
						undef $phylipTradPrint;
						}
				}
				
				if ($phylipTradPrint == 0)	# Check for extended format
					{
					my $defaultPhylip = ($phylipExtPrint) ? "y" : "n";
						undef $phylipExtPrint;
					until (defined $phylipExtPrint or $phylipExtPrint)
						{
						print "\tOutput results in extended phylip format (y|n) [$defaultPhylip]: ";
						$phylipExtPrint = <stdin>;
							chomp ($phylipExtPrint);
							exit(0) if ($phylipExtPrint =~ /^q/i);
							if (substr($phylipExtPrint, 0, 1) =~ /^y/i)
								{
								$phylipExtPrint = 1;
								}
							elsif (substr($phylipExtPrint, 0, 1) =~ /^n/i or $phylipExtPrint eq "")
								{
								$phylipExtPrint = 0;
								}
							else
								{
								print "\t\tInvalid input ($phylipExtPrint)\n";
								undef $phylipExtPrint;
								}
						}
					}

			my $defaultSeal = ($sealPrint) ? "y" : "n";
				undef $sealPrint;
			until (defined $sealPrint)
				{
				print "\tOutput results in Se-Al format (y|n) [$defaultSeal]: ";
				$sealPrint = <stdin>;
					chomp ($sealPrint);
					exit(0) if ($sealPrint =~ /^q/i);
					if (substr($sealPrint, 0, 1) =~ /^y/i)
						{
						$sealPrint = 1;
						}
					elsif (substr($sealPrint, 0, 1) =~ /^n/i or $sealPrint eq "")
						{
						$sealPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($sealPrint)\n";
						undef $sealPrint;
						}
				}
				
			if ($nexusPrint or $phylipTradPrint or $phylipExtPrint)
				{
				my $defaultInterleave = ($outFormat) ? "y" : "n";
					undef $outFormat;
				until (defined $outFormat)
					{
					print "\tInterleave output for nexus- and/or phylip-formatted output with length (n|10-100) [$defaultInterleave]: ";
					$outFormat = <stdin>;
						chomp ($outFormat);
						exit(0) if ($outFormat =~ /^q/i);
						if ($outFormat =~ /(\d+)/)
							{
							if ($outFormat =~ /\D/ or $outFormat < 10 or $outFormat > 100)
								{
								print "\t\tInvalid input ($outFormat)\n";
								undef $outFormat;
								}
							else
								{
								$outFormat = $1;
								}
							}
						elsif (substr($outFormat, 0, 1) =~ /^n/i or $outFormat eq "")
							{
							$outFormat = "sequential";
							}
						else
							{
							print "\t\tInvalid input ($outFormat)\n";
							undef $outFormat;
							}
					}
				if ($outFormat ne "sequential")
					{
					$interleaveLength = $outFormat;
					$outFormat = "interleave";
					}
				}

		# Get verbose output mode
			my $defaultVerbose = ($verbose) ? "y" : "n";
				undef $verbose;
			until (defined $verbose)
				{
				print "\tOutput verbose results to screen (y|n) [$defaultVerbose]: ";
				$verbose = <stdin>;
					chomp ($verbose);
					exit(0) if ($verbose =~ /^q/i);
					if (substr($verbose, 0, 1) =~ /^y/i)
						{
						$verbose = 1;
						print "\n";
						}
					elsif (substr($verbose, 0, 1) =~ /^n/i or $verbose eq "")
						{
						$verbose = 0;
						}
					elsif (substr($verbose, 0, 1) =~ /^x/i or $verbose eq "")
						{
						$verbose = $debug = 1;
						}
					else
						{
						print "\t\tInvalid input ($verbose)\n";
						undef $verbose;
						}
				}
		}
	elsif (join(' ', @ARGV) =~ /\s-h/ or $ARGV[0] =~ /^-h/)	# Print help screen
		{
		print "Usage: seqConverter.pl -d<filename> -o<f|n|pc|pe|s> [-a] [-c<number>] [-g]\n";
		print "       [-G] [-H] [-i<f|g|n|p|s>] [-j] [-l<number>] [-n] [-r<a|i>] [-s] [-t]\n";
		print "       [-u] [-v] [-h]\n";
		print "Version: $version\n";
		print "Options: -a = print out accession numbers instead of taxon labels for nexus\n";
		print "              and phylip output\n";
		print "         -c<number> = global genetic code (to be used for translation unless\n";
		print "                      otherwise specified) (default = standard code; end\n";
		print "                      with * to override local settings)\n";
			foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
				{
				print "                      $code: $transTable{$code}";
				print "\n";
				}
		print "         -d<filename> = file containing raw sequence information; * = batch\n";
		print "                        convert all specified file types in working\n";
		print "                        directory (suffixes must be .fasta / .fas, .gb\n";
		print "                        .nexus / .nex, .phylip / .ph, or .seal as\n";
		print "                        appropriate)\n";
		print "         -g = do not convert ~ gap characters (e.g., from BioEdit) to -\n";
		print "         -G = convert flanking gaps to Ns\n";
		print "         -H = convert sequence data to haplotype data\n";
		print "         -i<f|g|n|p|s> = format of sequence file (fasta (f), GenBank (g),\n";
		print "                         nexus (n), phylip (p), or Se-Al (s)); default =\n";
		print "                         autodetect\n";
		print "         -j = produce jackknifed data sets each with a single sequence\n";
		print "              deleted\n";
		print "         -l<number> = interleave nexus- and/or phylip-formatted output\n";
		print "                      (default = sequential) with specified sequence length\n";
		print "                      (between 10 and 100 inclusive; default = 80);\n";
		print "                      NOTE: also sets default interleave length for fasta\n";
		print "                            output\n";
		print "         -n = convert ambiguous nucleotides to Ns\n";
		print "         -o<f|n|pc|pe|s> = output results in any or all of fasta (f), nexus\n";
		print "                           (n), classic or extended phylip (pc or pe),\n";
		print "                           and/or Se-Al (s) formats (NOTE: flag can be\n";
		print "                           invoked multiple times)\n";
		print "         -r<a|i> = order sequences in final output alphabetically by name\n";
		print "                   (a; default) or in input order from file (i)\n";
		print "         -s = split data set into individual partitions according to nexus\n";
		print "              charset statements\n";
		print "         -t = translate sequences into amino acids\n";
		print "         -u = interactive user-input mode\n";
		print "         -h = print this message and quit\n";
		print "         -v = verbose output\n";
		exit(0);
		}
	else	# Process switches
		{
		for (my $i = 0; $i <= $#ARGV; $i++)
			{
			if ($ARGV[$i] eq "-a")
				{
				$accPrint = 1;
				}
			elsif ($ARGV[$i] =~ /-c(\d+)/)
				{
				$globalGenCode = $1;
				if ($globalGenCode != 0 and not defined $genCodes{$globalGenCode})
					{
					print "Don't understand argument: $ARGV[$i]\n";
					print "Usage: seqConverter.pl -d<filename> -o<f|n|pc|pe|s> [-a] [-c<number>] [-g]\n";
					print "       [-H] [-i<f|g|n|p|s>] [-j] [-l<number>] [-n] [-r<a|i>] [-s] [-t] [-u]\n";
					print "       [-v] [-h]\n";
					print "Version: $version\n";
					exit(1);
					}
				$globalOverride = 1 if ($ARGV[$i] =~ /\*$/);
				}
			elsif ($ARGV[$i] =~ /^-d(.*)/)
				{
				$readFile = $1;
				unless ($readFile eq "*")
					{
					die "ERROR: Data file $readFile does not exist.\n" unless (-e $readFile);
					}
				push @seqFiles, $readFile;
				}
			elsif ($ARGV[$i] eq "-g")
				{
				$gapChange = 0;
				}
			elsif ($ARGV[$i] eq "-G")
				{
				$flankGap = 1;
				}
			elsif ($ARGV[$i] eq "-H")
				{
				$haploTyping = 1;
				}
			elsif ($ARGV[$i] eq "-if")
				{
				$inputType = "fasta";
				}
			elsif ($ARGV[$i] eq "-ig")
				{
				$inputType = "GenBank";
				}
			elsif ($ARGV[$i] eq "-in")
				{
				$inputType = "nexus";
				}
			elsif ($ARGV[$i] eq "-ip")
				{
				$inputType = "phylip";
				}
			elsif ($ARGV[$i] eq "-is")
				{
				$inputType = "Se-Al";
				}
			elsif ($ARGV[$i] eq "-j")
				{
				$jackknife = 1;
				}
			elsif ($ARGV[$i] eq "-k")
				{
				$nameClean = 0;
				}
			elsif ($ARGV[$i] =~ /-l(\d+)?/)
				{
				$interleaveLength = $1 if (defined $1);
				$outFormat = "interleave";
				}
			elsif ($ARGV[$i] eq "-n")
				{
				$ambigChange = 1;
				}
			elsif ($ARGV[$i] eq "-of")
				{
				$fastaPrint = 1;
				}
			elsif ($ARGV[$i] eq "-on")
				{
				$nexusPrint = 1;
				}
			elsif ($ARGV[$i] eq "-opc")
				{
				$phylipTradPrint = 1;
				}
			elsif ($ARGV[$i] eq "-ope")
				{
				$phylipExtPrint = 1;
				}
			elsif ($ARGV[$i] eq "-os")
				{
				$sealPrint = 1;
				}
			elsif ($ARGV[$i] eq "-ra")
				{
				$seqOrder = "alphabetical";
				}
			elsif ($ARGV[$i] eq "-ri")
				{
				$seqOrder = "input";
				}
			elsif ($ARGV[$i] eq "-s")
				{
				$seqSplit = 1;
				}
			elsif ($ARGV[$i] eq "-t")
				{
				$translateSeqs = 1;
				$globalGenCode = 1 unless ($globalGenCode);
				}
			elsif ($ARGV[$i] eq "-v")
				{
				$verbose = 1;
				}
			elsif ($ARGV[$i] eq "-x")
				{
				$debug = 1;
				$verbose = 1;
				}
			else
				{
				print "Don't understand argument: $ARGV[$i]\n";
			print "Usage: seqConverter.pl -d<filename> -o<f|n|pc|pe|s> [-a] [-c<number>] [-g]\n";
			print "       [-G] [-H] [-i<f|g|n|p|s>] [-j] [-l<number>] [-n] [-r<a|i>] [-s] [-t]\n";
			print "       [-u] [-v] [-h]\n";
				print "Version: $version\n";
				exit(1); 
				}
			}
		}
	$accPrint = 0 if ($accPrint eq "n");
	geneticCoder() if ($globalGenCode);
	$fastaLength = $interleaveLength;

# Check for I/O errors
	die "ERROR: Must supply name of file containing sequence data.\n" if (not @seqFiles);
	die "ERROR: Must supply at least one output format.\n" unless ($fastaPrint or $nexusPrint or $phylipTradPrint or $phylipExtPrint or $sealPrint);
	die "ERROR: Sequence length for interleaved format must be between 10 and 100 inclusive.\n" if (($fastaPrint or (($nexusPrint or $phylipTradPrint or $phylipExtPrint) and $outFormat eq "interleave")) and ($interleaveLength < 10 or $interleaveLength > 100));

# Read in sequence data
	if ($seqFiles[0] eq "*")	# Batch convert all appropriate files in working directory
		{
		if ($inputType)
			{
			undef @seqFiles;	

			system("ls -d * > convertList.txt");
				setLineBreak("convertList.txt");
				open (LIST, "<convertList.txt") or die "Cannot open file containing names of all sequence files, convertList.txt\n";
					while (<LIST>)
						{
						chomp;
						next unless ($_);
						push @seqFiles, $_ if ($inputType eq "Se-Al" and $_ =~ /\.seal$/);
						push @seqFiles, $_ if ($inputType eq "phylip" and ($_ =~ /\.phylip$/ or $_ =~ /\.ph$/));
						push @seqFiles, $_ if ($inputType eq "nexus" and ($_ =~ /\.nexus$/ or $_ =~ /\.nex$/));
						push @seqFiles, $_ if ($inputType eq "GenBank" and ($_ =~ /\.gb$/));
						push @seqFiles, $_ if ($inputType eq "fasta" and ($_ =~ /\.fasta$/ or $_ =~ /\.fas$/));
						}
				close LIST;
			unlink ("convertList.txt") unless ($debug);
	
			die "ERROR: No files of type $inputType found for batch conversion.\n" if (not @seqFiles);
			}
		else
			{
			die "ERROR: Must specify input file type for batch conversion\n";
			}
		}

	foreach my $seqFile (@seqFiles)
		{
		print "\nConverting file $seqFile ...\n";

		# Set output file names
			$dataSource = $seqFile;
				$dataSource =~ s/\.\w+$//;
		
			if ($fastaPrint)
				{
				$fastaOut = $dataSource . ".fasta";
					$fastaOut =~ s/\.fasta$/_new.fasta/ if ($fastaOut eq $seqFile);
					$fastaOut =~ s/\.fasta$/_haplo.fasta/ if ($haploTyping);
				}
			if ($nexusPrint)
				{
				$nexusOut = $dataSource . ".nex";
					$nexusOut =~ s/\.nex$/_new.nex/ if ($nexusOut eq $seqFile);
					$nexusOut =~ s/\.nex$/_haplo.nex/ if ($haploTyping);
				}
			if ($phylipTradPrint or $phylipExtPrint)
				{
				$phylipOut = $dataSource . ".phylip";
					$phylipOut =~ s/\.phylip$/_new.phylip/ if ($phylipOut eq $seqFile);
					$phylipOut =~ s/\.phylip$/_haplo.phylip/ if ($haploTyping);
				}
			if ($sealPrint)
				{
				$sealOut = $dataSource . ".seal";
					$sealOut =~ s/\.seal$/_new.seal/ if ($sealOut eq $seqFile);
					$sealOut =~ s/\.seal$/_haplo.seal/ if ($haploTyping);
				}
				
			$haploFile = $dataSource . "_haplotypeSeqs.txt";
			
			$phylipExtPrint = 0 if ($phylipTradPrint);
		
		# Read in sequence data
			# Clear variables
				undef @accNum;
				undef %nameLabel;
				undef %sequence;
				undef %geneticCode;
				undef %accPresent;
				undef %deletedSeq;
				undef %finalSeq;
				undef $seqCount;
				undef $ntax;
				undef @charsetList;
					undef %charsetStart;
					undef %charsetEnd;
				my (@haploList, %haploID, %haploSeqs);
				$maxLength = 0;

			seqRead($seqFile);
			if (not @accNum)
				{
				print "\tERROR: Could not read in sequences from file $seqFile; skipping to next file\n";
				next;
				}

		# Process for printing
			my $stopCodonCount;
			$seqType = "protein" if ($globalGenCode and $translateSeqs);
			foreach my $seq (@accNum)
				{
				$sequence{$seq} =~ s/\~/-/g if ($gapChange);
				$sequence{$seq} =~ s/R|Y|M|K|S|W|H|B|V|D/N/ig if ($ambigChange and $seqType eq "nucleotide");

				$finalSeq{$seq} = $sequence{$seq};
					$geneticCode{$seq} = $globalGenCode if ($globalOverride);
					if ($globalGenCode and $translateSeqs)
						{
						$finalSeq{$seq} = translate($sequence{$seq}, $geneticCode{$seq});
						$stopCodonCount++ if (($finalSeq{$seq} =~ tr/\*//) > 2);	# Check for stop codons
						}
				$maxLength = length($finalSeq{$seq}) if (length($finalSeq{$seq}) > $maxLength);
				$nameLabel{$seq} =~ s/\W+/_/g if ($nameClean);	# Clean sequence labels of non-alphanumerics
				}
			printf "\n\tWARNING: $stopCodonCount of %s sequences had more than two stop codons; check for\n\t\t1) proper reading frame,\n\t\t2) proper genetic code, or\n\t\t3) that sequences are coding DNA\n", scalar(@accNum) if ($globalGenCode and $stopCodonCount);
	
		# Add gaps to end of any sequence less than maximum length
			foreach my $seq (@accNum)
				{
				$finalSeq{$seq} .= "-" x ($maxLength - length($sequence{$seq}));
					$sequence{$seq} = $finalSeq{$seq};
				}

		# Convert flanking gaps to Ns if desired
			if ($flankGap)
				{
				foreach my $seq (@accNum)
					{
					if ($finalSeq{$seq} =~ /^(\-+)/)
						{
						my $startGap = $1;
							my $startN = "N" x length($startGap);
						$finalSeq{$seq} =~s/^$startGap/$startN/;
						}
					if ($finalSeq{$seq} =~ /(\-+)$/)
						{
						my $endGap = $1;
							my $endN = "N" x length($endGap);
						$finalSeq{$seq} =~s/$endGap$/$endN/;
						}
					}
				}
				
		# Determine haplotypes as needed
			if ($haploTyping)
				{
				foreach my $entry (@accNum)
					{
					$haploID{$entry} = 0;
					
					foreach my $haploType (@haploList)
						{
						if ($finalSeq{$entry} eq $finalSeq{$haploType})	# Matches existing haplotype; add to list
							{
							$haploID{$entry} = $haploType;
							push @{$haploSeqs{$haploType}}, $entry;
							last;
							}
						}
						
					if (not $haploID{$entry})	# No match to existing haplotype; define new haplotype
						{
						$haploID{$entry} = "haplo" . (scalar(@haploList) + 1);
						push @haploList, $haploID{$entry};
						$nameLabel{$haploID{$entry}} = $haploID{$entry};
						$finalSeq{$haploID{$entry}} = $finalSeq{$entry};
						push @{$haploSeqs{$haploID{$entry}}}, $entry;
						}
					}
					
				undef @accNum;
					@accNum = @haploList;
					
				open (HAPLO, ">$haploFile") or die "Cannot print to file $haploFile\n";
					foreach my $haploType (@haploList)
						{
						print HAPLO "$haploType:";
						print HAPLO "\t$nameLabel{$_}" foreach @{$haploSeqs{$haploType}};
						print HAPLO "\n";
						}
				close HAPLO;
				}
	
		# Print results!
			print "\nPrinting results ...\n";
			@accNum = sort { $nameLabel{$a} cmp $nameLabel{$b} } keys %nameLabel if ($seqOrder eq "alphabetical" and not $haploTyping);

			# Print full data set
				$ntax = scalar(@accNum);
					seqPrint($seqFile);

			# Print jackknifed data sets
				if ($jackknife)
					{
					my $delseqCount = 0;
					foreach my $seq (@accNum)
						{
						$deletedSeq{$seq} = 1;
							$delseqCount++;

						# Change output file names
							$fastaOut = $dataSource . "_jack$delseqCount.fasta";
							$nexusOut = $dataSource . "_jack$delseqCount.nex";
							$phylipOut = $dataSource . "_jack$delseqCount.phylip";
							$sealOut = $dataSource . "_jack$delseqCount.seal";

						$ntax = scalar(@accNum) - 1;
						seqPrint($seqFile);

						$deletedSeq{$seq} = 0;	# Reset deleted sequence
						}
					}

			# Print data set partitions
				if (@charsetList and $seqSplit)
					{
					my $delCount = 0;
					foreach my $partition (@charsetList)
						{
						$delCount = 0;
						# Change output file names
							$fastaOut = $dataSource . "_$partition.fasta";
							$nexusOut = $dataSource . "_$partition.nex";
							$phylipOut = $dataSource . "_$partition.phylip";
							$sealOut = $dataSource . "_$partition.seal";
						
						# Restrict sequence to partition limits
							foreach my $seq (@accNum)
								{
								$deletedSeq{$seq} = 0;	# Reset all deleted sequences
								$finalSeq{$seq} = substr($sequence{$seq}, $charsetStart{$partition} - 1, $charsetEnd{$partition} - $charsetStart{$partition} + 1);
								# Check that sequence remains informative
									unless ($finalSeq{$seq} =~ /a/i or $finalSeq{$seq} =~ /c/i or $finalSeq{$seq} =~ /g/i or $finalSeq{$seq} =~ /t/i)
										{
										$delCount++;
										$deletedSeq{$seq} = 1;
										}
								}
						$ntax = scalar(@accNum) - $delCount;
						$maxLength = $charsetEnd{$partition} - $charsetStart{$partition} + 1;
						seqPrint($seqFile);
						}
					}
		}

exit(0);

### Subroutines used in the program

sub setLineBreak	# Check line breaks of input files and set input record separator accordingly
	{
	my $inFile = shift;
	$/ ="\n";
	open (IN, "<$inFile") or die "Cannot open $inFile to check form of line breaks.\n";
		while (<IN>)
			{
			if ($_ =~ /\r\n/)
				{
				print "\tDOS line breaks detected ...\n" if ($verbose);
				$/ ="\r\n";
				last;
				}
			elsif ($_ =~ /\r/)
				{
				print "\tMac line breaks detected ...\n" if ($verbose);
				$/ ="\r";
				last;
				}
			else
				{
				print "\tUnix line breaks detected ...\n" if ($verbose);
				$/ ="\n";
				last;
				}
			}
	close IN;
	}

sub seqRead
	{
	my $seqFile = shift;
	undef %sequence;

	print "\nReading in sequence data from file $seqFile (type is $inputType) ...\n" if ($inputType);
	setLineBreak($seqFile);
	open (SEQ, "<$seqFile") or die "Cannot open file containing sequences, $seqFile\n";
		my ($header, $tempAcc, $tempName, $tempSeq);
		my $fastaAcc;
		my (%nexusSpecies, %nexusAcc, $nexusRead, $commentFlag);
		my ($phylipLineCount, $phylipTaxa, $phylipChars, %phylipSeq);
		my $sealCode;
		my ($sealDelFlag, $owner) = (0, 0);
		my ($gbAcc, $gbRead);
		my $macBlock = 0;
		my %accCount;

		while (<SEQ>)
			{
			chomp;
			my $lineRead = $_;
			next unless ($lineRead);
			
			# Autodetect sequence format
				if (not $inputType)
					{
					$inputType = "fasta" if ($lineRead =~ /^>/);
					$inputType = "nexus" if ($lineRead =~ /\#nexus/i);
					$inputType = "phylip" if ($lineRead =~ /^\s*\d+\s+\d+/);
					$inputType = "Se-Al" if ($lineRead =~ /^\s*Database=\{/i);
					$inputType = "GenBank" if ($lineRead =~ /^\s*LOCUS/);
					print "\nReading in sequence data from file $seqFile (type determined to be $inputType) ...\n" if ($inputType);
					}
			
			if ($inputType eq "nexus")
				{
				# Check if charset statement present
					if ($lineRead =~ /charset/i)
						{
						$lineRead =~ s/\s+//g;
						$lineRead =~ s/;$//;

						my ($charsetName, $charsetBounds) = split('=', $lineRead);
							$charsetName =~ s/charset//i;
								push @charsetList, $charsetName;
							if ($charsetBounds =~ /(\d+)-*(\d*)/)
								{
								$charsetStart{$charsetName} = $1;
								if ($2)
									{
									$charsetEnd{$charsetName} = $2;
									}
								else
									{
									$charsetEnd{$charsetName} = $1;
									}
								}
						}
						
				# Otherwise block out MacClade / PAUP blocks
					$macBlock = 1 if ($lineRead =~ /begin macclade;/i or $lineRead =~ /begin paup;/i);
					$macBlock = 0 if ($macBlock and $lineRead =~ /end;/i);
					next if ($macBlock);
						
				# Otherwise only read in data lines or charset statements
					if ($lineRead =~ /^\s*matrix/i)
						{
						$nexusRead = 1;
						next;
						}
					$nexusRead = 0 if ($lineRead =~ /;\s*$/);
					next unless ($nexusRead);
					
					# Remove MacClade sequence lengths
						$lineRead =~ s/\[\d+\]$//;

					$commentFlag = 1 if ($lineRead =~ /\[/);
					if ($lineRead =~ /\]/)
						{
						$commentFlag = 0;
						next;
						}
					next if ($commentFlag);

					next unless ($lineRead =~ /a/i or $lineRead =~ /c/i or $lineRead =~ /g/i or $lineRead =~ /t/i or $lineRead =~ /n/i or $lineRead =~ /\?/ or $lineRead =~ /\-/);

				# Clean up input line
					$lineRead =~ s/^\s+//;
					$lineRead =~ s/\'//g;

				my @nexusLine = split(/\s+/, $lineRead);
					my $species = shift(@nexusLine);
						$species =~ s/\s+/_/g;
						$species =~ s/\_+/_/g;
					my $seq = join('', @nexusLine);
						$seq =~ s/\s+//g;
						$seqType = "protein" if ($seq =~ /E/i or $seq =~ /Q/i or $seq =~ /I/i or $seq =~ /L/i or $seq =~ /F/i or $seq =~ /P/i);
				if (not defined $nexusSpecies{$species})
					{
					$nexusSpecies{$species} = 1;
					$seqCount++;
					$nexusAcc{$species} = "tAlign_".$seqCount;
					push @accNum, $nexusAcc{$species};
						$nameLabel{$nexusAcc{$species}} = $species;
						$sequence{$nexusAcc{$species}} = uc($seq);
						$geneticCode{$nexusAcc{$species}} = $globalGenCode;
					}
				else	# Sequences are in interleaved format; append sequence
					{
					$sequence{$nexusAcc{$species}} .= uc($seq);
					}
				}

			if ($inputType eq "fasta")
				{
				if ($lineRead =~/^\s*>/)
					{
					my $species;
					$seqCount++;
					(my $tempSpecies = $lineRead) =~ s/^\s*>//;
					
						if ($tempSpecies =~ /^Mit\.\s+/)	# Entry comes from European RNA project
							{
							$tempSpecies =~ s/^Mit\.\s+//i;	# To fix entries from European RNA project
							my @speciesInfo = split(/\s+/, $tempSpecies);
								$species = join('_', $speciesInfo[0], $speciesInfo[1]);
							if (defined $speciesInfo[2])
								{
								$fastaAcc = $speciesInfo[2];
								$accCount{$fastaAcc}++;
									if ($accCount{$fastaAcc} > 1)
										{
										print "\nWARNING: Accession number $fastaAcc used more than once";
											print "; skipping this and subsequent entries" if ($skipDuplicates);
										print "\n";
										}
								}
							else
								{
								$fastaAcc = "tAlign_".$seqCount;
								$accCount{$fastaAcc}++;
								}
							}
						else
							{
							my @speciesLine = split(/\s+/, $tempSpecies);
							if ($speciesLine[$#speciesLine] =~ /^\(?\w+\d+\)?$/ and scalar(@speciesLine) > 1)	# Check whether last entry is an accession number
								{
								$fastaAcc = pop (@speciesLine);
								$fastaAcc =~ s/^\(//g;
								$fastaAcc =~ s/\)$//g;
								$accCount{$fastaAcc}++;
									if ($accCount{$fastaAcc} > 1)
										{
										print "\nWARNING: Accession number $fastaAcc used more than once";
											if ($skipDuplicates)
												{
												print "; skipping this and subsequent entries\n";
												}
											else
												{
												print "; assigning temporary accession number\n";
												$fastaAcc = "tAlign_" . $seqCount;
												}
										}
								}
							else
								{
								$fastaAcc = "tAlign_".$seqCount;
								$accCount{$fastaAcc}++;
								}
							$species = join('_', @speciesLine);
								$species = "Sequence_".$seqCount if ($species eq "");
							}
					push @accNum, $fastaAcc unless ($accCount{$fastaAcc} > 1 and $skipDuplicates);
						$geneticCode{$fastaAcc} = $globalGenCode;
					($nameLabel{$fastaAcc} = $species) =~ s/\s+/_/;
						$nameLabel{$fastaAcc} =~ s/\_+/_/;
					}
				else
					{
					next if ($accCount{$fastaAcc} > 1 and $skipDuplicates);
					$lineRead =~ s/\s+//g;
						$seqType = "protein" if ($lineRead =~ /E/i or $lineRead =~ /Q/i or $lineRead =~ /I/i or $lineRead =~ /L/i or $lineRead =~ /F/i or $lineRead =~ /P/i);
					$sequence{$fastaAcc} .= uc($lineRead);
					}
				}

			if ($inputType eq "Se-Al")
				{
				my $header;
				$sealDelFlag = 1 if ($lineRead =~/MCoL/);	# Se-Al sometimes places deleted species at end of file; do not read in remainder of file
					next if ($sealDelFlag == 1);
				next unless ($lineRead =~/NumSites/i or $lineRead =~/Owner/i or $lineRead =~/Name/i or $lineRead =~/Accession/i or $lineRead =~/Sequence/i or $lineRead =~/GeneticCode/i or $lineRead =~/Frame/i);
				if ($lineRead =~/Owner\s*\=\s*(\d+)/i)
					{
					$owner = $1;
					}
				if ($lineRead =~/Accession/i and $owner == 2)
					{
					$seqCount++;
					if ($lineRead =~ /null/ or $lineRead =~ /\"\"/)
						{
						$tempAcc = "tAlign_" . $seqCount;
						$accCount{$tempAcc}++;
						}
					else
						{
						($header, $tempAcc) = split (/=/, $lineRead);
							$tempAcc =~ s/\"//g;
							$tempAcc =~ s/;//g;
							$accCount{$tempAcc}++;
								if ($accCount{$tempAcc} > 1)
									{
									print "\nWARNING: Accession number $fastaAcc used more than once";
										if ($skipDuplicates)
											{
											print "; skipping this and subsequent entries\n";
											}
										else
											{
											print "; assigning temporary accession number\n";
											$tempAcc = "tAlign_" . $seqCount;
											}
									print "\n";
									}
						}
					push @accNum, $tempAcc unless ($accCount{$tempAcc} > 1 and $skipDuplicates);
					}
				if ($lineRead =~/Name/i and $owner == 2)
					{
					($header, $tempName) = split (/=/, $lineRead);
						$tempName =~ s/\"//g;
						$tempName =~ s/\s*;//g;
						$tempName =~ s/\s+/_/g;
						$tempName =~ s/\_+/_/g;
					}
				if ($lineRead =~/GeneticCode=(\d)/i and $owner == 2)
					{
					$sealCode = $1;
					$geneticCode{$tempAcc} = $seal2gb{$sealCode};
					}
				if ($lineRead =~/Sequence/i and $owner == 2)
					{
					next if ($accCount{$tempAcc} > 1 and $skipDuplicates);
					($header, $tempSeq) = split (/=/, $lineRead);
						$tempSeq =~ s/\"//g;
						$tempSeq =~ s/;//g;
						$tempSeq =~ s/\s+//g;
					$nameLabel{$tempAcc} = $tempName;
					$sequence{$tempAcc} = uc($tempSeq);
						$seqType = "protein" if ($tempSeq =~ /E/i or $tempSeq =~ /Q/i or $tempSeq =~ /I/i or $tempSeq =~ /L/i or $tempSeq =~ /F/i or $tempSeq =~ /P/i);
					}
				if ($lineRead =~/Frame=(\d)/i)	# Correct for reading frame
					{
					my $readingFrame = $1;
					$sequence{$tempAcc} = "--" . $sequence{$tempAcc} if ($readingFrame == 2);
					$sequence{$tempAcc} = "-" . $sequence{$tempAcc} if ($readingFrame == 3);
					}
				}

			if ($inputType eq "phylip")
				{
				if ($lineRead =~ /^\s*(\d+)\s+(\d+)/)
					{
					$phylipTaxa = $1;
					$phylipChars = $2;
					$phylipLineCount = 0;
					}
				else
					{
					$phylipLineCount++;
					
					$lineRead =~ s/\s//g;
					
					$phylipSeq{$phylipLineCount} .= $lineRead;
					
					$phylipLineCount = 0 if ($phylipLineCount == $phylipTaxa);
					}
				}

			if ($inputType eq "GenBank")
				{
				# Get species name and accession number
					# Pure GenBank format
						$gbAcc = $1 if ($lineRead =~ /^\s*ACCESSION\s+(\w+\d+)/);
						if ($lineRead =~ /^\s*ORGANISM\s+/)
							{
							$seqCount++;
								$gbAcc = "tAlign_" . $seqCount if (not defined $gbAcc);
							$lineRead =~ s/^\s+//;
							my @orgLine = split(/\s+/, $lineRead);
								my $header = shift (@orgLine);
								$nameLabel{$gbAcc} = join('_', @orgLine);
								$accCount{$gbAcc}++;
							}
					# BioEdit format
						if ($lineRead =~ /^\s*TITLE/ and not defined ($gbAcc))
							{
							$seqCount++;
								$gbAcc = "tAlign_" . $seqCount;
							my @accLine = split (/\s+/, $lineRead);
								$gbAcc = $1 if ($accLine[2] =~ /^\((\w+\d+)\)/);
							$nameLabel{$gbAcc} = $accLine[1];
							$accCount{$gbAcc}++;
							}

				if ($lineRead =~ /^\s*ORIGIN/)
					{
					if ($accCount{$gbAcc} > 1)
						{
						print "\nWARNING: Accession number $gbAcc used more than once";
						if ($skipDuplicates)
							{
							print "; skipping this and subsequent entries\n";
							next;
							}
						else
							{
							print "; assigning temporary accession number\n";
							$gbAcc = "tAlign_" . $seqCount;
							$gbRead = 1;
							}
						}
					else
						{
						$gbRead = 1;
						$seqCount++;
						}
					next;
					}

				if ($lineRead =~ /^\s*\/\//)	# End of accession; process
					{
					$gbRead = 0;
					next unless ($gbAcc);
					
					push @accNum, $gbAcc unless ($accCount{$gbAcc} > 1);
					$geneticCode{$gbAcc} = $globalGenCode;
					$nameLabel{$gbAcc} =~ s/\s+/_/g;
						$nameLabel{$gbAcc} =~ s/\_+/_/g;
					$sequence{$gbAcc} =~ s/\d//g;
						$sequence{$gbAcc} =~ s/\s+//g;
						$sequence{$gbAcc} =~ s/\~/-/g if ($gapChange);
					undef $gbAcc;
					}
					
				next unless ($gbRead);
				
				if ($gbRead and $lineRead =~ /^\s+\d+/)
					{
					$sequence{$gbAcc} .= uc($lineRead);
					}
				}
			}
	close SEQ;
	
	if ($inputType eq "phylip")	# Postprocess input to derive taxon names and sequence; accounts for both sequential and extended formatting
		{
		for (my $i = 1; $i <= $phylipTaxa; $i++)
			{
			my $phylipAcc = "tAlign_" . $i;
			
			push @accNum, $phylipAcc;
			$geneticCode{$phylipAcc} = $globalGenCode;
			
			# Derive taxon name and sequence
				$sequence{$phylipAcc} = uc(substr($phylipSeq{$i}, 0 - $phylipChars));
					$seqType = "protein" if ($sequence{$phylipAcc} =~ /E/i or $sequence{$phylipAcc} =~ /Q/i or $sequence{$phylipAcc} =~ /I/i or $sequence{$phylipAcc} =~ /L/i or $sequence{$phylipAcc} =~ /F/i or $sequence{$phylipAcc} =~ /P/i);

				$nameLabel{$phylipAcc} = substr($phylipSeq{$i}, 0, length($phylipSeq{$i}) - $phylipChars);
					$nameLabel{$phylipAcc} =~ s/\s+/_/g;
					$nameLabel{$phylipAcc} =~ s/\_+/_/g;
			}
		}
	}
	
sub seqPrint
	{
	my $seqFile = shift;
	
	$interleaveLength = $maxLength if ($outFormat eq "sequential");

	# Print fasta-formatted file
		if ($fastaPrint)
			{
			print "\tWriting to fasta-formatted file $fastaOut ...\n";
			open (FASTA, ">$fastaOut") or die "Cannot open fasta file for aligned DNA sequences, $fastaOut";
				foreach my $entry (@accNum)
					{
#					next if ($deletedSeq{$entry} or not defined $sequence{$entry});
					print FASTA ">$nameLabel{$entry}";
						print FASTA "\t($entry)" unless ($entry =~ /^tAlign/);
						print FASTA "\n";
						
					# Print sequence
						my $fastaSeq = $finalSeq{$entry};
						for (my $breakpoint = 0; $breakpoint <= length($fastaSeq); $breakpoint += $fastaLength)
							{
							print FASTA substr($fastaSeq, $breakpoint, $fastaLength) . "\n";
							}
#					my $fastaSeq = $finalSeq{$entry};
#						my $breakPoint = 80;
#						my $breakCount = 0;
#						until ($breakPoint > length($fastaSeq))
#							{
#							$breakCount++;
#							my $replaceString = "\n" . substr($fastaSeq, $breakPoint, 1);
#							substr($fastaSeq, $breakPoint, 1) = $replaceString;
#							$breakPoint += 81;	# Latter accounts for addition of \n to string
#							}
#					print FASTA "\n$fastaSeq\n";
					}
			close FASTA;
			}

	# Print nexus-formatted file
		if ($nexusPrint)
			{
			print "\tWriting to nexus file $nexusOut";
				print " in interleaved format" if ($outFormat eq "interleave");
				print " ...\n";
			open (NEX, ">$nexusOut") or die "Cannot open nexus file for aligned DNA sequences, $nexusOut";
				print NEX "#NEXUS\n\n";
				print NEX "[File created from $seqFile using $perlScript v$version on ".localtime()."]\n\n";
				print NEX "begin data;\n";
				print NEX "\tdimensions ntax = $ntax nchar = $maxLength;\n";
				print NEX "\tformat datatype = $seqType gap = - missing = ?";
					print NEX " interleave" if ($outFormat eq "interleave");
					print NEX ";\n\n";
				print NEX "\tmatrix\n\n";
					for (my $interleaveCount = 1; $interleaveCount <= ceil($maxLength / $interleaveLength); $interleaveCount++)
						{
						my $seqStart = (($interleaveCount - 1) * $interleaveLength) + 1;
						my $seqEnd = $interleaveCount * $interleaveLength;
							$seqEnd = $maxLength if ($maxLength <= $seqEnd);
							my $printLength = $seqEnd - $seqStart + 1;
						
						foreach my $entry (@accNum)
							{
							next if ($deletedSeq{$entry});
							my $nexusName = $nameLabel{$entry};
								$nexusName = $entry if ($accPrint);
							if ($nexusName =~ /\W/)
								{
								print NEX "'$nexusName'";
								}
							else
								{
								print NEX "$nexusName";
								}

							my $printSeq = substr($finalSeq{$entry}, $seqStart - 1, $printLength);
								if ($outFormat eq "interleave")	# Add gaps every 10 elements in sequence
									{
									for (my $gapAdd = 100; $gapAdd >= 10; $gapAdd-=10)
										{
										next if $gapAdd > $printLength;
										my $bp = substr($printSeq, $gapAdd - 1, 1) . " ";
										substr($printSeq, $gapAdd - 1, 1, $bp);
										}
									}
								print NEX "\t$printSeq\n";
							}
						print NEX "\n" if ($interleaveCount * $interleaveLength <= $maxLength);
						}
#				foreach my $entry (@accNum)
#					{
#					next if ($deletedSeq{$entry});
#					if ($nameLabel{$entry} =~ /\W/)
#						{
#						print NEX "'$nameLabel{$entry}'";
#						}
#					else
#						{
#						print NEX "$nameLabel{$entry}";
#						}
#					print NEX "\t$finalSeq{$entry}\n";
#					}
				print NEX "\t;\nend;\n";
			close NEX;

			if (-e "/Developer/Tools/SetFile")
				{
				system ("/Developer/Tools/SetFile -t 'TEXT' $nexusOut");
				system ("/Developer/Tools/SetFile -c 'PAUP' $nexusOut");
				}
			}

	# Print phylip-formatted file (on demand)
		if ($phylipTradPrint or $phylipExtPrint)
			{
			my $maxTaxLength = 50;
				$maxTaxLength = 10 if ($phylipTradPrint);
				my $blankName = " " x $maxTaxLength;
			my %shortNameCount;	
				
			print "\tWriting to phylip file $phylipOut ...\n";
			open (PHYLIP, ">$phylipOut") or die "Cannot open phylip file for aligned DNA sequences, $phylipOut";
				print PHYLIP "    $ntax    $maxLength\n";
				for (my $interleaveCount = 1; $interleaveCount <= ceil($maxLength / $interleaveLength); $interleaveCount++)
					{
					my $seqStart = (($interleaveCount - 1) * $interleaveLength) + 1;
					my $seqEnd = $interleaveCount * $interleaveLength;
						$seqEnd = $maxLength if ($maxLength <= $seqEnd);
						my $printLength = $seqEnd - $seqStart + 1;

#					foreach my $entry (@accNum)
#						{
#						my $trimmedName = substr($nameLabel{$entry}, 0, $maxTaxLength);
#							$shortNameCount{trimmedName} = 0;
#						}

					foreach my $entry (@accNum)
						{
						next if ($deletedSeq{$entry});
						
						my $phylipName = $nameLabel{$entry};
							$phylipName = $entry if ($accPrint);
	
						# Print name label as appropriate; also check label and adjust to proper length if needed
							if ($interleaveCount == 1)
								{
								if (length($phylipName) < $maxTaxLength)
									{
									$shortNameCount{$phylipName}++;
									$phylipName .= " " x ($maxTaxLength - length($phylipName)) if ($phylipTradPrint);	# Pad end of name with spaces as needed
									}
								else
									{
									my $trimmedName = substr($phylipName, 0, $maxTaxLength);
									$shortNameCount{$trimmedName}++;
									if ($shortNameCount{$trimmedName} > 1)	# Check for duplicates among shortened names and make unique by adding numbers
										{
										$phylipName = substr($phylipName, 0, $maxTaxLength - length($shortNameCount{$trimmedName}));
											$phylipName .= $shortNameCount{$trimmedName};
											$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
										}
									else
										{
										$phylipName = $trimmedName;
										}
									}
								print PHYLIP "$phylipName";
								print PHYLIP " " if ($phylipExtPrint);
								}
							else
								{
								print PHYLIP "$blankName";
								}
								
						# Print sequence
							my $printSeq = substr($finalSeq{$entry}, $seqStart - 1, $printLength);
								if ($outFormat eq "interleave")	# Add gaps every 10 elements in sequence
									{
									for (my $gapAdd = 100; $gapAdd >= 10; $gapAdd-=10)
										{
										next if $gapAdd > $printLength;
										my $bp = substr($printSeq, $gapAdd - 1, 1) . " ";
										substr($printSeq, $gapAdd - 1, 1, $bp);
										}
									}
								print PHYLIP "$printSeq\n";
						}
					print PHYLIP "\n" if ($interleaveCount * $interleaveLength <= $maxLength);

#				foreach my $entry (@accNum)
#					{
#					next if ($deletedSeq{$entry});
#					
#					my $phylipName = $nameLabel{$entry};
#
#					# Check name label and adjust to proper length if needed
#						if (length($phylipName) < $maxTaxLength)
#							{
#							$shortNameCount{$phylipName}++;
#							$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
#							}
#						else
#							{
#							my $trimmedName = substr($phylipName, 0 , $maxTaxLength);
#							$shortNameCount{$trimmedName}++;
#							if ($shortNameCount{$trimmedName} > 1)	# Check for duplicates among shortened names and make unique by adding numbers
#								{
#								$phylipName = substr($phylipName, 0, $maxTaxLength - length($shortNameCount{$trimmedName}));
#									$phylipName .= $shortNameCount{$trimmedName};
#									$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
#								}
#							else
#								{
#								$phylipName = $trimmedName;
#								}
#							}
#						
#					print PHYLIP "$phylipName";
#						print PHYLIP " " if ($phylipExtPrint);
#					print PHYLIP "$finalSeq{$entry}\n";
					}
			close PHYLIP;
			}

	# Print Se-Al-formatted file (on demand)
		if ($sealPrint)
			{
			print "\tWriting to Se_Al file $sealOut ...\n";
			open (SEAL, ">$sealOut") or die "Cannot open Se-Al file for aligned DNA sequences, $sealOut\n";
				print SEAL "Database={\n";
				print SEAL "\tID='MLst';\n";
				print SEAL "\tOwner=null;\n";
				print SEAL "\tName=null;\n";
				print SEAL "\tDescription=null;\n";
				print SEAL "\tFlags=0;\n";
				print SEAL "\tCount=2;\n";
				print SEAL "\t{\n\t\t{\n";
				
				print SEAL "\t\t\tID='PAli';\n";
				print SEAL "\t\t\tOwner=1;\n";
				print SEAL "\t\t\tName=\"$seqFile\";\n";
				print SEAL "\t\t\tDescription=null;\n";
				print SEAL "\t\t\tFlags=0;\n";
				print SEAL "\t\t\tNumSites=$maxLength;\n";
				print SEAL "\t\t\tType=";
					if ($seqType eq "nucleotide")
						{
						print SEAL "\"Nucleotide\"";
						}
					else
						{
						print SEAL "\"Amino Acid\"";
						}
					print SEAL ";\n";
				print SEAL "\t\t\tFeatures=null;\n";
				print SEAL "\t\t\tColourMode=";
					if ($seqType eq "nucleotide")
						{
						print SEAL "1";
						}
					else
						{
						print SEAL "2";
						}
					print SEAL ";\n";
				print SEAL "\t\t\tLabelMode=0;\n";
				print SEAL "\t\t\ttriplets=false;\n";
				print SEAL "\t\t\tinverse=true;\n";
				print SEAL "\t\t\tCount=$ntax;\n";
				print SEAL "\t\t\t{\n";
				
				my $i = 0;
				foreach my $sequence (@accNum)
					{
					next if ($deletedSeq{$sequence});
					$i++;
					print SEAL "\t\t\t\t{\n";
					print SEAL "\t\t\t\t\tID='PSeq';\n";
					print SEAL "\t\t\t\t\tOwner=2;\n";
					print SEAL "\t\t\t\t\tName=\"$nameLabel{$sequence}\";\n";
					print SEAL "\t\t\t\t\tDescription=null;\n";
					print SEAL "\t\t\t\t\tFlags=0;\n";
					print SEAL "\t\t\t\t\tAccession=";
						if ($sequence =~/^tAlign_/)
							{
							print SEAL "null;\n";
							}
						else
							{
							print SEAL "$sequence;\n";
							}
					if ($seqType eq "nucleotide")
						{
						print SEAL "\t\t\t\t\tType=\"DNA\";\n";
						}
					else
						{
						print SEAL "\t\t\t\t\tType=\"AA\";\n";
						}
					print SEAL "\t\t\t\t\tLength=".length($finalSeq{$sequence}).";\n";
					print SEAL "\t\t\t\t\tSequence=\"$finalSeq{$sequence}\";\n";
					if (defined $geneticCode{$sequence} and $geneticCode{$sequence} != 0)
						{
						print SEAL "\t\t\t\t\tGeneticCode=$gb2seal{$geneticCode{$sequence}};\n";
						}
					else
						{
						print SEAL "\t\t\t\t\tGeneticCode=-1;\n";	# Default for Se-Al is non-coding
						}
					print SEAL "\t\t\t\t\tCodeTable=null;\n";
					print SEAL "\t\t\t\t\tFrame=1;\n";
					print SEAL "\t\t\t\t\tFeatures=null;\n";
					print SEAL "\t\t\t\t\tParent=null;\n";
					print SEAL "\t\t\t\t\tComplemented=false;\n";
					print SEAL "\t\t\t\t\tReversed=false;\n";
					print SEAL "\t\t\t\t}";
					print SEAL "," unless ($i == $ntax);
					print SEAL "\n";
					}
				
				print SEAL "\t\t\t};\n";
				print SEAL "\t\t},\n";
				print SEAL "\t\t{\n";
				print SEAL "\t\t\tID='MCoL';\n";
				print SEAL "\t\t\tOwner=1;\n";
				print SEAL "\t\t\tName=\"Genetic Codes\";\n";
				print SEAL "\t\t\tDescription=\"Custom Genetic Codes\";\n";
				print SEAL "\t\t\tFlags=0;\n";
				print SEAL "\t\t\tCount=0;\n";
				print SEAL "\t\t}\n";
				print SEAL "\t};\n";
				print SEAL "};\n";
			close SEAL;

			if (-e "/Developer/Tools/SetFile")
				{
				system ("/Developer/Tools/SetFile -t 'TEXT' $sealOut");
				system ("/Developer/Tools/SetFile -c 'SEAL' $sealOut");
				}
			}
	}

sub geneticCoder	# Create translation tables for all genetic codes
	{
	my %geneticCode = ('1' => 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '2' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
					   '3' => 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '4' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '5' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
					   '6' => 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '9' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
					   '10' => 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '11' => 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '12' => 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '13' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
					   '14' => 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
					   '15' => 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '16' => 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '21' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
					   '22' => 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '23' => 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG');
	
	foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
		{
		# Establish basic translation table for each genetic code
#			print "\nEstablishing \"$transTable{$code}\" genetic code ...\n" if ($debug);
			my $position = 0;
			foreach my $base1 (qw (T C A G))
				{
				foreach my $base2 (qw (T C A G))
					{
					foreach my $base3 (qw (T C A G))
						{
						my $codon = $base1.$base2.$base3;
						$DNAtoAA{$code}{$codon} = substr($geneticCode{$code}, $position, 1);
#							print "\t$codon = $DNAtoAA{$code}{$codon}\n" if ($debug);
						$position++;
						}
					}
				}
	
		# Extend translation table to account for ambiguity codes (note: does not account for gaps)
#			print "\nExtending translation table to account for ambiguity codes ...\n" if ($debug);
			foreach my $firstPos (@ambigList)
				{
				foreach my $secondPos (@ambigList)
					{
					foreach my $thirdPos (@ambigList)
						{
						my $codon = $firstPos.$secondPos.$thirdPos;
						next if (defined $DNAtoAA{$code}{$codon});
						my $refAA = "";
						foreach my $firstNT (@ {$constitNTlist{$firstPos} })
							{
							last if (defined $DNAtoAA{$code}{$codon});
							foreach my $secondNT (@ {$constitNTlist{$secondPos} })
								{
								last if (defined $DNAtoAA{$code}{$codon});
								foreach my $thirdNT (@ {$constitNTlist{$thirdPos} })
									{
									my $testCodon = $firstNT.$secondNT.$thirdNT;
									if (not $refAA)
										{
										$refAA = $DNAtoAA{$code}{$testCodon};
										}
									else
										{
										if ($DNAtoAA{$code}{$testCodon} ne $refAA)
											{
											$DNAtoAA{$code}{$codon} = "?";
											last;
											}
										}
									}
								}
							}
						$DNAtoAA{$code}{$codon} = $refAA if (not defined $DNAtoAA{$code}{$codon});
#						print "\t$codon = $DNAtoAA{$code}{$codon}\n" if ($debug);
						}
					}
				}
		}
	return;
	}
	
sub translate	# Translate a DNA sequence to an AA sequence (note: does not account for gaps)
	{
	my $DNAseq = shift;
	my $userCode = shift;
	
	my $protSeq;
	for (my $codonStart = 0; $codonStart < length($DNAseq); $codonStart += 3)
		{
		if (length($DNAseq) - $codonStart >= 3)	# Codon is complete; translate
			{
			my $codon = substr($DNAseq, $codonStart, 3);
			if ($codon =~ /-/ or $codon =~ /\./)
				{
				$protSeq .= "?";
				}
			else
				{
				$protSeq .= $DNAtoAA{$userCode}{$codon};
				}
			}
		else	# Incomplete codon; automatically translates as ?
			{
			$protSeq .= "?";
			}
		}

	return $protSeq;
	}

# Version history
#
#	v1.2 (August 12, 2010)
#		- added switch to actually allow fasta output to be specified (or not)
#		- sets TYPE and CREATOR codes for nexus files on Mac systems when
#		  SetFile is present
#		- can now parse GenBank formatted output (both pure GenBank and BioEdit
#		  versions)
#		- error checking: warns if same accession number used more than once
#		  and skips subsequent entires 
#		- added ability to:
#			- detect sequence type of input (nucleotide vs protein) and set as
#			  appropriate in output files
#			- convert nucleotide input to proteins
#			- convert sequence input to haplotypes
#			- interleave nexus- and phylip-formatted output (request by Michael
#			  Craige) and to use inputted value to specify line lengths in fasta
#			  output
#			- output individual data partitions specified according to nexus-
#			  formatted charset statements
#			- output jackknifed data sets, each missing a single taxon
#			- clean sequence labels of non-alphanumeric characters (on by
#			  default)
#			- convert all ~ gap characters (e.g., from BioEdit) to -
#			- convert all ambiguous nucleotides to Ns
#			- change flnaking gaps to Ns (indirect request by Simon Creer)
#		- fixed classic phylip output such that it now conforms 100% to the
#		  phylip guidelines
#		- improved parsing of MacClade generating files, including blocking out
#		  MacClade and PAUP blocks when reading in nexus-formatted files
#		- fixed translation between GenBank and Se-Al genetic codes
#		- changed all instances of #nexus to #NEXUS in output files for
#		  compatibility with TNT (thanks to Douglas Chester for spotting this)
#		- minor bug fixes
#
#	v1.1 (March 2, 2006)
#		- added ability to batch convert all specified file types in working
#		  directory (use -d*)
#		- updated to seqRead module 1.1.1 (includes autodetection of sequence
#		  format)
#		- checks that necessary input file(s) exists before proceeding
#		- added GNU GPL statement
#		- sets TYPE and CREATOR codes for Se-Al files on Mac systems when
#		  SetFile is present
#		- minor bug fixes
#
#	v1.0 (May 30, 2005)
#		- initial release
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.gnu.org/copyleft/gpl.html or by writing to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Fifth Floor, Boston, MA, 02110-1301, USA.

