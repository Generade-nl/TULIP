#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Storable;
use Scalar::Util qw(looks_like_number);

# The Uncorrected Long-read Integration Process, bundling long reads and basic scaffolds

# Bug: now writes NNNs if there is no direct evidence between previosuly terminal seeds (i.e. connected via subterminal seeds)

my $version = "0.4 late 2016 (European eel)";
my @log = ("TULIP bundling $version");
push @log, "Command line ",join(" ",($0,@ARGV));
push @log, "Directory ",cwd();


#### INPUT PARAMETERS
my @mummercoords; my @blasr; my @daligner; my @sam;
my @readsfasta; my @seedsfiles;
my $seedlength = 0;
my $fbuffer = 1000; my $fbuffermax = 5000; 
# my $lengtherror = 0.1; # 10% tolerance	<- currently not used
my $alignmentquality = 0.95; # interpreted as min 95% coverage of the seed
my $alignmentidentity = 0; # %identity, (e.g. 0.95 = 95%). Default 0 = no filtering
my $haploid = 1; my $diploid = 0;	# <- future option?
# my $ploidycutoff = 0.25;
# my $coverage = 1;
my $usedalignerseeds = 0;
my $unusedreads = 0;
my $nobundles = 0;

my $maxreadsinmemory = 100000; # writing starts after filling it up
my $readmappingcutoff = 0.8; # 80% of seeds on a read must belong to the same scaffold
my $extend = 0; # if set, add read sequence beyond the seed boundaries

my $help = 0;
my $outputprefix = "out";
my $configfile;
my $infile;


GetOptions ('coords=s' => \@mummercoords,
			'blasr=s' => \@blasr,
			'daligner=s' => \@daligner,
			'sam=s' => \@sam,
			'seeds=s' => \@seedsfiles,			
			'reads=s' => \@readsfasta,
			'seedlength=s' => \$seedlength,
			'config=s' => \$configfile,		
			'out=s' => \$outputprefix,
			'input=s' => \$infile,
			'filter=s' => \$alignmentquality,
			'identity=s' => \$alignmentidentity,
			#'diploid' => \$diploid,
			'haploid' => \$haploid,
			'memory=s' => \$maxreadsinmemory,
			'mapping=s' => \$readmappingcutoff,
			'extend' => \$extend,
			'allreads' => \$unusedreads,
			'nobundles' => \$nobundles,
			'help' => \$help);
			
if ($help) { printHelp(); } 
if (not $infile) { printHelp("Please specify an --input file prefix"); }
if ($alignmentquality >= 1) { printHelp("The --filter should be between 0 and 1"); }
if ($alignmentquality <0 ) { printHelp("The --filter should be between 0 and 1"); }
if ($alignmentidentity <0) { printHelp("The --identity should be 0 (no filtering) or higher"); }
#if ($coverage < 0) { printHelp("The --coverage should be at least 1"); }
#if ($ploidycutoff < 0) { printHelp("The --minorpath should be between 0 and 1"); }
#if ($ploidycutoff >= 1) { printHelp("The --minorpath should be between 0 and 1"); }
if ($haploid && $diploid) { printHelp("--haploid and --diploid are not compatible"); }
if ((not @mummercoords) && (not @blasr) && (not @daligner) && (not @sam) && (not $configfile)) { printHelp("Please specify alignment coordinates"); }
if (scalar @mummercoords + scalar @blasr + scalar @daligner + scalar @sam > 1) { printHelp("Please use the --config option to specify multiple alignment files"); }
if (not @seedsfiles) { printHelp("Please use --seeds to specify seed sequences"); }
if ((not $configfile) && (not @readsfasta)) { printHelp("Please use --reads or --config to specify read sequences"); }
if ($maxreadsinmemory <1) { $maxreadsinmemory = 1; }
if ($readmappingcutoff <0) { printHelp("The --mapping should be between 0 and 1"); }
if ($readmappingcutoff >1) { printHelp("The --mapping should be between 0 and 1"); }

#### DATA STRUCTURES
my $adjacency;	# note: this is now a reference, so called slightly differently in subroutines than in the sparse_layout equivalents
#my %adjacency;		# two keys seed.exit. Also reciprocal, so $adjacency{1o}{2i} and $adjacency{2i}{1o} both exist
					# entry:	[0] direct count
					#			[1] 0 hypothetical, 1 evidence
					#			[2] minimum gap
					#			[3] mean gap (based on direct counts)
					#			[4]	maximum gap
					#			[5] gap standard deviation
					#			[6] scaffold nr
					# add		[7] sequence
my %scaffoldbundle;	# key counter, entry fasta-format string of reads
my %seedsequence; # key seed, entry sequence
my $seeduse;	# ref to: key seed, entry scaffold nr
my $scaffold;	# ref to: key counter, entry (R|read|fr|from|to) or (M|name.io); replace (MR|...|...) by (MR|seq) later
my %allele;		# key path, entry alternative path
my %fastabuffer; my %recentfastabuffer;

#### MAIN
printLog("Run started");

determineseedlength();
readLayout();
graphStats(); 

bundleReads(); 

readSeeds($usedalignerseeds);
writeScaffolds();

printLog("Run ended");
writeLog();
exit;


#### DEBUG
sub writeGraph { # 0 ignore;
	if ($_[0] == 0) { return(); }
	my $output = $outputprefix.".graph";
	open OUTPUT, ">$output";
	printLog("Writing graph to $output");
	foreach my $k1 (sort keys %$adjacency) {
		foreach my $k2 (sort keys %{$$adjacency{$k1}}) {
			$k1 =~ /^(.+)([io])$/; my @kk1 = ($1,$2);
			$k2 =~ /^(.+)([io])$/; my @kk2 = ($1,$2);
			next if ($kk2[0] lt $kk1[0]);
			print OUTPUT "$kk1[0]\t$kk1[1]$kk2[1]\t$kk2[0]\t",join("\t",@{$$adjacency{$k1}{$k2}}),"\n";
			}
		}
	close OUTPUT;
	return();
}


#### IO SUBROUTINES
sub printHelp {
	if (@_) { print "Error:         $_[0]\n"; }
	print "Usage:         ./bundlesequences.perl --seeds seeds.fasta --reads reads.fasta --sam bwamem.sam --input prefix\n";
	print "               ./bundlesequences.perl --seeds seeds.fasta --config config.txt --input prefix\n\n";
	print "Options:\n";
	print "--out          Output prefix, default \'out\'\n";
	print "--input        Input prefix, required files are .graph_tmp, .seeds_tmp, .scaffolds_tmp\n";
	print "--filter       Ignore alignments with less seed coverage, default = 0.95 (95%)\n";
	print "--allreads     Export unused reads to a bundle (scaffold 0), default off\n";
	print "--nobundles    Do not export read bundles, only scaffolds, default off\n";
	print "--identity     If set, ignore alignments with lesser sequence identity (e.g. 0.9 = 90%), default off (0)\n";
	print "--haploid      Treat graph bubbles as errors to be fixed (default)\n";
	print "--diploid      Treat graph bubbles as alleles\n";
	print "--memory       Number of reads to keep in memory before writing to disk, default = 100000\n";	
	print "--mapping      Fraction of seeds on a read that should belong to the same scaffold, default = 0.8 (80%)\n\n";
	print "Author:        Christiaan Henkel, Leiden University & University of Applied Sciences Leiden (c.v.henkel\@biology.leidenuniv.nl)\n";
	print "Version:       $version\n";
	exit;
}

sub fastaOrFastq {
	my $openfile = $_[0];
	open (TEST, "<".$openfile) or die "Could not open $openfile for reading\n"; 
	my @input;
	for (my $i = 0; $i < 3; $i++) {
		$input[$i] = <TEST>;
		chomp($input[$i]);
		}
	seek (TEST,0,0);
	if ($input[0] =~ /^>/) { return "fasta"; }
	if ($input[0] =~ /^@/) {
		if ($input[2] =~ /^\+/) {
			return ("fastq");
			}
		}
	print "$openfile does not appear to be in either FASTA or FASTQ format\n";
	exit;
}

sub writeLog {
	open LOG, ">$outputprefix".".bundle_log";
	print LOG join("\n",@log),"\n";
	return();
	}

sub printLog {
	my $message = join("",@_);
	chomp($message);
	my @time = localtime();	
	print "$message\n";
	$message = join(":",(sprintf("%02d",$time[2]),sprintf("%02d",$time[1]),sprintf("%02d",$time[0])))."\t".$message; 
	push @log, $message;
	return();
}

sub readLayout {
	my $graphin = $infile.".graph_tmp";
	printLog("Retrieving graph from $graphin");
	$adjacency = retrieve($graphin) || die "Could not read $graphin\n";
	#
	my $scin = $infile.".scaffolds_tmp";
	printLog("Retrieving scaffolds from $scin");
	$scaffold = retrieve($scin) || die "Could not read $scin\n";
	#
	my $seedin = $infile.".seeds_tmp";
	printLog("Retrieving scaffold seeds from $seedin");
	$seeduse = retrieve($seedin) || die "Could not read $seedin\n";
	return();
}

sub graphStats { # note: modified to use a reference to %adjacency
	my %incount;
	my %outcount;
	my %connectioncount;
	my %seeds;
	my @length;
	foreach my $k1 (keys %{$adjacency}) {
		$k1 =~ /^(.+)([io])$/;
		$seeds{$1} = 1;
		if ($2 eq "i") { $incount{scalar keys $$adjacency{$k1}}++; $connectioncount{scalar keys $$adjacency{$k1}}++; }
		if ($2 eq "o") { $outcount{scalar keys $$adjacency{$k1}}++; $connectioncount{scalar keys $$adjacency{$k1}}++; }
		foreach my $k2 (keys %{$$adjacency{$k1}}) {
			$length[0] += $$adjacency{$k1}{$k2}[2];
			$length[1] += $$adjacency{$k1}{$k2}[3];
			$length[2] += $$adjacency{$k1}{$k2}[4];
			}
		}
	printLog(scalar keys %seeds, " seeds in the graph");	
	return ();
}

sub writeReadBundle {
	if ($nobundles) { return(); }
	#
	foreach (keys %scaffoldbundle) {
		my $outbundle = $outputprefix."_readbundle_".$_.".fasta";
		open BUNDLE, ">>$outbundle" || die "Error opening $outbundle\n";
		print BUNDLE $scaffoldbundle{$_};
		close BUNDLE;
		}
	return();
}

sub readSeeds {
	my $dalignermode = shift @_;
	printLog("Adding seed sequences from $seedsfiles[0]");
	open FASTA, "<$seedsfiles[0]\n" || die "Error opening $_\n";
	my $totalseedcount = 0;
	my $counter = 0; 
	my @fasta;
	while (<FASTA>) {
		chomp;
		if (/>/) { 
			($totalseedcount,$counter) = storeSeed($fasta[0],$fasta[1],$dalignermode,$totalseedcount,$counter);
			$fasta[0] = $_;
			$fasta[1] = "";
			}
		else { $fasta[1] .= $_; }
		}
	($totalseedcount,$counter) = storeSeed($fasta[0],$fasta[1],$dalignermode,$totalseedcount,$counter);
	printLog("Sequences found for $totalseedcount seeds");
	return();
}

sub storeSeed {
	(my $header, my $ssequence, my $dmode, my $totcount, my $scount) = @_;
	if (not defined $header) { return ($totcount,$scount); }
	my $seedname;
	if ($dmode) {
		$seedname = $scount++;
		}
	else {	
		$header =~ s/^>//;
		$header =~ s/(\d+)$//;
		$seedname = $1;
		}
	if (exists $$seeduse{$seedname}) {
		$seedsequence{$seedname} = $ssequence;
		$totcount++;
		}
	return($totcount,$scount);  
}

sub writeScaffolds {
	my $outscaffold = $outputprefix."_scaffold_";
	my $sccount = 0; my $seqcount = 0;
	if ($nobundles) { printLog("Writing scaffolds to $outputprefix"."_scaffolds.fasta"); }
	else { printLog("Writing scaffold sequences to $outscaffold"."nr.fasta"); }
	if ($nobundles) { 
		my $currentscaffold = $outputprefix."_scaffolds.fasta";
		open OUTPUT,">$currentscaffold" || die "Error opening $currentscaffold\n";
		}
	foreach my $scid (sort {$a <=> $b} keys %$scaffold) {
		my $currentscaffold = $outscaffold.$scid.".fasta";
		if (not $nobundles) {
			open OUTPUT,">$currentscaffold" || die "Error opening $currentscaffold\n";
			}
		my @scafstack = @{$$scaffold{$scid}}; #split /\t/,$$scaffold{$scid}; # 1003	linear	160942501i 160942501o 160961165o 160961165i 160961501o 160961501i  
		my $scafmode = shift @scafstack;
		my $outseq;
		for (my $i = 1;$i < $#scafstack; $i +=2) {
			my @seed1 = nodeOrientation($scafstack[$i]);
			if ($seed1[1] eq "o") {	
				$outseq .= uc($seedsequence{$seed1[0]});
				}
			else {
				$outseq .= uc(revcomp($seedsequence{$seed1[0]}));
				}
			if ($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[1]) { #observed link
				if (looks_like_number($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[7])) {
					# no sequence between the seeds, they overlap so trim the seed
					if ($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[7] != 0) {
						$outseq = substr($outseq,0,$$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[7]);
						}
					}
				else {	
					if (scalar @{$$adjacency{$scafstack[$i]}{$scafstack[$i+1]}} == 7) { # what happened here...? <- bug
						$outseq .= "n" x int($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[3]+0.5) ;
	#					print "$scid $scafmode $scafstack[$i] $scafstack[$i+1]\t",join(",",@{$$adjacency{$scafstack[$i]}{$scafstack[$i+1]}}),"\n";
						}
					else {	
						$outseq .= lc($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[7]);	
						}
					}
				}
			else { #hypothetical link
				if ($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[3] > -1) {
					$outseq .= 'n' x int($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[3]+0.5);
					print "$scafstack[$i] $scafstack[$i+1]\t",join(",",@{$$adjacency{$scafstack[$i]}{$scafstack[$i+1]}}),"\n";
					}
				else {
					$outseq = substr($outseq,0,int($$adjacency{$scafstack[$i]}{$scafstack[$i+1]}[3]));
					}
				}
			}
		my @lastseed = nodeOrientation($scafstack[-1]);
		if ($lastseed[1] eq "o") {	
			$outseq .= uc($seedsequence{$lastseed[0]});
			}
		else {
			$outseq .= uc(revcomp($seedsequence{$lastseed[0]}));
			}			
		$sccount++;
		$seqcount += length($outseq);
		print OUTPUT ">scaffold_".$scafmode."_$scid\n";
		print OUTPUT "$outseq\n";
		if (not $nobundles) { close OUTPUT; }
		}
	printLog("$seqcount nt written to $sccount scaffolds");	
	if ($nobundles) { close OUTPUT; }
	return();	
}


#### MAIN ALGORITHM SUBROUTINES
sub determineseedlength {
	if (not $seedlength) {
		open FASTA, "<$seedsfiles[0]\n" || die "Error opening $_\n";
		my $firstline = (<FASTA>);
		if ($firstline !~ /^>/) { printHelp("$seedsfiles[0] is not a FASTA file\n"); }
		my $secondline;
		while (<FASTA>) {
			chomp;
			last if ($_ =~ /^>/);
			$secondline .= $_;
			}
		$seedlength = length($secondline);
		}
	printLog("Seed length is $seedlength");	
	return();
}

sub bundleReads {
	my $readcount = 1; my $usedcount = 0;  my $readcounttomax = 0;
	my $currentread = "";
	my $threshold = $seedlength*$alignmentquality; 
	my @filesstack;
	if ($configfile) {
		open CONF, "<$configfile";
		while (<CONF>) {
			my @cnf = split /[\s\t]+/,$_; # dalinger alignmentfile fasta
			if (scalar @cnf != 3) { die "Invalid configuration file: $_\n"; }
			if ($cnf[0] ne "daligner" && $cnf[0] ne "coords" && $cnf[0] ne "blasr" && $cnf[0] ne "sam") { die "Invalid configuration file: $_\n"; }
			push @filesstack, "$cnf[0]\t$cnf[1]\t$cnf[2]";
			}
		close CONF;	
		}
	else {	
		foreach (my $i = 0; $i <= $#mummercoords; $i++) { push @filesstack, "coords\t$mummercoords[$i]\t$readsfasta[$i]"; }
		foreach (my $i = 0; $i <= $#blasr; $i++) { push @filesstack, "blasr\t$blasr[$i]\t$readsfasta[$i]"; }
		foreach (my $i = 0; $i <= $#daligner; $i++) { push @filesstack, "daligner\t$daligner[$i]\t$readsfasta[$i]"; }
		foreach (my $i = 0; $i <= $#sam; $i++) { push @filesstack, "sam\t$sam[$i]\t$readsfasta[$i]"; }
		}
	my $filecounter = 0;
	foreach (@filesstack) {
		my %stack;
		my @openfile = split /\t/,$_;
		$filecounter++;
		my $trackreadnr = 0;		
		printLog("Bundling $openfile[0] coordinates from $openfile[1] with reads from $openfile[2]");
		open COORDS, "<$openfile[1]" || die "Error opening $openfile[1]\n";
		open FASTA, "<$openfile[2]" || die "Error opening $openfile[2]\n";
		if ($openfile[0] eq "daligner") { my $devnull = <COORDS>; $devnull = <COORDS>; $usedalignerseeds = 1; }
		%fastabuffer = ();
		%recentfastabuffer = ();
		while (<COORDS>) {
			chomp;
			my @data;
			if ($openfile[0] eq "daligner") { 
				@data = daligner2Coords($_); # not tab-sep, parse as string
				next if (not @data);
				}
			elsif ($openfile[0] eq "blasr") {
				@data = split /[\t\s]/,$_;
				@data = blasr2Coords(@data);
				}
			elsif ($openfile[0] eq "sam") {
				next if ($_ =~ /^@/);
				@data = split /[\t\s]/,$_;
				@data = sam2Coords(@data);
				}
			else {
				@data = split /[\t\s]/,$_;
				$data[6] = $data[6]/100; # pct identitity scale
				}
			next if (scalar @data == 0);	
			$data[7] =~ /(\d+)$/; # remove sometime, included for a legacy alignment file
			$data[7] = $1; # remove sometime, included for a legacy alignment file
			next if ((not exists $$adjacency{$data[7]."i"}) && (not exists $$adjacency{$data[7]."o"}));
			next if ($data[4] < $threshold);
			@data = extrapolateSeed(@data);
			next if ($alignmentidentity && $data[6] < $alignmentidentity);
			my @positions = sort {$a <=> $b} ($data[2],$data[3]);
			my $orientation;
			if ($positions[0] eq $data[2]) { $orientation = "f"; } else { $orientation = "r"; }
			if ($data[8] eq $currentread) {
				if (not exists $stack{$positions[0]}) { ###########!!!!!!!!!!!!!!!!!!!!!! fix this
					push @{$stack{$positions[0]}}, ($data[0],$data[1],$positions[1],$orientation,$data[7],$data[8]);  # seedstart seedend position orientation seedname readname
					}
				next;
				}
			else {
				if (keys %stack) {
					($usedcount,$trackreadnr) = storeRead(\%stack,$usedcount,$filecounter,$openfile[0],$trackreadnr);
					$readcounttomax++;
					if ($readcounttomax >= $maxreadsinmemory) {
						writeReadBundle(); 
						%scaffoldbundle = ();
						$readcounttomax = 0;
						}
					}
				$currentread = $data[8];
				$readcount++;
				%stack = ();
				push @{$stack{$positions[0]}}, ($data[0],$data[1],$positions[1],$orientation,$data[7],$data[8]); # seedstart seedend position orientation seedname readname
				} # wrap up read 
			}
		if (keys %stack) {	
			($usedcount,$trackreadnr) = storeRead(\%stack,$usedcount,$filecounter,$openfile[0],$trackreadnr);
			}
		writeReadBundle();
		%scaffoldbundle = ();
		close COORDS;
		close FASTA;
		} # each mummer
	printLog("$usedcount long reads bundled for scaffold correction");
	return();
}

sub readFastaDaligner {
	my @inarg = @_; # $readname, $readtrack
	while ($inarg[1] < $inarg[0]-1) {
		my $fheader = <FASTA>; 
		my $fseq = <FASTA>; 
		$inarg[1]++;
		}
	my $fheader = <FASTA>; chomp ($fheader); $fheader =~ s/^>//;
	my $fseq = <FASTA>; chomp ($fseq);
	$inarg[1]++;
	return ($fheader,$fseq,$inarg[1]);
}

sub readFasta {
	my @inarg = @_; # $readname
	if (exists $fastabuffer{$inarg[0]}) { 
		return ($inarg[0],$fastabuffer{$inarg[0]},0);
		}
	if (scalar keys %fastabuffer > $fbuffermax) {
		%fastabuffer = %recentfastabuffer;
		%recentfastabuffer = ();
		}
	for (my $i = 1; $i < $fbuffer; $i++) {
		my $headerin = <FASTA>;
		my $seqin = <FASTA>;
		chomp ($headerin);
		chomp ($seqin);
		$headerin =~ s/^>//;
		$headerin =~ s/\s.+$//;
		if ($headerin eq $inarg[0]) {
			return ($headerin,$seqin,0);
			}
		$fastabuffer{$headerin} = $seqin;
		if (scalar keys %fastabuffer > ($fbuffermax-$fbuffer)) {
			$recentfastabuffer{$headerin} = $seqin;
			}
		}
}

sub storeRead {
	(my $stackref, my $usedcount, my $filecounter, my $aligner, my $readtrack) = @_;
	# @{$stack{$positions[0]}}  ($data[0],$data[1],$positions[1],$orientation,$data[7],$data[8]); # seedstart seedend position orientation seedname readname
	# simply determines whether a read and a scaffold belong together based on seed content
	my %scaffoldsinread; # key scaffold, entry count of seeds
	my $favscaf;
	my @sortkeys = sort {$a <=> $b} keys %$stackref;  
	my $readname = $$stackref{$sortkeys[0]}[5];
	my @storeread; # name/seq/nr
	my @orientations1; my @orientations2;
	# add an entire read to %scaffoldbundle (if criteria met)
	foreach (@sortkeys) {
		$scaffoldsinread{$$seeduse{$$stackref{$_}[4]}}++;
		if ($$stackref{$_}[3] eq "f") {
			push @orientations1,"o"; # for later use
			push @orientations2,"i";
			}
		else {
			push @orientations1,"i";
			push @orientations2,"o";
			}
		}
	foreach (keys %scaffoldsinread) {
		if ($scaffoldsinread{$_} >= $readmappingcutoff) {
			$favscaf = $_;
			$usedcount++;
			if ($aligner eq "daligner") {
				@storeread = readFastaDaligner($readname,$readtrack);
$readtrack = $storeread[2];				
				}
			else {
				@storeread = readFasta($readname);
				}
			last;
			}
		}
	if (@storeread) {
		$scaffoldbundle{$favscaf} .= ">".$storeread[0]."\n".$storeread[1]."\n"; 
		}
	# now add partial sequences to $%adjacency
	INODE: for (my $i = 0; $i < $#sortkeys; $i++) {
		my $nod1 = $$stackref{$sortkeys[$i]}[4].$orientations1[$i];
		if (exists $$adjacency{$nod1}) {
			JNODE: for (my $j = $i+1; $j <= $#sortkeys; $j++) {
				my $nod2 = $$stackref{$sortkeys[$j]}[4].$orientations2[$j];
				if (exists $$adjacency{$nod1}{$nod2}) {
					if ((scalar @{$$adjacency{$nod1}{$nod2}}) < 8) {
						my @newadj = parsePairs($sortkeys[$i],@{$$stackref{$sortkeys[$i]}},$sortkeys[$j],@{$$stackref{$sortkeys[$j]}}); # return: seed.or seed.or gap name from to
						addScaffoldSequence(@newadj,$storeread[1]);
						next INODE;
						}
					}
				}
			}
		}
	return ($usedcount,$readtrack);
}

sub addScaffoldSequence {
	my $sequence = pop @_;
	my @aligninfo = @_; # 0 seed.or 1 seed.or 2 gap 3 name 4 from 5 to
	if ($aligninfo[2] >0) {
		my $insertseq = substr($sequence,$aligninfo[4],$aligninfo[5]-$aligninfo[4]);
		$$adjacency{$aligninfo[0]}{$aligninfo[1]}[7] = $insertseq;
		$$adjacency{$aligninfo[1]}{$aligninfo[0]}[7] = revcomp($insertseq);
		}
	else {
		$$adjacency{$aligninfo[0]}{$aligninfo[1]}[7] = $aligninfo[2];
		$$adjacency{$aligninfo[1]}{$aligninfo[0]}[7] = $aligninfo[2];
		}
	return();
}
	
sub blasr2Coords {
	#	[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[seedS]
	#4 250	12632	12377	247	256	91.19	gi|512322365|gb|CM001609.2|-81703	S1_100
	# qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV	
#[A ID] [B ID] [Jaccard score] [# shared min-mers] [0=A fwd, 1=A rc] [A start] [A end] [A length] [0=B fwd, 1=B rc] [B start] [B end] [B length]
# m140930_013011_sidney_c100699772550000001823139903261597_s1_p0/9/0_17359 seed_20024 95.955511 10.000000 0 6854 7120 17359 0 147 412 1000
	my @mummercoords = ($_[9],$_[10],0,0,$_[10]-$_[9]+1,$_[6]-$_[5]+1,$_[3]/100,$_[1],$_[0]);
	if ($_[8]) { # rev
		($mummercoords[2],$mummercoords[3]) = ($_[6],$_[5]);
		}
	else {
		($mummercoords[2],$mummercoords[3]) = ($_[5],$_[6]);
		}
	return(@mummercoords);
}

sub sam2Coords {
	my @flags = split //,unpack("b11",pack("i*",$_[1]));
	my $clip1 = 0; my $clip2 = 0;
	if ($flags[2]) { return (); } # unmapped
	if ($_[5] =~ s/^(\d+)[HS]//) { $clip1 = $1; } # clipping at beginning of CIGAR
	if ($_[5] =~ s/(\d+)[HS]$//) { $clip2 = $1; } # clipping at end of CIGAR
	my @splitcigar = split /[MDISH]/,$_[5];
	my $cigarseed = 0; my $cigarread = 0; my $nomatchcount = 0;
	while ($_[5]) {
		$_[5] =~ s/(\d+)([MID])//;
		if ($2 eq "M") { $cigarseed += $1; $cigarread += $1; }
		if ($2 eq "I") { $cigarread += $1; $nomatchcount += $1; }
		if ($2 eq "D") { $cigarseed += $1; $nomatchcount += $1; }
		}
	if ($flags[4]) { # reverse
		return($_[3],$_[3]+$cigarseed,$clip2+$cigarread,$clip2+1,$cigarseed,$cigarread,1-($nomatchcount/$cigarseed),$_[2],$_[0]);
		}
	else { # [S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	seed read
		return($_[3],$_[3]+$cigarseed,$clip1+1,$clip1+$cigarread,$cigarseed,$cigarread,1-($nomatchcount/$cigarseed),$_[2],$_[0]);
		}
}

sub daligner2Coords {
	if ($_[0] eq "") { return (); }
	#if ($_[0] =~ /^\S/) { return (); } 
	$_[0] =~ s/,//g; # thousands separators
	$_[0] =~ /(\d+)\s+(\d+)\s+([ncNC])[\s\[]+(\d+)[\.\s]+(\d+)[x\s\]\[]+(\d+)[\.\s]+(\d+)[\]\s\<\:]+(\d+)\sdiffs/;  #    1    9 c   [     0.. 1,876] x [ 9,017..10,825] :   <    398 diffs  ( 18 trace pts)
	# 1seed 2 read 3orientation 4seedstart 5seedend 6readstart 7readend 8diff
	if (lc($3) eq "n") { 
		return ($6+1,$7,$4+1,$5,$7-$6,$5-$4,1-($8/($5-$4)),$2,$1);
		}
	else {
		return ($6+1,$7,$5,$4+1,$7-$6,$5-$4,1-($8/($5-$4)),$2,$1);
		}
}

sub extrapolateSeed {
	if ($_[4] == $seedlength) { return (@_); }
	if ($_[0] > 1) {
		if ($_[2] < $_[3]) {
			$_[2] -= ($_[0]-1);
			$_[0] = 1;
			}
		else {
			$_[2] += ($_[0]-1);
			$_[0] = 1;
			}
		}
	if ($_[1] < $seedlength) {
		if ($_[3] > $_[2]) {
			$_[3] += ($seedlength-$_[1]);
			$_[1] = $seedlength;
			}
		else {
			$_[3] -= ($seedlength-$_[1]);
			$_[1] = $seedlength;
			}
		}
	$_[5] += ($seedlength-$_[4]);
	$_[6] = $_[6] * ($_[4]/$seedlength); # updated quality
	$_[4] = $seedlength;
	return(@_);
}

sub parsePairs {
	#my @newadj = parsePairs($firstalign,@{$stack{$firstalign}},$secondalign,@{$stack{$secondalign}}); # return: seed.or seed.or gap readname from to 
	# 1seedstart 2seedend 3position 4orientation 5seedname 6readname
	# 8seedstart 9seedend 10position 11orientation 12seedname 13readname
	my @returnval;
	if ($_[4] eq "f") {
		if ($_[11] eq "f") { #ff
			@returnval = ($_[5]."o",$_[12]."i");
			}
		else { # fr
			@returnval = ($_[5]."o",$_[12]."o");
			}
		}
	elsif ($_[11] eq "f") { # rf
		@returnval = ($_[5]."i",$_[12]."i");	
		}
	else { # rr
		@returnval = ($_[5]."i",$_[12]."o");
		}
	push @returnval, $_[7]-$_[3]-1;	# gap
	push @returnval, ($_[6],$_[3]+1,$_[7]-1);	# readname from to
	return (@returnval);	
}


#### GRAPH MANIPULATION & GENERIC SUBROUTINES

sub revcomp {
	my $seq = shift(@_);
	if ((length $seq) eq 0) {return "";}
	$seq =~ tr/ATCGN/TAGCN/;
	return (scalar reverse($seq));
}

sub opposite {
	$_[0] =~ /^(.*)([io]$)/;
	if ($2 eq "i") { return ($1."o"); } else { return ($1."i"); }
}

sub nodeOrientation {
	$_[0] =~ /^(.+)([io]$)/;
	return ($1,$2);
}

sub node {
	$_[0] =~ /^(.+)[io]$/;
	return ($1);
}

sub	orientation {
	$_[0] =~ /^.+([io]$)/;
	return ($1);
}

sub min {
	@_ = sort {$a <=> $b} @_;
	return ($_[0]);
}

sub max {
	@_ = sort {$b <=> $a} @_;
	return ($_[0]);
}
