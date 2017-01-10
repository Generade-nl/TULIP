#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Storable;

# The Uncorrected Long-read Integration Process, seed layout

# To do (not urgent):
# Implement diploid layout (? perhaps fix this after read bundling with a scaffold)
# Better determination of typical coverage instead of median/upper quartile
# Fix sub trimNodeLong
# Merge subs collapseLinks & collapseScaffoldLinks


my $version = "0.4 late 2016 (European eel)";
my @log = ("TULIP seed layout $version");
push @log, "Command line ",join(" ",($0,@ARGV));
push @log, "Directory ",cwd();


#### INPUT PARAMETERS

my @mummercoords; my @blasr; my @daligner; my @sam;
my $lengtherror = 0.1; # 10% tolerance
my $alignmentquality = 0.95; # interpreted as min 95% coverage of the seed
my $alignmentidentity = 0; # %identity, (e.g. 0.95 = 95%). Default 0 = no filtering
my $haploid = 0; # currently not used
my $diploid = 1; # currently only used in scaffolding: suspect links of uncertain span are allowed if they are unambiguous
my $ploidycutoff = 0.25; # currently not used
my $seedlength = 0;
my $help = 0;
my $coverage = 1;
my $outputprefix = "out";
my @seedsfiles;
my $configfile;
my $graphmode = 0; # 0 write; 1 ignore
my $scaffolds = 0; # 0 write; 1 ignore
my $debuggraph;


#### HIDDEN PARAMETERS

my $coverageratio = 0.1; # if a-b in two orientations but one has only 10% evidence, delete it (could be e.g. a missed PacBio spacer) # also used in branch decisions
my $lengthinterval = 500; # (more between min and max is highly suspect) # re-interpreted: the SD of all links should be less than this number.
my $subterminaldepth = 5; # for scaffolding: connect terminals based on links to $subterminaldepth subterminal nodes into the scaffold
my %pathbranchingandlength; # ( path complexity x search depth) # complexity = nr of fork/join nodes, excl the initial fork but incl the terminal join. E.g. 1x4 = no branching except initial & terminal
@{$pathbranchingandlength{1}{"cautious"}} = (1,7); # simplifymode 1
@{$pathbranchingandlength{1}{"aggressive"}} = (2,15); # simplifymode 1 
@{$pathbranchingandlength{2}{"cautious"}} = (2,7); # try harder
@{$pathbranchingandlength{2}{"aggressive"}} = (5,15); # try harder
my $linearization = 5; # search depth
my @aggressivetrimming = (6,10,20); # if paths b <+6 and c >10 forking from a, delete b; [2] = search depth
my @loopparameters = (3,10); # maxbranching, maxlength
my $minimumscaffoldlength = 20; # in nodes # not used
my $highcoverage = 10; # nodes with links at $highcoverage*$mediancoverage will be eliminated before even attempting untangling
my %evidenceratio; # if a decision is necessary between two paths, choose one if it has $evidenceratio more evidence than the alternative
$evidenceratio{1} = 10; 
$evidenceratio{2} = 5;
my $coveragepercentile = 0.75; # used for calculating the 'median' coverage (which would be at 0.5)

GetOptions ('coords=s' => \@mummercoords,
			'blasr=s' => \@blasr,
			'daligner=s' => \@daligner,
			'sam=s' => \@sam,
			'graph=s' => \$debuggraph,
			'seeds=s' => \@seedsfiles,			
			'config=s' => \$configfile,		
			'out=s' => \$outputprefix,			
			'error=s' => \$lengtherror,
			'filter=s' => \$alignmentquality,
			'identity=s' => \$alignmentidentity,
			'diploid' => \$diploid,
			'haploid' => \$haploid,
			'minorpath=s' => \$ploidycutoff,
			'coverage=s' => \$coverage,
			'seedlength=s' => \$seedlength,
			'help' => \$help,
			'nograph' => \$graphmode,
			'noscaffolds' => \$scaffolds);
			
if ($help) { printHelp(); } 
if ($lengtherror >= 1) { printHelp("The --error should be between 0 and 1"); }
if ($lengtherror <0 ) { printHelp("The --error should be between 0 and 1"); }
if ($alignmentquality >= 1) { printHelp("The --filter should be between 0 and 1"); }
if ($alignmentquality <0 ) { printHelp("The --filter should be between 0 and 1"); }
if ($alignmentidentity <0) { printHelp("The --identity should be 0 (no filtering) or higher"); }
if ($coverage < 0) { printHelp("The --coverage should be at least 1"); }
if ($ploidycutoff < 0) { printHelp("The --minorpath should be between 0 and 1"); }
if ($ploidycutoff >= 1) { printHelp("The --minorpath should be between 0 and 1"); }
if ($haploid && $diploid) { printHelp("--haploid and --diploid are not compatible"); }

if ($haploid == 0 && $diploid == 0) { $haploid = 1; }
$lengtherror = 1- $lengtherror;	
if (not $debuggraph) {
	if ((not @mummercoords) && (not @blasr) && (not @daligner) && (not @sam) && (not $configfile)) { printHelp("Please specify alignment coordinates"); }
	if (scalar @mummercoords + scalar @blasr + scalar @daligner + scalar @sam > 1) { printHelp("Please use the --config option to specify multiple alignment files"); }
	}
if ((not $seedlength) && (not @seedsfiles)) { printHelp("Please specify a seed length (--seeds or --seedlength)"); }


#### DATA STRUCTURES

my %adjacency;		# two keys seed.exit. Also reciprocal, so $adjacency{1o}{2i} and $adjacency{2i}{1o} both exist
					# entry:	[0] direct count
					#			[1] 0 hypothetical, 1 evidence, 2+ targeted for deleting
					#			[2] minimum gap
					#			[3] mean gap (based on direct counts)
					#			[4]	maximum gap
					#			[5] gap standard deviation
					#			[6] scaffold nr
my %linklength;	# accessory to %adjacency (same keys), holds an array of lenghts for calculating the SD
my %blacklist;	# key seed, entry something
my %terminal;	# key seed, entry seed.orientation
my %subterminal;	# key seed.orientation, (real terminal, entry minimum distance to end (<$maxreadlength))
my $maxreadlength = 0; # determined from alignments
my $mediancoverage = 0; # determined from alignments

my %seeduse;	# key seed, entry scaffold nr
my %scaffold;	# key counter, linear|circular nodeio nodeio ...
my %scaffoldlen;	# key counter, entry length
my %scafadjacency; # analogous to %adjacency, but with scaffold nrs instead of seed nodes; 
my %scafend; # key node.or, entry scaffold io
my %salvagedadjacency; # adjacency that has been deleted but might still be nice to have

my %suspectscaffold;	# key seed.or, entry (all nodes belonging to that scaffold) # endpoints of small scaffolds
my %allele;		# key path, entry alternative path
my %nodestack;	# key type, entry @stack of adjacencies. Used for prioritizing simplification routines
my %loopprotect;	# specific for bubble merging


#### MAIN

printLog("Run started");

determineSeedLength();														# The seed length is fixed, specified either by a FASTA file or the command line
initGraph();																# Nothing much

if ($debuggraph) { readGraph($debuggraph) } else { constructGraph(1); }		# Add seeds aligning adjacently on a long read to the graph
calculateLinkDistances();													# Gap standard deviations
determineCoverage(); 														# Based on the links in the graph, determine a 'typical'  coverage (not pretty)

	
graphStats();																# Node count, edge length
removeSuspectLinks("nodes"); 												# Links between two nodes in two orientations; fuzzy size ranges; very low coverage
removeProbableRepeatNodes();												# Nodes with links that have evidence >> the typical coverage

graphStats();
simplifyGraph(1);															# Collapse triangles, trim tiny tips, collapse simple bubbles, delete implausible connections & repeats
graphStats();

extractScaffolds();															# Used here to get an N50/scaffold count, and as input for the scaffold graph
markTerminalNodes();														# Only terminal nodes are connected in the next round
markSubterminalNodes();														# But subterminal connections are allowed as evidence

constructGraph(2);															# Add all possible links, but only between terminal nodes
calculateLinkDistances();
graphStats();
buildScaffoldGraph();														# An additional graph, using the scaffolds as nodes

simplifyScaffoldGraph(); 													# Collapsing (and...) on the scaffold graph; this deletes links in the seed graph
resetScaffolds();															# 
	
determineCoverage();	
simplifyGraph(2);															# As before, but also longer tips & bubbles, plus different repeat handling (tries to salvage more nodes)
graphStats();

extractScaffolds();
resetScaffolds();	
addSalvagedLinks();															# During scaffold graph simplification, some repeat-spanning links may have been deleted
extractScaffolds();
writeScaffolds($scaffolds); 
writeGraph($outputprefix,$graphmode); # 0 write; 1 ignore
	
printLog("Run ended");
writeLog();
exit;


#### IO SUBROUTINES

sub printHelp {
	print "TULIP-seeds:    The Uncorrected Long-read Intergration Process, seed layout stage\n\n";
	if (@_) { print "Error:         $_[0]\n\n"; }
	print "Usage:         ./tulipseed.perl --seeds seeds.fasta --sam bwamem.sam\n";
	print "               ./tulipseed.perl --seedlength 250 --sam bwamem.sam\n";
	print "               ./tulipseed.perl --seedlength 250 --config configfile.txt\n\n";
	print "Options:\n";
	print "--out          Output prefix, default \'out\'\n";
	print "--error        Tolerance of alignment length differences, default = 0.1 (10%)\n";
	print "--filter       Ignore alignments with less seed coverage, default = 0.95 (95%)\n";
	print "--identity     If set, ignore alignments with lesser sequence identity (e.g. 0.9 = 90%), default off (0)\n";
	print "--haploid      Treat graph bubbles as errors to be fixed (default)\n";
	print "--diploid      Treat graph bubbles as alleles\n";
	print "--minorpath    The minimum path evidence to be treated as allele, default = 0.25 (25%)\n";
	print "--coverage     Ignore links with less evidence, default = 1\n";
	print "Authors:       Christiaan Henkel & Michael Liem, Leiden University & University of Applied Sciences Leiden\n";
	print "Contact:       c.v.henkel\@biology.leidenuniv.nl\n";
	print "Version:       $version\n";
	exit;
}

sub writeLog {
	open LOG, ">$outputprefix".".layout_log";
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

sub graphStats {
	my %incount;
	my %outcount;
	my %connectioncount;
	my %seeds;
	my @length;
	foreach my $k1 (keys %adjacency) {
		$k1 =~ /^(.+)([io])$/;
		$seeds{$1} = 1;
		my $k1con = connectionCount($k1);
		if ($2 eq "i") { $incount{$k1con}++; $connectioncount{$k1con}++; }
		if ($2 eq "o") { $outcount{$k1con}++; $connectioncount{$k1con}++; }
		foreach my $k2 (keys %{$adjacency{$k1}}) {
			$length[0] += $adjacency{$k1}{$k2}[2];
			$length[1] += $adjacency{$k1}{$k2}[3];
			$length[2] += $adjacency{$k1}{$k2}[4];
			}
		}
	printLog(scalar keys %seeds, " seeds in the graph");	
	printLog("Total length between seeds = ".int(0.5+$length[1]/2)." nt \n");
	return ();
}

sub writeGraph { # 1 ignore;
	if ($_[1] == 1) { return(); }
	my $prefix = $_[0];
	my $output = $prefix.".graph";
	my $output2 = $prefix.".graph_tmp";
	open OUTPUT, ">$output";
	printLog("Writing graph to $output and $output2");
	foreach my $k1 (sort keys %adjacency) {
		foreach my $k2 (sort keys %{$adjacency{$k1}}) {
			$k1 =~ /^(.+)([io])$/; my @kk1 = ($1,$2);
			$k2 =~ /^(.+)([io])$/; my @kk2 = ($1,$2);
			next if ($kk2[0] lt $kk1[0]);
			print OUTPUT "$kk1[0]\t$kk1[1]$kk2[1]\t$kk2[0]\t",join("\t",@{$adjacency{$k1}{$k2}}),"\n";
			}
		}
	close OUTPUT;
	store \%adjacency, $output2;
	return();
}

sub writeScaffoldGraph { # 1 ignore;
	if ($_[1] == 1) { return(); }
	my $prefix = $_[0];
	my $output = $prefix.".graph";
	open OUTPUT, ">$output";
	printLog("Writing scaffold graph to $output");
	foreach my $k1 (sort keys %scafadjacency) {
		foreach my $k2 (sort keys %{$scafadjacency{$k1}}) {
			$k1 =~ /^(.+)([io])$/; my @kk1 = ($1,$2);
			$k2 =~ /^(.+)([io])$/; my @kk2 = ($1,$2);
			next if ($kk2[0] lt $kk1[0]);
			print OUTPUT "$kk1[0]\t$kk1[1]$kk2[1]\t$kk2[0]\t",join("\t",@{$scafadjacency{$k1}{$k2}}),"\n";
			}
		}
	close OUTPUT;
	return();
}

sub writeScaffolds {
	if ($_[0] == 1) { return(); }
	my $scafoutput = $outputprefix.".scaffolds";
	my $scafstats = $outputprefix.".scaffolds_stats";
	open SCOUTPUT, ">$scafoutput";
	open SCOUTPUTSTATS, ">$scafstats";
	printLog("Writing scaffold information to $scafoutput and $scafstats");
	foreach my $k1 (sort keys %scaffold) {
#		print SCOUTPUT "$k1\t",join(" ",@{$scaffold{$k1}}),"\n";
		print SCOUTPUT "$k1\t";
#		if (exists $terminal{node($scaffold{$k1}[1])}) { print SCOUTPUT "*"; }
		print SCOUTPUT join(" ",@{$scaffold{$k1}});
#		if (exists $terminal{node($scaffold{$k1}[-1])}) { print SCOUTPUT "*"; }
		print SCOUTPUT "\n";
		print SCOUTPUTSTATS "$k1\t$scaffold{$k1}[0]\t",sprintf("%.0f",$scaffoldlen{$k1}),"\n";
		}
	close SCOUTPUT;	
	my $output2 = $outputprefix.".scaffolds_tmp";
	my $output3 = $outputprefix.".seeds_tmp";
	printLog("Writing temporary scaffold files to $output2 and $output3");	
	store \%scaffold, $output2;	
	store \%seeduse, $output3;		
	return();
}


#### MAIN ALGORITHM SUBROUTINES

sub determineSeedLength {
	if (not $seedlength) {
		open (FASTA, "<$seedsfiles[0]\n") or die "Error opening $_\n";
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

sub initGraph {
	%adjacency = ();
	%linklength = ();
	@{$adjacency{"0i"}{"0i"}} = (0,0,0,0,0);
	$blacklist{0}++;
	return();
}	

sub deleteNullNode {
	my @adjlist = keys %{$adjacency{"0i"}};
	while (@adjlist) {
		my $del = shift @adjlist;
		deleteLink("0i",$del);
		}
	return();
}

sub readGraph { # format: 1000	io	999	15	1	908	972.466666666667	1059	61	ok
	open INPUT, "<$debuggraph";
	while (<INPUT>) {
		chomp;
		my @x = split /\t/,$_;
		$x[0] .= substr($x[1],0,1);
		$x[2] .= substr($x[1],1,1);
		my $a1 = $x[0];
		my $a2 = $x[2];
		shift @x; shift @x; shift @x;
		push @{$adjacency{$a1}{$a2}},@x;
		push @{$adjacency{$a2}{$a1}},@x;
		}
	close INPUT;
	return();
}
	
sub constructGraph {
	my $readmode = shift @_; # 1: initial, add all adjacent alignments as links; 2: final, ignore blacklisted and non-terminal nodes
	my %seedcount; my $readcount = 1; my $usedcount = 0; my $mode2count = 0; my $readignorecount = 0;
	my $currentread = "";
	my $threshold = $seedlength*$alignmentquality; my $alignok;
	my $linksadded = 0; my $newlinksadded = 0; my $previouslyadded = 0;
	my %stack;
	my @filesstack;
	if ($configfile) {
		open CONF, "<$configfile";
		while (<CONF>) {
			my @cnf = split /[\s\t]+/,$_; # dalinger alignmentfile fasta
			if (scalar @cnf != 3) { die "Invalid configuration file: $_\n"; }
			if ($cnf[0] ne "daligner" && $cnf[0] ne "coords" && $cnf[0] ne "blasr" && $cnf[0] ne "sam") { die "Invalid configuration file: $_\n"; }
			push @filesstack, "$cnf[0]\t$cnf[1]";
			}
		close CONF;	
		}
	else {	
		foreach (@mummercoords) { push @filesstack, "coords\t$_"; }
		foreach (@blasr) { push @filesstack, "blasr\t$_"; }
		foreach (@daligner) { push @filesstack, "daligner\t$_"; }
		foreach (@sam) { push @filesstack, "sam\t$_"; }
		}
	my $filecounter = 0;
	foreach (@filesstack) {
		my @openfile = split /\t/,$_;
		$filecounter++;
		printLog("Reading alignment coordinates from $openfile[1]");
		open (COORDS, "<$openfile[1]") or die "Error opening $openfile[1]\n";
		seek (COORDS,0,0); # not sure if necessary, reset to first byte on second opening
		if ($openfile[0] eq "daligner") { my $devnull = <COORDS>; $devnull = <COORDS>; }
		$currentread = "";
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
			next if ($data[7] == 0); # reserved node name
			next if (exists $blacklist{$data[7]});
			if ($data[4] < $threshold) { $alignok = 0; } else { $alignok = 1; }
			@data = extrapolateSeed(@data);
			if ($alignmentidentity && $data[6] < $alignmentidentity) { $alignok = 0; }
			if ($currentread eq "") { $currentread = $data[8]; } # first line in alignment file
			my @positions = sort {$a <=> $b} ($data[2],$data[3]);
			my $orientation;
			if ($positions[0] eq $data[2]) { $orientation = "f"; } else { $orientation = "r"; }
			if ($data[8] eq $currentread) {
				if (exists $stack{$positions[0]}) {	# same alignment position already occupied
					if ($stack{$positions[0]}[4] ne $data[7]) {
						}
					next;
					}
				if ($alignok) {	push @{$stack{$positions[0]}}, ($data[0],$data[1],$positions[1],$orientation,$data[7],$data[8]); } # seedstart seedend position orientation seedname readname
				next;
				}
			else {
				($usedcount,$mode2count,$readignorecount,$newlinksadded,$linksadded) = addLinks($readmode,\%seedcount,\%stack,$usedcount,$mode2count,$readignorecount,$filecounter,$newlinksadded,$linksadded);
				$currentread = $data[8];
				$readcount++;
				%stack = ();
				if ($alignok) { push @{$stack{$positions[0]}}, ($data[0],$data[1],$positions[1],$orientation,$data[7],$data[8]); } # seedstart seedend position orientation seedname readname
				} # wrap up read 
			}
		($usedcount,$mode2count,$readignorecount,$newlinksadded,$linksadded) = addLinks($readmode,\%seedcount,\%stack,$usedcount,$mode2count,$readignorecount,$filecounter,$newlinksadded,$linksadded);
		close COORDS;
		} # each file
	if ($readmode == 1) {
		printLog("$readcount long reads, ",scalar keys %seedcount," seeds, ",sprintf("%.2f",$usedcount/$readcount)," seeds per read");
		if (not $newlinksadded) { return (0); }
		return($linksadded/$newlinksadded); # return coverage
		}
	elsif ($readmode == 2) {
		printLog("$mode2count new links added, ",scalar keys %blacklist," seeds ignored");
		return();
		}
	return();
}

sub calculateLinkDistances {
	print "Calculating link distance metrics\n"; # currently the standard deviation
	foreach my $k1 (keys %linklength) {
		foreach my $k2 (keys %{$linklength{$k1}}) {
			my @ll = sort {$a <=> $b} @{$linklength{$k1}{$k2}};
			my $stdev = 0; my $squares = 0;
			if (scalar @ll > 1) {
				foreach(@ll) {
					$squares += ($adjacency{$k1}{$k2}[3]-$_)**2;
					}
				$stdev = ($squares/(scalar @ll))**0.5;
				}
			setAdjacencyField($k1,$k2,5,$stdev);
			}
		}
	%linklength = ();
	@{$adjacency{"0i"}{"0i"}} = (0,0,0,0,0,0);
	return();
}

sub removeSuspectLinks {
	my $delmode = $_[0]; # "nodes", "scaffolds"
	my $graphref = \%adjacency;
	if ($delmode eq "scaffolds") { $graphref = \%scafadjacency; }
	my @deleteme;
	my $lowcov = 0; my $ratio = 0; my $interval = 0; my $connect = 0; my $self = 0; # my $overlapping = 0;
	print "Finding suspect links\n";
	foreach my $first (sort keys %{$graphref}) {
		foreach my $second (sort keys %{$$graphref{$first}}) {
			next if (node($first) < node($second));
			next if ($first eq "0i");
			next if ($second eq "0i");
			if ($$graphref{$first}{$second}[0] < $coverage) {
				push @deleteme, ($first,$second);
				$lowcov++;
				next;
				}
			#if ($$graphref{$first}{$second}[4]-$$graphref{$first}{$second}[2] > $lengthinterval) { # arbitrary cutoff?
			if ($$graphref{$first}{$second}[5] > $lengthinterval) { # standard deviation of the gap length should be modest
				if (not ($diploid && $delmode eq "scaffolds"&& connectionCountS($first) == 1 && connectionCountS($second) == 1)) { # not if unambiguous
					push @deleteme, ($first,$second);
					$interval++;
					next;
					}
				}
			if (exists $$graphref{$first}{opposite($second)}) {
				if ($$graphref{$first}{opposite($second)}[0]*$coverageratio >= $$graphref{$first}{$second}[0]) { # arbitrary but necessary, it happens!
					push @deleteme, ($first,$second);
					$ratio++;
					next;
					}
				}
			if (exists $$graphref{opposite($first)}) {
				if (exists $$graphref{opposite($first)}{opposite($second)}) {
					if ($$graphref{opposite($first)}{opposite($second)}[0]*$coverageratio >= $$graphref{$first}{$second}[0]) {
						push @deleteme, ($first, $second);
						$ratio++;
						next;
						}
					}
				}
			if (node($first) eq node($second)) {
				push @deleteme, ($first, $second);
				$self++;
				}
			}
		}
	printLog("Deleting links: $lowcov for low coverage, $ratio for existing alternative connection, $interval for uncertain span, $self for self-linking");  
	while (@deleteme) {
		my $first = shift @deleteme; my $second = shift @deleteme;
		if ($delmode eq "nodes") {
			setAdjacencyField($first,$second,1,2);
			}
		else {
			setAdjacencyFieldS($first,$second,1,2);
			}
		}
	removeCollapsedLinks(0);
	return();
}

sub determineCoverage {
	my @evidence;
	my $sample = 0;
	my $maxsample = 10000; # do not look beyond 10000 random adjacencies
	K1LOOP: foreach my $k1 (keys %adjacency) {
		foreach my $k2 (keys %{$adjacency{$k1}}) {
			push @evidence, $adjacency{$k1}{$k2}[0];
			$sample++;
			if ($sample >= $maxsample) { last K1LOOP; }
			}
		}
	$sample = int($#evidence*$coveragepercentile);
	@evidence = sort {$a <=> $b} @evidence;
	$mediancoverage = sprintf("%.0f",$evidence[$sample]);
	#print "Median coverage $mediancoverage\n";
	return($mediancoverage);
}

sub removeProbableRepeatNodes {
	my %delnomination;
	print "Using high coverage to identify repeat nodes\n";
	foreach my $k1 (keys %adjacency) {
		foreach my $k2 (keys %{$adjacency{$k1}}) {
			$delnomination{$k1} += $adjacency{$k1}{$k2}[0];
			$delnomination{$k2} += $adjacency{$k1}{$k2}[0];
			}
		}
	my $actuallydeleted = 0;	
	foreach (keys %delnomination) {
		next if ($delnomination{$_} <= $highcoverage*$mediancoverage);
		$actuallydeleted += deleteNode(node($_));
		$blacklist{node($_)}++;
		}
	printLog("Removed $actuallydeleted putative repeat nodes based on high coverage\n");	
	return();
}

sub simplifyGraph {
	# Visits all adjacencies (several times) to see if any ambiguous connections can be explained away. Start easy and conservative, and finally end by aggressively
	# deleting remaining ambiguous connections. The end result is a completely linearized graph (but also fragmented). 
	# Two modes, each with a cautious and an aggressive pass
	my $simplifymode = shift @_; # 1: initial; 2: final
	my $aggression = "cautious"; # second pass is aggressive
	my $needscleanup = 0;
	my $repeatstatus = 1;
	print "Removing superfluous connections\n";
	my %simplifycounter;
	%loopprotect = (); # merge shortcut
	my $sortmergestack = 1; 
	@{$nodestack{"tip"}} = (); @{$nodestack{"longtip"}} = (); @{$nodestack{"simplebubble"}} = (); @{$nodestack{"bubble"}} = (); @{$nodestack{"connection"}} = (); @{$nodestack{"repeat"}} = (); @{$nodestack{"realrepeat"}} = (); @{$nodestack{"mayberepeat"}} = (); @{$nodestack{"loop"}} = (); @{$nodestack{"decision"}} = ();
	my %adjlength;
	foreach my $nodei (keys %adjacency) {
		foreach my $nodej (keys %{$adjacency{$nodei}}) {
			if (exists $adjlength{$nodei}) {
				if ($adjlength{$nodei} > $adjacency{$nodei}{$nodej}[3]) {
					$adjlength{$nodei} = $adjacency{$nodei}{$nodej}[3];
					}
				}
			else { $adjlength{$nodei} = $adjacency{$nodei}{$nodej}[3]; }
			}
		}
	@{$nodestack{"collapse"}} = sort {$adjlength{$a} <=> $adjlength{$b}} keys %adjlength; # shortest adjacencies first. sorting is clumsy. Plus it introduces differences between runs...
	MAINLOOP: while (@{$nodestack{"collapse"}} || @{$nodestack{"bubble"}} || @{$nodestack{"tip"}} || @{$nodestack{"longtip"}} || @{$nodestack{"connection"}} || @{$nodestack{"repeat"}} || @{$nodestack{"simplebubble"}} || @{$nodestack{"realrepeat"}} || @{$nodestack{"mayberepeat"}} || @{$nodestack{"loop"}} || @{$nodestack{"decision"}}) {
		my $status;
		my $anode1; my $anode2; my $a1con = 0; my $a2con = 0;
		if (@{$nodestack{"collapse"}}) { $anode1 = shift @{$nodestack{"collapse"}}; $status = "Collapsing"; }
		elsif ($needscleanup) { $status = "Removing collapsed";  }
		elsif (@{$nodestack{"tip"}}) { $anode1 = shift @{$nodestack{"tip"}}; $status = "Trimming"; }
		elsif (@{$nodestack{"longtip"}}) { $anode1 = shift @{$nodestack{"longtip"}}; $status = "Trimming more"; }
		elsif (@{$nodestack{"simplebubble"}}) { $anode1 = shift @{$nodestack{"simplebubble"}}; $status = "Collapsing more"; }
		elsif (@{$nodestack{"connection"}}) { $anode1 = shift @{$nodestack{"connection"}}; $status = "Linearizing"; }
		elsif (@{$nodestack{"loop"}}) { $anode1 = shift @{$nodestack{"loop"}}; $status = "Removing loops"; }
		elsif (@{$nodestack{"bubble"}}) { $anode1 = shift @{$nodestack{"bubble"}}; $status = "Merging"; }
		elsif (@{$nodestack{"decision"}}) { $anode1 = shift @{$nodestack{"decision"}}; $status = "Deciding"; }
		elsif (@{$nodestack{"repeat"}}) { 
			$status = "Marking repeats";
			if ($aggression eq "cautious") { 
				print "Repeat classification: original ",scalar @{$nodestack{"repeat"}};
				@{$nodestack{"repeat"}} = keys %adjacency; # not-so-elegant fix because a few repeat nodes have apparently escaped!
				markRepeatNodes(\@{$nodestack{"repeat"}},\@{$nodestack{"realrepeat"}},\@{$nodestack{"tip"}}); # divide repeats in two priority classes
#				print ", real ",scalar @{$nodestack{"realrepeat"}},", maybe ",scalar @{$nodestack{"tip"}},"\n";
				}
			else {
				markRepeatNodes(\@{$nodestack{"repeat"}},\@{$nodestack{"realrepeat"}},\@{$nodestack{"mayberepeat"}}); # divide repeats in two priority classes
				}
			$aggression = "aggressive";
			@{$nodestack{"repeat"}} = ();
			next;
			}
		elsif (@{$nodestack{"realrepeat"}}) { 
			$anode1 = shift @{$nodestack{"realrepeat"}}; $status = "Deleting"; 
			}
		elsif (@{$nodestack{"mayberepeat"}}) { 
			if ($repeatstatus == 1 && $simplifymode == 2) { deleteNullNode; $repeatstatus++; } 	
			$anode1 = shift @{$nodestack{"mayberepeat"}}; 
			$status = "Deleting";  
			}
		if ($status ne "Removing collapsed") {
			$anode2 = opposite($anode1);
			$a1con = connectionCount($anode1);
			if (not $a1con) { next; }
			$a2con = connectionCount($anode2);
			if ($anode1 eq "0i") { next; }	
			}
		# actual action starts here	
		if ($status eq "Collapsing") { 
			$needscleanup += collapseLinks($anode1,$anode2,$a1con,$a2con); 
			push @{$nodestack{"tip"}}, $anode1; 
			}
		elsif ($status eq "Removing collapsed") { 
			$simplifycounter{"collapse"} += removeCollapsedLinks(0); 
			$needscleanup = 0; 
			}
		elsif ($status eq "Trimming") { 
			$simplifycounter{"trim"} += trimNode($anode1,$anode2,$a1con,$a2con,$aggression); 
			push @{$nodestack{"longtip"}}, $anode1;
			}
		elsif ($status eq "Trimming more") { 
			$simplifycounter{"trim"} += trimNodeLong($anode1,$anode2,$a1con,$a2con,$aggressivetrimming[0],$aggressivetrimming[1],$aggressivetrimming[2],$aggression); # path b, path c, search depth
			push @{$nodestack{"simplebubble"}}, $anode1;
			}
		elsif ($status eq "Collapsing more") { 
			$simplifycounter{"bubble"} += collapseSimpleBubble($anode1,$anode2,$a1con,$a2con,$aggression); 
			push @{$nodestack{"connection"}}, $anode1; 
			}
		elsif ($status eq "Linearizing") { 
			my @linearizationsuccess = linearizeNode($anode1,$anode2,$a1con,$a2con,$linearization);
			$simplifycounter{"connection"} += shift @linearizationsuccess;
			foreach (@linearizationsuccess) {
				push @{$nodestack{"collapse"}}, $_;
				}
			push @{$nodestack{"loop"}}, $anode1;
			}
		elsif ($status eq "Removing loops") {
			$simplifycounter{"loop"} += collapseLoop($anode1,$anode2,$a1con,$a2con,$loopparameters[0],$loopparameters[1]);	# branching and length
			push @{$nodestack{"bubble"}}, $anode1;
			}
		elsif ($status eq "Merging") {
			my @bubblesuccess = collapseBubble($anode1,$anode2,$a1con,$a2con,$pathbranchingandlength{$simplifymode}{$aggression}[0],$pathbranchingandlength{$simplifymode}{$aggression}[1],$simplifymode,$diploid); 
			if ($bubblesuccess[0] == -1) { # identified a node that needs to be tried asap
				unshift @{$nodestack{"bubble"}}, $bubblesuccess[1];
				$loopprotect{$bubblesuccess[1]} = 1;
				}
			elsif ($bubblesuccess[0] > 0) { # something happened
				$simplifycounter{"bubble"} += shift @bubblesuccess;
				push @{$nodestack{"decision"}}, $anode1;
				push @{$nodestack{"collapse"}}, @bubblesuccess;
				}
			else { # try again or give up?
				push @{$nodestack{"decision"}}, $anode1;
				}
			}
		elsif ($status eq "Deciding") {	
			my $lowcoverage = 1; my $cratio = $mediancoverage; # minor path at most 1 evidence, major path median (= lower than typical)
			if ($simplifymode == 2) { $lowcoverage = $mediancoverage/$evidenceratio{$simplifymode}; }
			$cratio = $evidenceratio{$simplifymode}; # 25/75 path evidence if $highcoverage == 3
			$simplifycounter{"decision"} += decideBranchByEvidence($anode1,$anode2,$a1con,$a2con,$lowcoverage,$cratio); 
			if ($aggression eq "cautious") { push @{$nodestack{"repeat"}}, $anode1; } else { push @{$nodestack{"mayberepeat"}}, $anode1; } 
			}
		elsif ($status eq "Deleting") {
			my @nbnodes;
			if ($a1con) {
				foreach (keys %{$adjacency{$anode1}}) { push @nbnodes, ($_,opposite($_)); } # second chance for the neighbouring nodes
				}
			if ($a2con) {
				foreach (keys %{$adjacency{$anode2}}) { push @nbnodes, ($_,opposite($_)); }
				}
			my $repeatsuccess = deleteRepeatNode($anode1,$anode2,$a1con,$a2con,"check");	
			$simplifycounter{"repeat"} += $repeatsuccess;
			if ($repeatsuccess) { push @{$nodestack{"tip"}}, @nbnodes; }
			}
		} # MAINLOOP
	my $totalsimplification = 0;
	foreach (values %simplifycounter) {
		$totalsimplification += $_;
		}
	printLog("$totalsimplification graph simplifications");
#	foreach (sort keys %simplifycounter) { print "$_\t$simplifycounter{$_}\n"; }
	return();
}		

sub markTerminalNodes {
	print "Indexing terminal nodes... ";
	%terminal = ();
	foreach my $candidate (keys %adjacency) {
		my $anode = node($candidate);
		next if (exists $blacklist{$anode});
		next if (exists $adjacency{opposite($candidate)});
		$terminal{$anode} = opposite($candidate);
		}
	print scalar keys %terminal, "\n";
	return ();
}

sub markSubterminalNodes {
	print "Indexing subterminal nodes... ";
	%subterminal = ();
	foreach my $startpoint (values %terminal) { # values are seed.orientation dead ends
		my $pathlength = 0;
		my $pathnodes = 0;
		my $currentnode = opposite($startpoint);
		next if ($currentnode eq "0i");
		FOLLOWPATH: while (exists $adjacency{$currentnode}) {
			my @nextnode = keys %{$adjacency{$currentnode}}; # should be only one (if this sub is called after simplification)
			last if (scalar @nextnode > 1);
			last if (exists $subterminal{$nextnode[0]});
			last if (exists $terminal{node($nextnode[0])}); # 
			last if (exists $subterminal{opposite($nextnode[0])}); #
			$pathlength += $seedlength + $adjacency{$currentnode}{$nextnode[0]}[2]; # minimum gap estimate
			$pathnodes++;
			if ($pathlength > $maxreadlength) { last FOLLOWPATH; }
			if ($pathnodes > $subterminaldepth) { last FOLLOWPATH; }
			@{$subterminal{$nextnode[0]}} = ($startpoint,$pathlength);
			$currentnode = opposite($nextnode[0]);
			}
		}
	print scalar keys %subterminal, "\n";
	return ();
}

sub extractScaffolds {
	my @sclen;
	my @scnodes;
	my $sctotallen = 0;
	my $sctotalnodes = 0;
	my $scaffoldcounter = 0;
	%scaffold = ();	
	print "Indexing scaffolds\n";
	foreach my $scafmode("linear","circular") {
		foreach my $startadj (sort keys %adjacency) {
			next if ($startadj eq "0i");
			if ($scafmode eq "linear") { # try all linear scaffolds first
				next if (exists $adjacency{opposite($startadj)}); # not a terminal
				}
			my @tmp = keys %{$adjacency{$startadj}};
			next if (scalar @{$adjacency{$startadj}{$tmp[0]}} >6); # already counted
			my $currentadj = $startadj; #      Sa  xNy  xLy
			$scaffoldcounter++;
			my $currentnodes = 1;
			my $currentlen = $seedlength;
			$seeduse{node($currentadj)} = $scaffoldcounter;
			push @{$scaffold{$scaffoldcounter}}, (opposite($currentadj),$currentadj);
			while (exists $adjacency{$currentadj}) {
				my @nextadj = keys %{$adjacency{$currentadj}}; # assume linearity, so only [0] exists
				last if (scalar @{$adjacency{$currentadj}{$nextadj[0]}} >6);
				setAdjacencyField($currentadj,$nextadj[0],6,$scaffoldcounter);
				$currentnodes++;
				$currentlen += $adjacency{$currentadj}{$nextadj[0]}[3] + $seedlength;
				$currentadj = opposite($nextadj[0]);
				push @{$scaffold{$scaffoldcounter}}, (opposite($currentadj),$currentadj);
				$seeduse{node($currentadj)} = $scaffoldcounter;
				}
			unshift @{$scaffold{$scaffoldcounter}}, $scafmode;
			push @sclen, $currentlen;
			push @scnodes, $currentnodes;
			$sctotallen += $currentlen;
			$scaffoldlen{$scaffoldcounter} = $currentlen;
			$sctotalnodes += $currentnodes;
			}
		}
	printLog("$scaffoldcounter scaffolds found");
	@sclen = sort {$a <=> $b} @sclen;
	my $cumulen = 0;
	while (@sclen) {
		my $xlen = shift @sclen;
		$cumulen += $xlen;
		if ($cumulen >= 0.5*$sctotallen) {
			printLog("Scaffold N50 = ".sprintf("%.0f",$xlen));
			last;
			}
		}
	@scnodes = sort {$a <=> $b} @scnodes;
	my $cumunodes = 0;
	while (@scnodes) {
		my $xn = shift @scnodes;
		$cumunodes += $xn;
		if ($cumunodes >= 0.5*$sctotalnodes) {
			printLog("Scaffold N50 seeds = $xn");
			last;
			}
		}
	printLog("Total scaffold length = ".sprintf("%.0f",$sctotallen)." nt containing $sctotalnodes seeds");
	return();
}

sub buildScaffoldGraph {
	print "Constructing the scaffold graph\n";
	foreach my $sc (keys %scaffold) {
		next if ($scaffold{$sc}[0] eq "circular");
		$scafend{$scaffold{$sc}[1]} = $sc."i"; 
		$scafend{$scaffold{$sc}[-1]} = $sc."o";
		}
	foreach my $sc (keys %scaffold) {	
		if (exists $adjacency{$scaffold{$sc}[1]}) {
			foreach my $target (keys %{$adjacency{$scaffold{$sc}[1]}}) {
				if (exists $scafend{$target}) {
					@{$scafadjacency{$scafend{$target}}{$sc."i"}} = @{$adjacency{$scaffold{$sc}[1]}{$target}};
					@{$scafadjacency{$sc."i"}{$scafend{$target}}} = @{$adjacency{$scaffold{$sc}[1]}{$target}};
					}
				}
			}
		if (exists $adjacency{$scaffold{$sc}[-1]}) {
			foreach my $target (keys %{$adjacency{$scaffold{$sc}[-1]}}) {
				if (exists $scafend{$target}) {
					@{$scafadjacency{$scafend{$target}}{$sc."o"}} = @{$adjacency{$scaffold{$sc}[-1]}{$target}};
					@{$scafadjacency{$sc."o"}{$scafend{$target}}} = @{$adjacency{$scaffold{$sc}[-1]}{$target}};
					}
				}
			}
		}
	print scalar keys %scafadjacency, " entries in the scaffold graph\n";
	return();
}	

sub simplifyScaffoldGraph {
	print "Simplifying the scaffold graph\n";
	my $delcounter = 0;
	foreach my $currentscaffold (sort keys %scaffold) {
		my $anode1 = $currentscaffold."i"; my $anode2 = $currentscaffold."o";
		my $a1con = connectionCountS($anode1); my $a2con = connectionCountS($anode2);
		if ($a1con) {
			collapseScaffoldLinks($anode1,$anode2,$a1con,$a2con);
			if (exists $scafadjacency{$anode1}{$anode2}) { setAdjacencyFieldS($anode1,$anode2,1,30); }
			}
		if ($a2con) {
			collapseScaffoldLinks($anode2,$anode1,$a2con,$a1con);
			if (exists $scafadjacency{$anode2}{$anode1}) { setAdjacencyFieldS($anode2,$anode1,1,30); }
			}
		}
	$delcounter = removeCollapsedLinks(1); # 1 salvages links
	printLog("Deleted $delcounter connections between scaffolds");
	return();
}

sub resetScaffolds {
	%scaffold = ();
	%seeduse = ();
	%scafend = ();
	%scafadjacency = ();
	foreach my $xnode (keys %adjacency) {
		foreach my $ynode (keys %{$adjacency{$xnode}}) {
			while (scalar @{$adjacency{$xnode}{$ynode}} >6) {
				pop @{$adjacency{$xnode}{$ynode}};
				}
			}
		}
	return ();
}

sub addSalvagedLinks {
	print "Adding salvaged links (",scalar keys %salvagedadjacency," available)\n";
	my $scounter = 0;
	foreach my $anode1 (keys %adjacency) {
		my $anode2 = opposite($anode1);
		next if (connectionCount($anode2));
		if (exists $salvagedadjacency{$anode2}) {
			if (scalar keys %{$salvagedadjacency{$anode2}} == 1) {
				my @bnode = keys %{$salvagedadjacency{$anode2}};
				next if (connectionCount($bnode[0]));
				next if (scalar keys %{$salvagedadjacency{$bnode[0]}} > 1);
				next if (not connectionCount(opposite($bnode[0]))); # otherwise we're adding dead nodes
				next if (exists $blacklist{node($bnode[0])});
#print "$anode2 $bnode[0] : ",join(" ",@{$salvagedadjacency{$anode2}{$bnode[0]}}),"\n";
				@{$adjacency{$anode2}{$bnode[0]}} = @{$salvagedadjacency{$anode2}{$bnode[0]}};
				@{$adjacency{$bnode[0]}{$anode2}} = @{$salvagedadjacency{$anode2}{$bnode[0]}};
				$scounter++;
				}
			}
		}
	printLog("$scounter salvaged links added");
	%salvagedadjacency = ();
	return();
}


#### GRAPH CONSTRUCTION SUBROUTINES

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

sub addLinks {
	(my $readmode, my $seedcountref, my $stackref, my $usedcount, my $mode2count, my $readignorecount, my $alignfile, my $discretelinks, my $alllinks) = @_;
	# @{$stack{$positions[0]}}  -> seedstart seedend position orientation seedname readname
	my @sortkeys = sort {$a <=> $b} keys %$stackref; 
	if ($readmode ==1 && $#sortkeys >1) {
		my $readlength = $sortkeys[-1] - $sortkeys[0];
		if ($readlength > $maxreadlength) { $maxreadlength = $readlength; }
		}
	my %newadjlist; 
	$usedcount += scalar @sortkeys;
	FIRST: for (my $i = 0; $i < $#sortkeys; $i++) {
		my $firstalign = $sortkeys[$i];
		next if (exists $blacklist{$$stackref{$firstalign}[4]});
		if ($readmode == 2) { 
			my $orientation1 = $$stackref{$firstalign}[3];
			if ($orientation1 eq "f") { $orientation1 = "o"; } else { $orientation1 = "i"; }
			if (not exists $terminal{$$stackref{$firstalign}[4]}) {
				if (not exists $subterminal{$$stackref{$firstalign}[4].$orientation1}) {
					next;
					}
				}	
			}
		my @secondlist;
		if ($readmode != 2) { push @secondlist, $sortkeys[$i+1]; } # only the next alignment
		elsif ($i+2 > $#sortkeys) { @secondlist = (); }
#		else { @secondlist = @sortkeys[$i+2..$#sortkeys]; } # all next alignments, but excluding the directly adjacent
		else { @secondlist = @sortkeys[$i+1..$#sortkeys]; } # all next alignments
		SECOND: foreach my $secondalign (@secondlist) {
			next if (exists $blacklist{$$stackref{$secondalign}[4]});
			next if ($firstalign eq $secondalign);
			if ($readmode == 2) {
				my $orientation2 = $$stackref{$secondalign}[3];
				if ($orientation2 eq "f") { $orientation2 = "i"; } else { $orientation2 = "o"; }			
				next if (not exists $terminal{$$stackref{$secondalign}[4]} && not exists $subterminal{$$stackref{$secondalign}[4].$orientation2}); 
				}
			my @newadj = parsePairs($firstalign,@{$$stackref{$firstalign}},$secondalign,@{$$stackref{$secondalign}}); # return: seed.or seed.or gap name from to
#			if ($readmode == 1 && ($newadj[2] < (-1*$lengtherror*$seedlength))) { # more or less overlapping seeds, ignore at least one
			if ($newadj[2] < (-1*$lengtherror*$seedlength)) { # more or less overlapping seeds, ignore at least one
				if (not (exists $adjacency{$newadj[0]} && exists $adjacency{$newadj[0]}{$newadj[1]})) { # never mind if the connection already exists, then it can't be really bad
					next;
					}
				} #overlapping seeds
			if ($readmode == 2) {
				my $pathdepth = 0;
				if (exists $subterminal{$newadj[0]}) { 
					$pathdepth += $subterminal{$newadj[0]}[1]; 
					$newadj[0] = $subterminal{$newadj[0]}[0]; # attach to the terminal instead of subterminal
					}
				if (exists $subterminal{$newadj[1]}) { 
					$pathdepth += $subterminal{$newadj[1]}[1]; 
					$newadj[1] = $subterminal{$newadj[1]}[0]; # attach to the terminal instead of subterminal
					}
				next if ($pathdepth > $newadj[2] + $seedlength - 1);
				$newadj[2] -= $pathdepth;
				# make sure the correct end of the terminal node is used
				if (exists $terminal{node($newadj[0])}) { next if ($newadj[0] ne $terminal{node($newadj[0])}); }
				if (exists $terminal{node($newadj[1])}) { next if ($newadj[1] ne $terminal{node($newadj[1])}); }
				$mode2count++;
				}
			@{$newadjlist{join(".",@newadj)}} = @newadj;
			} # second
		} # first $i
	foreach (sort keys %newadjlist) {
		if (node($newadjlist{$_}[0]) eq node($newadjlist{$_}[1])) { 
			$readignorecount++;
			return ($usedcount,$mode2count,$readignorecount,$discretelinks,$alllinks); # possible palindromic read
			}
		}
	foreach (sort keys %newadjlist) {
		my @newadj = @{$newadjlist{$_}};
		next if (exists $blacklist{node($newadj[0])});
		next if (exists $blacklist{node($newadj[1])});
		$$seedcountref{node($newadj[0])} = 1;
		$$seedcountref{node($newadj[1])} = 1;
		# adjacency:	[0] direct count [1] 0 hypothetical, 1 evidence [2] minimum gap [3] mean gap (based on direct counts) [4]	maximum gap [5] sd [6] scaffold nr
		if (exists $adjacency{$newadj[0]}) {
			if (exists $adjacency{$newadj[0]}{$newadj[1]}) {
				$alllinks++;
				if ($adjacency{$newadj[0]}{$newadj[1]}[0]) {
					setAdjacencyField($newadj[0],$newadj[1],3,(($adjacency{$newadj[0]}{$newadj[1]}[3]*$adjacency{$newadj[0]}{$newadj[1]}[0])+$newadj[2])/(1+$adjacency{$newadj[0]}{$newadj[1]}[0]));
					}
				else {
					setAdjacencyField($newadj[0],$newadj[1],3,$newadj[2]);
					}
				setAdjacencyField($newadj[0],$newadj[1],0,$adjacency{$newadj[0]}{$newadj[1]}[0]+1);
				if ($adjacency{$newadj[0]}{$newadj[1]}[2] > $newadj[2]) { setAdjacencyField($newadj[0],$newadj[1],2,$newadj[2]);  }
				if ($adjacency{$newadj[0]}{$newadj[1]}[4] < $newadj[2]) { setAdjacencyField($newadj[0],$newadj[1],4,$newadj[2]);  }
				}
			else { # new adjacency
				newAdjacency($newadj[0],$newadj[1],1,1,$newadj[2],$newadj[2],$newadj[2]);
				$alllinks++;
				$discretelinks++
				}
			}
		else { # new adjacency
			newAdjacency($newadj[0],$newadj[1],1,1,$newadj[2],$newadj[2],$newadj[2]);
			$alllinks++;
			$discretelinks++
			}
		my @sortkeys = sort {$a cmp $b} ($newadj[0],$newadj[1]);
		push @{$linklength{$sortkeys[0]}{$sortkeys[1]}}, $newadj[2]; # temporary storage
		} 	
	return ($usedcount,$mode2count,$readignorecount,$discretelinks,$alllinks);
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


#### PRIMARY SIMPLIFICATION SUBROUTINES

sub collapseScaffoldLinks {
	(my $anode1, my $anode2, my $a1con, my $a2con) = @_;
	my $dellist = 0;
	my @anb = sort {$scafadjacency{$anode1}{$a}[3] <=> $scafadjacency{$anode1}{$b}[3]} (keys %{$scafadjacency{$anode1}}); # shortest links first
	my $anode = node($anode1);
	BNODE: foreach my $bnode1 (@anb) {
		next if ($bnode1 eq "0i");
		my $bnode2 = opposite($bnode1);
		my $bnode = node($bnode1);
		next if ($anode == $bnode);
		CNODE: foreach my $cnode1 (@anb) { # shortest first
			my $cnode = node($cnode1);
			next if ($cnode == $bnode);
			next if ($cnode == $anode);
			next if ($cnode1 eq "0i");
			next if ($scafadjacency{$anode1}{$cnode1}[1] >= 2); # already targeted for deletion
			my $cnode2 = opposite($cnode1);
			if (exists $scafadjacency{$bnode2} && exists $scafadjacency{$bnode2}{$cnode1} && exists $scafadjacency{$anode1} && exists $scafadjacency{$anode1}{$cnode1} && exists $scafadjacency{$anode1}{$bnode1} ) { # A-B, A-C, B-C -> A-B-C ?
				if (fitsBetweenS($anode1,$bnode1,$cnode1)) {
					my @bmultiple = (connectionCountS($bnode1),connectionCountS($bnode2));
					if ($bmultiple[0] == 1 && $bmultiple[1] == 1) { # really cautious, unambiguously connected? <- very slow!
						setAdjacencyFieldS($anode1,$cnode1,1,2); #$adjacency{$anode1}{$cnode1}[1] value 2 = targeted for deletion (default = 1)
						$dellist = 1;
						next CNODE;	
						} # unambiguously connected
					elsif ($bmultiple[0] <= 5 && $bmultiple[1] <= 5) { 	# not unambiguous, but if we find some additional evidence to place B we'll delete A-C anyway.
						my @comnb = commonNeighboursS($bnode2,$cnode2); # A-B-C-X
						foreach (@comnb) {
							if (fitsBetweenS($bnode2,$cnode1,$_)) {
								setAdjacencyFieldS($_,$bnode2,1,5);
								setAdjacencyFieldS($anode1,$cnode1,1,6);
								if (exists $scafadjacency{$anode1}{$_}) {
									setAdjacencyFieldS($_,$anode1,1,7);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursS($anode2,$bnode1); # X-A-B-C
						foreach (@comnb) {
							if (fitsBetweenS($_,$anode2,$bnode1)) {
								setAdjacencyFieldS($bnode1,$_,1,8);
								setAdjacencyFieldS($anode1,$cnode1,1,9);
								if (exists $scafadjacency{$cnode1}{$_}) {
									setAdjacencyFieldS($cnode1,$_,1,10);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursBetweenS($anode1,$bnode1); # returns X in orientation opposite $bnode1 (therefore $cnode1)
						foreach (@comnb) { # A-X-B-C
							if (fitsBetweenS($anode1,opposite($_),$bnode1)) {
								setAdjacencyFieldS($anode1,$bnode1,1,11);
								setAdjacencyFieldS($anode1,$cnode1,1,12);
								if (exists $scafadjacency{$cnode1}{$_}) {
									setAdjacencyFieldS($_,$cnode1,1,13);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursBetweenS($cnode1,$bnode2); # returns X in orientation opposite $bnode2 (therefore $anode1)
						foreach (@comnb) { # A-B-X-C
							if (fitsBetweenS($bnode2,$_,$cnode1)) {
								setAdjacencyFieldS($cnode1,$bnode2,1,14);
								setAdjacencyFieldS($anode1,$cnode1,1,15);
								if (exists $scafadjacency{$anode1}{$_}) {
									setAdjacencyFieldS($_,$anode1,1,16);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursS($anode1,$bnode2); # A-B-C...X	
						foreach (@comnb) {
							next if ($_ eq $cnode1);
							if (possibleFitS($bnode2,$cnode1,$_)) {
								setAdjacencyFieldS($anode1,$cnode1,1,17);
								#setAdjacencyField($anode1,$_,1,2);
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursS($cnode1,$bnode1); # X...A-B-C
						foreach (@comnb) {
							next if ($_ eq $anode1);
							if (possibleFitS($bnode1,$anode1,$_)) {
								setAdjacencyFieldS($anode1,$cnode1,1,18);
								#setAdjacencyField($cnode1,$_,1,2);
								$dellist = 1;
								next CNODE;
								}
							}
						# now the most complicated ones, with a (simple) bubble adjacent to the A-B-C triple
						# 1) A-C, A-B-C, B-X-Z-C
						# 2) A-C, A-B-C, B-X, C-Y-X
						# 3) A-C, A-B-C, B-X-Y, C-Y
						my @xnodes = sort keys %{$scafadjacency{$bnode2}};
						my @ynodes; if (exists $scafadjacency{$cnode2}) { @ynodes = sort keys %{$scafadjacency{$cnode2}}; }
						my @znodes = sort keys %{$scafadjacency{$cnode1}};
						foreach my $xnode1 (@xnodes) {
							my $xnode2 = opposite($xnode1);
							if (exists $scafadjacency{$xnode2}) {
								foreach my $znode2 (@znodes) {
									my $znode1 = opposite($znode2);
									if (exists $scafadjacency{$znode1}) {
										if (exists $scafadjacency{$xnode2}{$znode1}) { # B-X-Z-C confirmed, therefore A-B-C somewhat -> delete A-C
											setAdjacencyFieldS($anode1,$cnode1,1,19);
											$dellist = 1;
											next CNODE;
											}
										}
									}
								}
							foreach my $ynode1 (@ynodes) {
								my $ynode2 = opposite($ynode1);
								if (exists $scafadjacency{$ynode2}) {
									if (exists $scafadjacency{$xnode1}{$ynode2}) { # B-C-Y-X / B-X confirmed, therefore A-B-C somewhat -> delete A-C
										setAdjacencyFieldS($anode1,$cnode1,1,20);
										$dellist = 1;
										next CNODE;
										}
									}
								if (exists $scafadjacency{$xnode2}) {
									if (exists $scafadjacency{$ynode1}{$xnode2}) { # B-X-Y / C-Y confirmed, therefore A-B-C somewhat -> delete A-C
										setAdjacencyFieldS($anode1,$cnode1,1,21);
										$dellist = 1;
										next CNODE;
										}
									}
								}
							}
						# plus the symmetrical counterpart, with a bubble on the A-B side
						@xnodes = sort keys %{$scafadjacency{$bnode1}};
						if (exists $scafadjacency{$anode2}) { @ynodes = sort keys %{$scafadjacency{$anode2}}; } else { @ynodes = (); }
						@znodes = sort keys %{$scafadjacency{$anode1}};
						foreach my $xnode1 (@xnodes) {
							my $xnode2 = opposite($xnode1);
							if (exists $scafadjacency{$xnode2}) {
								foreach my $znode2 (@znodes) {
									my $znode1 = opposite($znode2);
									if (exists $scafadjacency{$znode1}) {
										if (exists $scafadjacency{$xnode2}{$znode1}) { 
											setAdjacencyFieldS($anode1,$cnode1,1,22);
											$dellist = 1;
											next CNODE;
											}
										}
									}
								}
							foreach my $ynode1 (@ynodes) {
								my $ynode2 = opposite($ynode1);
								my $xnode2 = opposite($xnode1);
								if (exists $adjacency{$ynode2}) {
									if (exists $scafadjacency{$xnode1}{$ynode2}) { 
										setAdjacencyFieldS($anode1,$cnode1,1,23);
										$dellist = 1;
										next CNODE;
										}
									}
								if (exists $scafadjacency{$xnode2}) {
									if (exists $scafadjacency{$ynode1}{$xnode2}) { 
										setAdjacencyFieldS($anode1,$cnode1,1,24);
										$dellist = 1;
										next CNODE;
										}
									}
								}
							}									
						} # B has additional connections
					} # fit
				} # A-B, A-C, B-C -> A-B-C ?
			} # cnode
		} # bnode
	return ($dellist);
}

sub collapseLinks {
	(my $anode1, my $anode2, my $a1con, my $a2con) = @_;
	my $dellist = 0;
	my @anb = sort {$adjacency{$anode1}{$a}[3] <=> $adjacency{$anode1}{$b}[3]} (keys %{$adjacency{$anode1}}); # shortest links first
	my $anode = node($anode1);
	BNODE: foreach my $bnode1 (@anb) {
		next if ($bnode1 eq "0i");
		my $bnode2 = opposite($bnode1);
		my $bnode = node($bnode1);
		next if ($anode == $bnode);
		CNODE: foreach my $cnode1 (@anb) { # shortest first
			my $cnode = node($cnode1);
			next if ($cnode == $bnode);
			next if ($cnode == $anode);
			next if ($cnode1 eq "0i");
			next if ($adjacency{$anode1}{$cnode1}[1] >= 2); # already targeted for deletion
			my $cnode2 = opposite($cnode1);
			if (exists $adjacency{$bnode2} && exists $adjacency{$bnode2}{$cnode1} && exists $adjacency{$anode1} && exists $adjacency{$anode1}{$cnode1} && exists $adjacency{$anode1}{$bnode1} ) { # A-B, A-C, B-C -> A-B-C ?
				if (fitsBetween($anode1,$bnode1,$cnode1)) {
					my @bmultiple = (connectionCount($bnode1),connectionCount($bnode2));
					if ($bmultiple[0] == 1 && $bmultiple[1] == 1) { # really cautious, unambiguously connected? <- very slow!
						if (not exists $blacklist{node($bnode1)}) {
							setAdjacencyField($anode1,$cnode1,1,2); #$adjacency{$anode1}{$cnode1}[1] value 2 = targeted for deletion (default = 1)
							$dellist = 1;
							}
						else { # B already blacklisted
							setAdjacencyField($anode1,$bnode1,1,3);
							setAdjacencyField($bnode2,$cnode1,1,4);
							$dellist = 1;
							}
						next CNODE;	
						} # unambiguously connected
					elsif ($bmultiple[0] <= 5 && $bmultiple[1] <= 5) { 	# not unambiguous, but if we find some additional evidence to place B we'll delete A-C anyway.
						my @comnb = commonNeighbours($bnode2,$cnode2); # A-B-C-X
						foreach (@comnb) {
							if (fitsBetween($bnode2,$cnode1,$_)) {
								setAdjacencyField($_,$bnode2,1,5);
								setAdjacencyField($anode1,$cnode1,1,6);
								if (exists $adjacency{$anode1}{$_}) {
									setAdjacencyField($_,$anode1,1,7);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighbours($anode2,$bnode1); # X-A-B-C
						foreach (@comnb) {
							if (fitsBetween($_,$anode2,$bnode1)) {
								setAdjacencyField($bnode1,$_,1,8);
								setAdjacencyField($anode1,$cnode1,1,9);
								if (exists $adjacency{$cnode1}{$_}) {
									setAdjacencyField($cnode1,$_,1,10);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursBetween($anode1,$bnode1); # returns X in orientation opposite $bnode1 (therefore $cnode1)
						foreach (@comnb) { # A-X-B-C
							if (fitsBetween($anode1,opposite($_),$bnode1)) {
								setAdjacencyField($anode1,$bnode1,1,11);
								setAdjacencyField($anode1,$cnode1,1,12);
								if (exists $adjacency{$cnode1}{$_}) {
									setAdjacencyField($_,$cnode1,1,13);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighboursBetween($cnode1,$bnode2); # returns X in orientation opposite $bnode2 (therefore $anode1)
						foreach (@comnb) { # A-B-X-C
							if (fitsBetween($bnode2,$_,$cnode1)) {
								setAdjacencyField($cnode1,$bnode2,1,14);
								setAdjacencyField($anode1,$cnode1,1,15);
								if (exists $adjacency{$anode1}{$_}) {
									setAdjacencyField($_,$anode1,1,16);
									}
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighbours($anode1,$bnode2); # A-B-C...X	
						foreach (@comnb) {
							next if ($_ eq $cnode1);
							if (possibleFit($bnode2,$cnode1,$_)) {
								setAdjacencyField($anode1,$cnode1,1,17);
								#setAdjacencyField($anode1,$_,1,2);
								$dellist = 1;
								next CNODE;
								}
							}
						@comnb = commonNeighbours($cnode1,$bnode1); # X...A-B-C
						foreach (@comnb) {
							next if ($_ eq $anode1);
							if (possibleFit($bnode1,$anode1,$_)) {
								setAdjacencyField($anode1,$cnode1,1,18);
								#setAdjacencyField($cnode1,$_,1,2);
								$dellist = 1;
								next CNODE;
								}
							}
						# now the most complicated ones, with a (simple) bubble adjacent to the A-B-C triple
						# 1) A-C, A-B-C, B-X-Z-C
						# 2) A-C, A-B-C, B-X, C-Y-X
						# 3) A-C, A-B-C, B-X-Y, C-Y
						my @xnodes = sort keys %{$adjacency{$bnode2}};
						my @ynodes; if (exists $adjacency{$cnode2}) { @ynodes = sort keys %{$adjacency{$cnode2}}; }
						my @znodes = sort keys %{$adjacency{$cnode1}};
						foreach my $xnode1 (@xnodes) {
							my $xnode2 = opposite($xnode1);
							if (exists $adjacency{$xnode2}) {
								foreach my $znode2 (@znodes) {
									my $znode1 = opposite($znode2);
									if (exists $adjacency{$znode1}) {
										if (exists $adjacency{$xnode2}{$znode1}) { # B-X-Z-C confirmed, therefore A-B-C somewhat -> delete A-C
											setAdjacencyField($anode1,$cnode1,1,19);
											$dellist = 1;
											next CNODE;
											}
										}
									}
								}
							foreach my $ynode1 (@ynodes) {
								my $ynode2 = opposite($ynode1);
								if (exists $adjacency{$ynode2}) {
									if (exists $adjacency{$xnode1}{$ynode2}) { # B-C-Y-X / B-X confirmed, therefore A-B-C somewhat -> delete A-C
										setAdjacencyField($anode1,$cnode1,1,20);
										$dellist = 1;
										next CNODE;
										}
									}
								if (exists $adjacency{$xnode2}) {
									if (exists $adjacency{$ynode1}{$xnode2}) { # B-X-Y / C-Y confirmed, therefore A-B-C somewhat -> delete A-C
										setAdjacencyField($anode1,$cnode1,1,21);
										$dellist = 1;
										next CNODE;
										}
									}
								}
							}
						# plus the symmetrical counterpart, with a bubble on the A-B side
						@xnodes = sort keys %{$adjacency{$bnode1}};
						if (exists $adjacency{$anode2}) { @ynodes = sort keys %{$adjacency{$anode2}}; } else { @ynodes = (); }
						@znodes = sort keys %{$adjacency{$anode1}};
						foreach my $xnode1 (@xnodes) {
							my $xnode2 = opposite($xnode1);
							if (exists $adjacency{$xnode2}) {
								foreach my $znode2 (@znodes) {
									my $znode1 = opposite($znode2);
									if (exists $adjacency{$znode1}) {
										if (exists $adjacency{$xnode2}{$znode1}) { 
											setAdjacencyField($anode1,$cnode1,1,22);
											$dellist = 1;
											next CNODE;
											}
										}
									}
								}
							foreach my $ynode1 (@ynodes) {
								my $ynode2 = opposite($ynode1);
								my $xnode2 = opposite($xnode1);
								if (exists $adjacency{$ynode2}) {
									if (exists $adjacency{$xnode1}{$ynode2}) { 
										setAdjacencyField($anode1,$cnode1,1,23);
										$dellist = 1;
										next CNODE;
										}
									}
								if (exists $adjacency{$xnode2}) {
									if (exists $adjacency{$ynode1}{$xnode2}) { 
										setAdjacencyField($anode1,$cnode1,1,24);
										$dellist = 1;
										next CNODE;
										}
									}
								}
							}									
						} # B has additional connections
					} # fit
				} # A-B, A-C, B-C -> A-B-C ?
			} # cnode
		} # bnode
	return ($dellist);
}

sub removeCollapsedLinks {
	my $savemode = $_[0]; # 0 save nothing; 1 store removed links in %salvagedadjacency
	my @fkeys = keys %adjacency;
	my $collapsecount = 0;
	foreach my $fnode (@fkeys) {
		next if (not exists $adjacency{$fnode});
		my @skeys = keys %{$adjacency{$fnode}};
		foreach my $snode (@skeys) {
			if (exists $adjacency{$fnode} && exists $adjacency{$fnode}{$snode}) {
				if ($adjacency{$fnode}{$snode}[1] >= 2) {
					my $tmp = $adjacency{$fnode}{$snode}[1];
					my @salvage = @{$adjacency{$fnode}{$snode}};
					if (deleteLink($fnode,$snode)) {
						$collapsecount++;
						if ($savemode == 1) {
							@{$salvagedadjacency{$fnode}{$snode}} = @salvage;
							@{$salvagedadjacency{$snode}{$fnode}} = @salvage;
							}
						}
					}
				}
			}
		}
	return ($collapsecount);
}

sub trimNode {
	(my $anode1, my $anode2, my $a1con, my $a2con, my $trimmode) = @_;
	my @anb = sort keys %{$adjacency{$anode1}};
	my $trimcount = 0;
	BNODE: foreach my $bnode1 (@anb) {
		next if ($bnode1 eq "0i");
		my $bnode2 = opposite($bnode1);
		foreach my $cnode1 (@anb) { 
			next if ($cnode1 eq $bnode1);
			my $cnode2 = opposite($cnode1);		
			if (exists $adjacency{$cnode2}) {
				if (connectionCount($bnode2) == 0)  { 
					if (connectionCount($bnode1) == 1) { # B is a tip, bye, clip
						deleteLink($anode1,$bnode1);
						$trimcount++;
						$blacklist{node($bnode1)}++;
						next BNODE;
						}
					}
				} # exists adj $cnode2
			elsif ($trimmode eq "aggressive") {
				if (connectionCount($bnode2) == 0 && connectionCount($bnode1) == 1) {
					if (connectionCount($cnode2) == 0 && connectionCount($cnode1) == 1) {
						if (exists $adjacency{$anode1}{$bnode1} && exists $adjacency{$anode1}{$cnode1}) {
							my $deletenode;
							if ($adjacency{$anode1}{$bnode1}[3] > $adjacency{$anode1}{$cnode1}[3]) { $deletenode = $cnode1; }
							else { $deletenode = $bnode1; }
							deleteLink($anode1,$deletenode);
							$trimcount++;
							$blacklist{node($deletenode)}++;
							next BNODE;
							}
						}
					}
				}
			}		
		} # B is a tip
	# another type of tip: A is itself a dead end, with maybe multiple incoming. Delete A, there's no way to resolve it otherwise. 	
	# only treat this in the second round, if it cannot be resolved as a bubble
	if ($trimmode eq "aggressive") { 
		if ($a1con >1) {
			if ($a2con == 0) {
				if (not exists $adjacency{$anode1}{"0i"}) {
					foreach (@anb) {
						deleteLink($anode1,$_);
						}
					$trimcount++;
					$blacklist{node($anode1)}++;
					}
				}
			}
		}
	return($trimcount);
}
			
sub trimNodeLong {			# needs updating to use generic sub followPath
	(my $anode1, my $anode2, my $a1con, my $a2con, my $trimb, my $trimc, my $pathdepth, my $trimmode) = @_; # @aggressivetrimming: (minor path max length, major path min length, search depth)
	my $trimcount = 0;
	my @anb = sort keys %{$adjacency{$anode1}}; 
	my @oppanb;
	if (exists $adjacency{$anode2}) { @oppanb = keys %{$adjacency{$anode2}}; }
	if (scalar @anb == 2 && scalar @oppanb == 1 && $anb[0] ne "0i" && $anb[1] ne "0i" && $oppanb[0] ne "0i") { 
		my @pathb = followPathToEnd(opposite($anb[0]),$pathdepth); # max|fork|end o($anb[0]) node1.or node1.or node2.or node2.or
		my @pathc = followPathToEnd(opposite($anb[1]),$pathdepth); # $aggressivetrimming[2] = search depth
		my $verdictb = shift @pathb; my $blength = 0.5*((scalar @pathb)+1); #pop @pathb; 
		my $verdictc = shift @pathc; my $clength = 0.5*((scalar @pathc)+1); #pop @pathc;
		if ($verdictb eq "end") {
			if ($blength <= $trimb) { 
				if ($clength > $trimc || $verdictc eq "fork" || $verdictc eq "max") { # clip path b
					deleteLink($anode1,$anb[0]);
					$trimcount++;
					deletePath(\@pathb);
					}
				}
			}
		elsif ($verdictc eq "end") {
			if ($clength <= $trimb) {
				if ($blength > $trimc || $verdictb eq "fork" || $verdictb eq "max") { # clip path c
					deleteLink($anode1,$anb[1]);
					$trimcount++;
					deletePath(\@pathc);
					}
				}
			}
		} 
	return ($trimcount);
}

sub collapseSimpleBubble {
	# Removes bubbles in which the bubble nodes (1+1, 2+0, or 2+1) have no additional branches
	(my $anode1, my $anode2, my $a1con, my $a2con, my $bubblemode) = @_;
	my $bubblecount = 0;
	my @anb = sort keys %{$adjacency{$anode1}};
	my $bnode1; my $bnode2; my $cnode1; my $cnode2; 
	my $b1con = 0; my $b2con = 0; my $c1con = 0; my $c2con = 0;
	if (scalar @anb != 2 && $bubblemode eq "cautious") { return(0); }
	for (my $i = 0; $i < $#anb; $i++) {
		for (my $j = $i+1; $j <= $#anb; $j++) {
			if (connectionCount($anb[$i]) >= connectionCount($anb[$j])) {
				$bnode1 = $anb[$i]; $bnode2 = opposite($bnode1);
				$cnode1 = $anb[$j]; $cnode2 = opposite($cnode1);
				}
			else {
				$bnode1 = $anb[$j]; $bnode2 = opposite($bnode1);
				$cnode1 = $anb[$i]; $cnode2 = opposite($cnode1);
				}
			$b1con = connectionCount($bnode1);
			$c1con = connectionCount($cnode1);
			$b2con = connectionCount($bnode2); 
			$c2con = connectionCount($cnode2); 
			if ($b1con == 1 && $b2con == 1 && $c1con == 1 && $c2con == 1) { # B or C is not the endpoint of the bubble
				my @comnb = commonNeighbours($bnode2,$cnode2); # A-B-X and A-C-X
				if (scalar @comnb == 1) { # delete the least evidence path
					if (($adjacency{$anode1}{$bnode1}[0] + $adjacency{$bnode2}{$comnb[0]}[0]) > ($adjacency{$anode1}{$cnode1}[0] + $adjacency{$cnode2}{$comnb[0]}[0])) {
						deleteLink($anode1,$cnode1);
						deleteLink($cnode2,$comnb[0]);
						$bubblecount++;
						$blacklist{node($cnode1)}++;
						}
					else {
						deleteLink($anode1,$bnode1);
						deleteLink($bnode2,$comnb[0]);
						$bubblecount++;
						$blacklist{node($bnode1)}++;
						}
					}
				else {	
					my @tmp = keys %{$adjacency{$bnode2}}; my $xnode1 = $tmp[0]; my $xnode2 = opposite($xnode1); my $x1con = connectionCount($xnode1); my $x2con = 0;
					$x2con = connectionCount($xnode2);
					@tmp = keys %{$adjacency{$cnode2}}; my $ynode1 = $tmp[0]; my $ynode2 = opposite($ynode1); my $y1con = connectionCount($ynode1); my $y2con = 0;
					$y2con = connectionCount($ynode2);
					if ($x1con == 2 && $y1con == 1 && $y1con == 1 && exists $adjacency{$xnode1}{$ynode2}) { # A-B-X and A-C-Y-X
						if (($adjacency{$anode1}{$bnode1}[0] + $adjacency{$bnode2}{$xnode1}[0]) > ($adjacency{$anode1}{$cnode1}[0] + $adjacency{$cnode2}{$ynode1}[0] + $adjacency{$ynode2}{$xnode1}[0])) { # delete A-C-Y-X
							deleteLink($anode1,$cnode1);
							deleteLink($cnode2,$ynode1);
							deleteLink($ynode2,$xnode1);
							$bubblecount++;
							$blacklist{node($cnode1)}++;
							$blacklist{node($ynode1)}++;
							}
						else {
							deleteLink($anode1,$bnode1);
							deleteLink($bnode2,$xnode1);
							$bubblecount++;
							$blacklist{node($bnode1)}++;
							}
						}
					elsif ($y1con == 2 && $x1con == 1 && $x2con == 1 && exists $adjacency{$xnode2}{$ynode1}) { # A-B-X-Y and A-C-Y
						if (($adjacency{$anode1}{$cnode1}[0] + $adjacency{$cnode2}{$ynode1}[0]) > ($adjacency{$anode1}{$bnode1}[0] + $adjacency{$bnode2}{$xnode1}[0] + $adjacency{$xnode2}{$ynode1}[0])) { # delete A-B-X-Y
							deleteLink($anode1,$bnode1);
							deleteLink($bnode2,$xnode1);
							deleteLink($xnode2,$ynode1);
							$bubblecount++;
							$blacklist{node($bnode1)}++;
							$blacklist{node($xnode1)}++;
							}
						else {
							deleteLink($anode1,$cnode1);
							deleteLink($cnode2,$ynode1);
							$bubblecount++;
							$blacklist{node($cnode1)}++;
							}	
						}
					}
				}
			elsif ($b1con == 2 && $c1con == 1 && $c2con == 1) { # B is the endpoint of the bubble (if multiple connections for B/C, B is defined as the one with most)
				if (exists $adjacency{$cnode2}{$bnode1}) { # simple A-B-C, should have collapsed before, so presumably no fit. Choose the highest evidence path
					if ($adjacency{$anode1}{$bnode1}[0] > ($adjacency{$anode1}{$cnode1}[0] + $adjacency{$cnode2}{$bnode1}[0])) {
						deleteLink($anode1,$cnode1);
						deleteLink($cnode2,$bnode1);
						$bubblecount++;
						$blacklist{node($cnode1)}++;
						}
					else {
						deleteLink($anode1,$bnode1);
						$bubblecount++;
						}
					}
				else {
					my @comnb = commonNeighboursBetween($cnode2,$bnode1); # returns X in orientation opposite $bnode1
					if (scalar @comnb == 1) {
						my $xnode1 = $comnb[0]; my $xnode2 = opposite($xnode1);
						if (connectionCount($xnode1) == 1 && connectionCount($xnode2) == 1) { # A-B and A-C-X-B, keep the longer path
							deleteLink($anode1,$bnode1);
							$bubblecount++;
							}
						}
					}
				}
			}
		} # scalar @anb = 2
	return ($bubblecount);
}

sub linearizeNode {
	(my $anode1, my $anode2, my $a1con, my $a2con, my $linearizationparameter) = @_;
	my $connectioncount = 0;
	my @returnstack;
	# Still more than one in, one out. If 2/1, it could be a spurious connection, try resolving that.
	if ($a1con == 2 && $a2con == 1) {
		# If A-B and A-C, with an additional incoming path X-B -> is A-B an artifact?
		# It cannot be a bubble anymore (well... if it misses the criteria)
		# If there is a long path ----A-C-----, delete A-B (and try something else with B)
		# Even better is there is also ----B-----
		my @anb = keys %{$adjacency{$anode1}};
		my $bnode1 = 0; my $cnode1 = 0; 
		if (connectionCount($anb[0]) == 2 && connectionCount($anb[1]) == 1) { $bnode1 = $anb[0]; $cnode1 = $anb[1]; }
		if (connectionCount($anb[0]) == 1 && connectionCount($anb[1]) == 2) { $bnode1 = $anb[1]; $cnode1 = $anb[0]; }
		if ($bnode1 && $cnode1) {
			my $bnode2 = opposite($bnode1);
			my $cnode2 = opposite($cnode1);
			my $xnode1;
			my @bnb = keys %{$adjacency{$bnode1}}; # already tested, there are 2: anode1 and xnode1
			if ($bnode1 eq "0i") { $xnode1 = "0i"; } # irrelevant placeholder
			elsif ($bnb[0] eq $anode1) { $xnode1 = $bnb[1]; }
			else { $xnode1 = $bnb[0]; }
			my $xnode2 = opposite($xnode1);
			if (exists $adjacency{$bnode2}) { 
				if (connectionCount($bnode2) == 1 || $bnode1 eq "0i") {
					if (exists $adjacency{$cnode2}) { # not testing existence of xnode2, not absolutely required
						my @a2paths = followPath($anode2,0,$linearizationparameter,"any"); # this should be one long, nonbranching path
						my @c2paths = followPath($cnode2,0,$linearizationparameter,"any"); # this should be another long, nonbranching path
						my @x2paths = followPath($xnode2,0,$linearizationparameter,"hub"); # just to get the first ambiguity to add to the stack again
						my @b2paths = followPath($bnode2,0,$linearizationparameter,"hub");
						if (scalar @a2paths == 1 && scalar @c2paths == 1) { # if unique paths from a2 & c2
							my @a2 = split /\,/,$a2paths[0]; # path lengths 3x, forks, joins, evidence, elements, last_type, 1o 2i 2o 3o 3i etc
							my @c2 = split /\,/,$c2paths[0];
							if ($a2[6] >= $linearizationparameter && $c2[6] >= $linearizationparameter) { # both long enough
								deleteLink($anode1,$bnode1);
								$connectioncount++;
								if ($bnode1 ne "0i") {
									push @returnstack, ($bnode1, $bnode2);
									foreach (@x2paths, @b2paths) {
										my @bx2 = split /\,/,$_;
										if ($bx2[7] ne "middle") { # the last node out from b2 or x2 needs resolution, maybe that is easier now that A-B is gone
											push @returnstack, ($bx2[-1], $bx2[-2]);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	unshift @returnstack, $connectioncount;
	return(@returnstack);
}

sub collapseLoop {
	(my $anode1, my $anode2, my $a1con, my $a2con, my $maxbranching, my $maxpathlength) = @_;
	my $bubblecount = 0;	
	# I observed a lot of bubble-like loops, in which A-B (minor evidence) but also A-...-...-...-...-B (most evidence)
	# in this case, snap the A-B link if the alternative path is not too long
	if ($a1con < 2) { return (0); }
	my %deleteloopstack; # direct connection between start and end
	my @pathstartnodes;
	my %loopendnodes;
	foreach my $bnode1 (keys %{$adjacency{$anode1}}) {
		next if (node($bnode1) == node($anode1));
		if (connectionCount($bnode1) >= 2) { $loopendnodes{$bnode1} = 1; }
		if (connectionCount(opposite($bnode1))) { push @pathstartnodes, opposite($bnode1); }
		}
	foreach my $pathstart (@pathstartnodes) {
		next if (exists $loopendnodes{opposite($pathstart)});
		my @newpaths = followPath($pathstart,$maxbranching,$maxpathlength,"hub");
		foreach (@newpaths) { 
			my @pathelements = split /\,/,$_;				
			if (exists $loopendnodes{$pathelements[-2]}) { $deleteloopstack{$pathelements[-2]} = 1; }
			}
		}
	foreach (keys %deleteloopstack) {
		$bubblecount++;
		deleteLink($anode1,$_);
		}
	return($bubblecount);
}

sub collapseBubble {
	(my $anode1, my $anode2, my $a1con, my $a2con, my $branching, my $plen, my $trymode, my $diploidpreserve) = @_;
	my $bubblecount = 0;	
	# After all easy simplification options have been exhausted, try following paths from $anode1 to the next (merging) fork node.
	# If the merge node is the same for two different paths, and both are ~the same length, delete the lowest evidence path
	# Alternatively ($diploid mode), clip off a path with no additional connections
	my @newpaths = followPath($anode1,$branching,$plen,"hub");
	my %endindexedpaths;
	foreach (@newpaths) { 
		my @pathelements = split /\,/,$_;
		next if ($pathelements[-1] eq $anode1); # initial fork reported back
		push @{$endindexedpaths{$pathelements[-1]}}, $_;
		}
	foreach (sort keys %endindexedpaths) {
		next if (scalar @{$endindexedpaths{$_}} <2);
		my @pathstack = sort { length($a) <=> length($b) } @{$endindexedpaths{$_}}; # shortest first
		while (@pathstack) {
			my $path1 = shift @pathstack;
			my @pathelements1 = split /\,/,$path1;
			foreach my $path2 (@pathstack) {
				next if ($path1 eq $path2);  
				my @pathelements2 = split /\,/,$path2; #
				next if ($pathelements1[-1] ne $pathelements2[-1]); # endpoints not equal [should not happen]
				my @commonnode = considerPaths(9,\@pathelements1,\@pathelements2); # returns first common + position in second path (default 0,0) # start at array pos 9 = the first after the initial forking node
				if ($commonnode[0]) {
					if ($commonnode[1] != $#pathelements2 -1) {
						if (not exists $loopprotect{opposite($commonnode[0])}) {
							return(-1,opposite($commonnode[0]));
							}
						next; 
						}
					} 
				next if ($pathelements1[2] < $lengtherror*$pathelements2[0]); # max < min
				next if ($pathelements2[2] < $lengtherror*$pathelements1[0]); # max < min
				$bubblecount++;
				my $whathappened = treatBubblePaths(\@pathelements1,\@pathelements2,$trymode,$diploidpreserve);
				if ($whathappened) {
					my @returnstack;
					foreach (keys %endindexedpaths) {
						push @returnstack, ($_, opposite($_));
						}
					#push @returnstack, ($pathelements1[-1], $pathelements1[-2]);
					push @returnstack, ($anode1, $anode2);
					return ($bubblecount, @returnstack);
					}
				}
			}
		}
	return (0);	
}

sub decideBranchByEvidence {
	(my $anode1, my $anode2, my $a1con, my $a2con, my $lowev, my $evratio) = @_;
	my @anb = sort keys %{$adjacency{$anode1}};
	my $trimcount = 0;
	my $bnode1; my $cnode1;
	if ($a1con == 2 && $a2con == 1) {
		if ($adjacency{$anode1}{$anb[0]}[0] >= $adjacency{$anode1}{$anb[1]}[0]) {
			$bnode1 = $anb[0]; $cnode1 = $anb[1];
			}
		else {
			$bnode1 = $anb[1]; $cnode1 = $anb[0];
			}
		if ($adjacency{$anode1}{$cnode1}[0] <= $lowev) {
			if ($adjacency{$anode1}{$bnode1}[0] >= $evratio*$adjacency{$anode1}{$cnode1}[0]) {
				if (connectionCount($bnode1) == 1 && connectionCount($cnode1) == 1) {
					if (connectionCount(opposite($bnode1)) == 1 && connectionCount(opposite($cnode1)) == 1) {
						deleteLink($anode1,$cnode1);
						$trimcount++;
						}
					}
				}
			}
		}
	return($trimcount);	
}

sub markRepeatNodes {	
	(my $stackref, my $realref, my $mayberef) = @_; # divide $stackref between the two classes
	while (@{$stackref}) {
		my $rnode1 = shift @{$stackref};
		my $rnode2 = opposite($rnode1);
		my $r1con = connectionCount($rnode1);
		my $r2con = connectionCount($rnode2);
		# maybe: one side multiple links, other side at most one - and not to a repeat
		if ($r1con >= 3 && $r2con >= 3) { 
			push @{$realref}, $rnode1;
			next;
			}
		if ($r1con == 2) {
			if ($r2con == 0) {
				push @{$mayberef}, $rnode1;
				next;
				}
			elsif ($r2con ==1) {
				my @r2nb = keys %{$adjacency{$rnode2}};
				if (connectionCount($r2nb[0]) <= 1 && connectionCount(opposite($r2nb[0])) <= 1 ) { 
					push @{$mayberef}, $rnode1;
					next;
					}
				}
			}
		if ($r2con == 2) {
			if ($r1con == 0) {
				push @{$mayberef}, $rnode2;
				next;
				}
			elsif ($r1con ==1) {
				my @r1nb = keys %{$adjacency{$rnode1}};
				if (connectionCount($r1nb[0]) <= 1 && connectionCount(opposite($r1nb[0])) <= 1 ) { 
					push @{$mayberef}, $rnode2;
					next;
					}
				}
			}
		if ($r1con > 1 || $r2con > 1) {	
			unshift @{$mayberef}, $rnode1; 	
			}
		}
	return (1);
}

sub deleteRepeatNode {
	(my $anode1, my $anode2, my $a1con, my $a2con, my $check) = @_;
	if ($check eq "check") {
		if (($a1con <= 1) && ($a2con <= 1)) { return (0); }
		}
	my @anb1; my @anb2;
	if ($a1con) { @anb1 = sort keys %{$adjacency{$anode1}}; }
	if ($a2con) { @anb2 = sort keys %{$adjacency{$anode2}}; }
	foreach (@anb1) {
		deleteLink($_,$anode1);
		}
	foreach (@anb2) {
		deleteLink($_,$anode2);
		}
	$blacklist{node($anode1)}++;
	return (1);
}


#### SECONDARY SIMPLIFICATION SUBROUTINES

sub treatBubblePaths {
	# choose what to do based on path characteristics
	# diploid assembly should be implemented somewhere here!
	my @patha; my @pathb;	# (min,mean,max,forks,joins,ev,nr,$from,...)
	@patha = @{$_[0]}; # \@{$pathtofork{$bnode1}{$i}}
	@pathb = @{$_[1]}; # \@{$pathtofork{$cnode1}{$j}}
	my $simplifymode = $_[2]; # 1: initial pass, delete simple paths only; 2: complex
	my $diploidpreserve = $_[3];
	if ($patha[3] == 0 && $patha[4] <= 1) { # no forks or joins # the initial fork is not counted; the final join is.
		if ($pathb[3] == 0 && $pathb[4] <= 1) { # no forks or joins
			if ($patha[6] >= $pathb[6]) { # choose the path that preserves the most seeds
				deletePath(\@pathb); # also delete all adjacencies within this path
				return("deleteb");
				}
			else {
				deletePath(\@patha);
				return("deletea");
				}
			}
		# joins in path b
		if ($simplifymode == 1) { return(0); }
		detachPath(\@pathb); # do not delete the adjacencies within the path
		return("detachb");
		}
	elsif ($pathb[3] == 0 && $pathb[4] <= 1) {	# no forks or joins in path b, but they exist in path a
		if ($simplifymode == 1) { return(0); }
		detachPath(\@patha);
		return("detacha");
		}
	# forking in both paths, detach everything
	if ($simplifymode == 1) { return(0); }
	detachPath(\@patha);
	detachPath(\@pathb);
	return("detach2");
}

sub detachPath {
	my @detachpath = @{$_[0]}; # (min,mean,max,forks,joins,ev,nr,lasttype,$from,...)
	if (exists $adjacency{$detachpath[8]} && exists $adjacency{$detachpath[8]}{$detachpath[9]}) { # could have been deleted already
		$adjacency{"0i"}{$detachpath[9]} = $adjacency{$detachpath[8]}{$detachpath[9]};
		$adjacency{$detachpath[9]}{"0i"} = $adjacency{"0i"}{$detachpath[9]};
		}
	if (exists $adjacency{$detachpath[-2]} && exists $adjacency{$detachpath[-2]}{$detachpath[-3]}) { # could have been deleted already	
		$adjacency{"0i"}{$detachpath[-2]} = $adjacency{$detachpath[-2]}{$detachpath[-3]};
		$adjacency{$detachpath[-2]}{"0i"} = $adjacency{"0i"}{$detachpath[-2]};
		}
	deleteLink($detachpath[8],$detachpath[9]);
	deleteLink($detachpath[-2],$detachpath[-3]); # the last connection ([-1] is the same node as [-2])
	for (my $i = 8; $i < $#detachpath-2; $i += 2) {
		$blacklist{node($detachpath[$i])}++; # do not try to reattach these; they are only kept for their connections
		}
	return();
}

sub deletePath {
	my @delpath = @{$_[0]};
	my $firstpos = 0;
	if ($delpath[0] !~ /[io]/) { $firstpos = 8; } # bubble path format
	for (my $i = $firstpos; $i <= $#delpath-1; $i += 2) {
		deleteLink($delpath[$i],$delpath[$i+1]);
		}
	for (my $i = $firstpos+1; $i < $#delpath-1; $i++) {
		deleteNode(node($delpath[$i]));
		$blacklist{node($delpath[$i])}++;
		}
	return();
}

sub considerPaths {
	my %path; 
	my $startat = $_[0];
	my @patha = @{$_[1]}; # \@{$pathtofork{$bnode1}{$i}}
	my @pathb = @{$_[2]}; # \@{$pathtofork{$cnode1}{$j}}
	my $firstcommon = 0; my $firstb = 0;
	for (my $i = $startat; $i < $#patha; $i++) { 
		if (exists $path{$patha[$i]}) { $firstcommon = $patha[$i]; }
		$path{$patha[$i]} = 1; 
		}
	for (my $i = $startat; $i < $#pathb; $i++) { 
		if (exists $path{$pathb[$i]}) { $firstcommon = $pathb[$i]; $firstb = $i; last;}
		$path{$pathb[$i]} = 1;
		}
	return ($firstcommon,$firstb);
}
	
sub followPath {
	(my $from, my $complexity, my $maxsteps, my $endpoint) = @_; # endpoints "fork" "join" "hub" (fork & join) "end" "any"  
	my $forkmode = 0; my $joinmode = 0; my $anymode = 0; my $endmode;
	if ($endpoint eq "join") { $joinmode = 1; }
	elsif ($endpoint eq "fork") { $forkmode = 1; }
	elsif ($endpoint eq "hub") { $joinmode = 1; $forkmode = 1;}
	elsif ($endpoint eq "end") { $endmode = 1; } # note that an endpoint can also be a joining node...
	else { $anymode = 1; }
	my %returnpaths; # indexing
	my @workingpathindex = (0);
	my %workingpaths;
	my $pathcounter = 0;
	@{$workingpaths{0}} = (0,0,0,0,0,0,0,"start",$from); # path lengths 3x, forks, joins, evidence, elements, last_type, 1o 2i 2o 3o 3i etc
	while (@workingpathindex) {
		my $currentpathindex = shift @workingpathindex;
		my @currentpath = @{$workingpaths{$currentpathindex}}; 
		if (exists $adjacency{$currentpath[-1]} && $currentpath[6] < $maxsteps && ($currentpath[3]+$currentpath[4]) <= $complexity) {
			my @nextnodes = sort keys %{$adjacency{$currentpath[-1]}};
			foreach my $nnode (@nextnodes) {
				next if ($nnode eq "0i");
				my @newpath = @currentpath;
				$pathcounter++;
				my $nnode2 = opposite($nnode);
				$newpath[0] += $adjacency{$currentpath[-1]}{$nnode}[2] + $seedlength;
				$newpath[1] += $adjacency{$currentpath[-1]}{$nnode}[3] + $seedlength;
				$newpath[2] += $adjacency{$currentpath[-1]}{$nnode}[4] + $seedlength;
				$newpath[5] += $adjacency{$currentpath[-1]}{$nnode}[0]; # evidence
				$newpath[6]++; # elements
				push @newpath, ($nnode,$nnode2);
				my $verdict = "middle";
				if (connectionCount($nnode) >1) { # join
					$newpath[4]++;
					$verdict = "join";
					if ($joinmode && ($newpath[3]+$newpath[4]) <= $complexity) { $returnpaths{$pathcounter} = 1; }
					}
				if (exists $adjacency{$nnode2}) {
					if (connectionCount($nnode2) >1) { # fork
						if ($verdict eq "join") { $verdict = "hub"; } else { $verdict = "fork"; }
						if ($forkmode && ($newpath[3]+$newpath[4]) <= $complexity) { $returnpaths{$pathcounter} = 1; }
						}
					}
				elsif ($verdict eq "middle") { 
					$verdict = "end"; 
					if ($endmode && ($newpath[3]+$newpath[4]) <= $complexity) { $returnpaths{$pathcounter} = 1; }
					}
				elsif ($verdict eq "join") { 
					$verdict = "endjoin"; 
					if ($forkmode && ($newpath[3]+$newpath[4]) <= $complexity) { $returnpaths{$pathcounter} = 1; }
					}
				if ($anymode && $newpath[6] >= $maxsteps && ($newpath[3]+$newpath[4]) <= $complexity) {	$returnpaths{$pathcounter} = 1;	}
				$newpath[7] = $verdict;
				push @workingpathindex, $pathcounter;	
				@{$workingpaths{$pathcounter}} = @newpath;
				} # foreach @nextnodes
			} # add nextnodes
		} # workingpathindex
	my @returnarray;
	foreach (sort keys %returnpaths) {
		push @returnarray, join(",",@{$workingpaths{$_}});
		}
	return (@returnarray);	
}

sub followPathToEnd {
	my $from = $_[0]; 
	my $pathcounter = 0;
	my $maxpathsteps = $_[1]; # all path content
	my @path = ($from);
	while (scalar @path < $maxpathsteps) {
		if (exists $adjacency{$path[-1]}) {
			my @next = sort keys %{$adjacency{$path[-1]}};
			if (scalar @next >1) { return ("fork",@path); }
			push @path, $next[0],opposite($next[0]);
			}
		else { return ("end",@path); }	
		}
	return ("max",@path);
}	

sub commonNeighbours {
	# returns an array of common neigbour nodes.orientation - excluding 0i - of $_[0] and $_[1]
	if (not exists $adjacency{$_[0]}) { return (); }
	if (not exists $adjacency{$_[1]}) { return (); }
	my @returnarray;
	if (node($_[0]) eq node($_[1])) { return (); }
	foreach (sort keys %{$adjacency{$_[0]}}) {
		next if ($_ eq "0i");
		next if (node($_) eq node($_[0]));
		next if (node($_) eq node($_[1]));
		if (exists $adjacency{$_[1]}{$_}) {
			push @returnarray, $_;
			}
		}
	return (@returnarray);
}

sub commonNeighboursS {
	# returns an array of common neigbour nodes.orientation - excluding 0i - of $_[0] and $_[1]
	if (not exists $scafadjacency{$_[0]}) { return (); }
	if (not exists $scafadjacency{$_[1]}) { return (); }
	my @returnarray;
	if (node($_[0]) eq node($_[1])) { return (); }
	foreach (sort keys %{$scafadjacency{$_[0]}}) {
		next if ($_ eq "0i");
		next if (node($_) eq node($_[0]));
		next if (node($_) eq node($_[1]));
		if (exists $scafadjacency{$_[1]}{$_}) {
			push @returnarray, $_;
			}
		}
	return (@returnarray);
}

sub commonNeighboursBetween {
	# returns an array of common neigbour nodes.orientation - excluding 0i - of $_[0] and $_[1]
	# return orientation: opposite $_[1]
	if (not exists $adjacency{$_[0]}) { return (); }
	if (not exists $adjacency{$_[1]}) { return (); }
	if (node($_[0]) eq node($_[1])) { return (); }
	my @returnarray;
	foreach (sort keys %{$adjacency{$_[0]}}) {
		if (exists $adjacency{$_[1]}{opposite($_)}) {
		next if ($_ eq "0i");
			next if (node($_) eq node($_[0]));
			next if (node($_) eq node($_[1]));
			push @returnarray, opposite($_);
			}
		}
	return (@returnarray);
}

sub commonNeighboursBetweenS {
	# returns an array of common neigbour nodes.orientation - excluding 0i - of $_[0] and $_[1]
	# return orientation: opposite $_[1]
	if (not exists $scafadjacency{$_[0]}) { return (); }
	if (not exists $scafadjacency{$_[1]}) { return (); }
	if (node($_[0]) eq node($_[1])) { return (); }
	my @returnarray;
	foreach (sort keys %{$scafadjacency{$_[0]}}) {
		if (exists $scafadjacency{$_[1]}{opposite($_)}) {
		next if ($_ eq "0i");
			next if (node($_) eq node($_[0]));
			next if (node($_) eq node($_[1]));
			push @returnarray, opposite($_);
			}
		}
	return (@returnarray);
}

sub fitsBetween {
	# returns 1 if A-B and opp(B)-C fit in A-C (i.e. A-B-C is consistent)
	my $aa = $_[0]; my $bb = $_[1]; my $cc = $_[2]; my $bb2 = opposite($bb);
	# does not test if adjacencies exist...!
	if (($adjacency{$aa}{$bb}[2]+$adjacency{$bb2}{$cc}[2]+$seedlength)*$lengtherror < $adjacency{$aa}{$cc}[4]) { # [2] = min
		if (($adjacency{$aa}{$bb}[4]+$adjacency{$bb2}{$cc}[4]+$seedlength) > $lengtherror*$adjacency{$aa}{$cc}[2]) { # [4] = max
			return (1);
			}
		}
	return (0);
}

sub fitsBetweenS {
	# returns 1 if A-B and opp(B)-C fit in A-C (i.e. A-B-C is consistent)
	my $aa = $_[0]; my $bb = $_[1]; my $cc = $_[2]; my $bb2 = opposite($bb);
	# does not test if adjacencies exist...!
	if (($scafadjacency{$aa}{$bb}[2]+$scafadjacency{$bb2}{$cc}[2]+$scaffoldlen{node($bb)})*$lengtherror < $scafadjacency{$aa}{$cc}[4]) { # [2] = min
		if (($scafadjacency{$aa}{$bb}[4]+$scafadjacency{$bb2}{$cc}[4]+$scaffoldlen{node($bb)}) > $lengtherror*$scafadjacency{$aa}{$cc}[2]) { # [4] = max
			return (1);
			}
		}
	return (0);
}

sub possibleFit {
	# returns 1 if A-B and opp(B)-C fit in A-C (i.e. A-B-C is consistent)
	# does not assume the direct existence of opp(B)-C
	my $aa = $_[0]; my $bb = $_[1]; my $cc = $_[2]; 
	if (($adjacency{$aa}{$bb}[2]+$seedlength)*$lengtherror < $adjacency{$aa}{$cc}[4]) { # [2] = min
		return (1);
		}
	return (0);
}	
	
sub possibleFitS {
	# returns 1 if A-B and opp(B)-C fit in A-C (i.e. A-B-C is consistent)
	# does not assume the direct existence of opp(B)-C
	my $aa = $_[0]; my $bb = $_[1]; my $cc = $_[2]; 
	if (($scafadjacency{$aa}{$bb}[2]+$scaffoldlen{node($bb)})*$lengtherror < $scafadjacency{$aa}{$cc}[4]) { # [2] = min
		return (1);
		}
	return (0);
}

	
#### GRAPH MANIPULATION & GENERIC SUBROUTINES

sub deleteLink {
	my $rval = 0;
	if (exists $adjacency{$_[0]}{$_[1]}) { $rval = 1; }
	$adjacency{$_[0]}{$_[1]} = ();
	$adjacency{$_[1]}{$_[0]} = ();
	delete ($adjacency{$_[0]}{$_[1]});
	delete ($adjacency{$_[1]}{$_[0]});
	if (connectionCount($_[0]) == 0) { $adjacency{$_[0]} = (); delete $adjacency{$_[0]}; }
	if (connectionCount($_[1]) == 0) { $adjacency{$_[1]} = (); delete $adjacency{$_[1]}; }
	return ($rval);
}

sub deleteNode {
	my $returnval = 0;
	if (exists $adjacency{$_[0]."i"}) {
		$returnval = 1;
		foreach (keys %{$adjacency{$_[0]."i"}}) {
			deleteLink($_,$_[0]."i");
			}
		}
	if (exists $adjacency{$_[0]."o"}) {
		$returnval = 1;
		foreach (keys %{$adjacency{$_[0]."o"}}) {
			deleteLink($_,$_[0]."o");
			}
		}
	return ($returnval);
}

sub setAdjacencyField {
	$adjacency{$_[0]}{$_[1]}[$_[2]] = $_[3];
	$adjacency{$_[1]}{$_[0]}[$_[2]] = $_[3];
	return();
}

sub setAdjacencyFieldS {
	$scafadjacency{$_[0]}{$_[1]}[$_[2]] = $_[3];
	$scafadjacency{$_[1]}{$_[0]}[$_[2]] = $_[3];
	my @sc0 = nodeOrientation($_[0]);
	my @sc1 = nodeOrientation($_[1]);
	my $nodeor0; my $nodeor1;
	if ($sc0[1] eq "i") { $nodeor0 = $scaffold{$sc0[0]}[1]; } else { $nodeor0 = $scaffold{$sc0[0]}[-1]; }
	if ($sc1[1] eq "i") { $nodeor1 = $scaffold{$sc1[0]}[1]; } else { $nodeor1 = $scaffold{$sc1[0]}[-1]; }
	$adjacency{$nodeor0}{$nodeor1}[$_[2]] = $_[3];
	$adjacency{$nodeor1}{$nodeor0}[$_[2]] = $_[3];
	return();
}

sub newAdjacency {
	push @{$adjacency{$_[0]}{$_[1]}}, @_[2..6];
	push @{$adjacency{$_[1]}{$_[0]}}, @_[2..6];
	return();
}

sub existsAdjacency {
	if (exists $adjacency{$_[0]}) { return (1); }
	return (0);
}

sub existsLink {
	if (exists $adjacency{$_[0]}) {
		if (exists $adjacency{$_[0]}{$_[1]}) {
			return (1);
			}
		}
	return (0);
}

sub connectionCount {
	if (exists $adjacency{$_[0]}) {
		return (scalar keys %{$adjacency{$_[0]}});
		}
	return (0);
}

sub connectionCountS {
	if (exists $scafadjacency{$_[0]}) {
		return (scalar keys %{$scafadjacency{$_[0]}});
		}
	return (0);
}

sub opposite {
	if ($_[0] eq "0i") { return ("0i"); }
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


#### JUNKYARD

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

sub daligner2Coords {
	if ($_[0] eq "") { return (); }
	$_[0] =~ s/,//g; # thousands separators
	$_[0] =~ /(\d+)\s+(\d+)\s+([ncNC])[\s\[]+(\d+)[\.\s]+(\d+)[x\s\]\[]+(\d+)[\.\s]+(\d+)[\]\s\<\:]+(\d+)\sdiffs/;  #    1    9 c   [     0.. 1,876] x [ 9,017..10,825] :   <    398 diffs  ( 18 trace pts)
	# 1seed 2 read 3orientation 4seedstart 5seedend 6readstart 7readend 8diff
	if (lc($3) eq "n") { 
		return ($6+1,$7,$4+1,$5,$7-$6,$5-$4,1-($8/($5-$4)),$2,$1);
		}
	else {
		return ($7,$6+1,$5,$4+1,$6-$7,$5-$4,1-($8/($5-$4)),$2,$1);
		#return ($6+1,$7,$5,$4+1,$7-$6,$5-$4,1-($8/($5-$4)),$2,$1);
		}
}