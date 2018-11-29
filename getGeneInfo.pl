#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use REST::Client;
#use LWP::UserAgent;
use Net::Ping;
use File::Basename;
use modules::transcript;
use modules::segment;
use modules::domain;

$Getopt::Std::STANDARD_HELP_VERSION = 1;

############################################################################################################
##	Script to retrieve all info necessary to create a new gene in ushvam2				  ##
##	david baux 09/2017										  ##
##	david.baux@inserm.fr										  ##
############################################################################################################


#######
# the idea is to:
#	-read a file extracted from biomart which gives from gene name => ENST, ENSP, HGNCid, RefSeq NM, UNIPROT
#	-then use a another file got from NCBI (LRG_RefSeqGene.txt) to get NM versions, NG, NP, main transcript
#	-then use togows to get via ucsc api chr, exon, tss, exons positions, strand...http://togows.org/api/ucsc/hg19/refGene/name2=actg1
#	-and UNIPROT for port name, size...
#######

my ($MART, $REFGENE, $HGNC_FILE, $LIFTOVER, $LIFTOVER_CHAIN);

if (-f 'data/mart_export.txt') {
	$MART = "data/mart_export.txt"; #http://www.ensembl.org/biomart/martview/ format: Gene stable ID	Transcript stable ID	Gene name	RefSeq mRNA ID	Protein stable ID	HGNC ID	UniProtKB Gene Name ID
}
else {'die no MART file, you should download a biomart file with ENST, ENSP, HGNCid, RefSeq NM, UNIPROT'}
if (-f 'data/LRG_RefSeqGene.txt') {
	$REFGENE = "data/LRG_RefSeqGene.txt"; #dowloaded from 	ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene
}
if (-f 'data/HGNC_coding.txt') {
	$HGNC_FILE = "data/HGNC_coding.txt"; #dowloaded from http://www.genenames.org/cgi-bin/statistics
}
else {'die no refSeq file, you should download ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene'}
if (-f 'liftover/liftOver_i386') {$LIFTOVER = 'liftover/liftOver_i386'}
else {'die no liftover binary, you should download one for your system at UCSC http://hgdownload.soe.ucsc.edu/admin/exe/'}
if (! -d 'results') {mkdir('results', '0755')}
if (! -d 'tmp') {mkdir('tmp', '0755')}

my $BEDTOOLS = '/usr/local/bin/bedtools';
if (!-f $BEDTOOLS) {undef $BEDTOOLS}

#check togows.org availability
my $p = Net::Ping->new();
if (!defined($p->ping("togows.org", 1))) {die "\n togows.org is not reachable, please check your internet connection\n"}

my (%opts, $filename, $path, $genome, $offset, @transcript, %transcript_hash, %segment_hash, %domain_hash, %uniprot_hash);#, $gene, $segment, $domaine);
getopts('sne:l:g:o:', \%opts);

if ((not exists $opts{'l'}) || ($opts{'l'} !~ /\.txt$/o) || (not exists $opts{'g'})  || $opts{'g'} !~ /hg(19|38)$/o) {
	&HELP_MESSAGE();
	exit
}

if ($opts{'l'}) {($filename, $path) = fileparse($opts{'l'}, qr/\.[^.]*/)}
#if ($opts{'l'} =~ /^(.+)([^\/]+)\.txt$/o) {$list = $2;$path = $1.$list} #get file path and prefix
#elsif ($opts{'l'} =~ /^([^\/]+)\.txt$/o) {$list = $1; $path = $list}
if ($opts{'g'} =~ /hg(19|38)/) {$genome = "hg$1"}

if ($genome eq 'hg19' && -f 'liftover/hg19ToHg38.over.chain.gz') {$LIFTOVER_CHAIN = 'liftover/hg19ToHg38.over.chain.gz'}
else {'die no liftover hg19238 chain, you should download one at UCSC http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'}
if ($genome eq 'hg38' && -f 'liftover/hg19ToHg38.over.chain.gz') {$LIFTOVER_CHAIN = 'liftover/hg38ToHg19.over.chain.gz'}
else {'die no liftover hg38219 chain, you should download one at UCSC http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg38ToHg19.over.chain.gz'}

if ($opts{'o'} && $opts{'o'} =~ /(\d+)/o) {$offset = $1}
else {$offset = 0}
my $gene_counter = 0;
my $email='';#email to put as header in UIPROT REST request
if ($opts{'e'} && $opts{'e'} =~ /^([\w\.]+@\w+\.\w{2,3})$/o) {$email = $1}
#prepare error file
#open E, ">".$list."_error.txt";
#close E;


&populate();
&main();

exit;

sub populate {
	open F, $path.$filename.'.txt' or die $!;
	#my ($i, $j) = (0, 0);
	my @genes;
	print "\nGenome: $genome\n";
	#my $client = REST::Client->new(timeout => 10);
	#my $ua = LWP::UserAgent->new('agent' => 'Mozilla/5.0');
	#$ua->agent('Mozilla/5.0');
	local $| = 1;
	my $iso = 1;#to count specific genes isoforms for protein names which should be unique un ushvam2
	my $current_gene;
	while (my $gene_name = <F>) {
		#if (/(\w+)/o) {push @genes, $1}	
		chomp($gene_name);
		if ($gene_name =~ /(\w+)/o) {
			###Todo check HGNC ok
			my $valid_hgnc = 0;
			open H, "<$HGNC_FILE" or die "Can't find $HGNC_FILE which is supposed to be the HGNC file\n";
			while (<H>) {
				#print $_;
				#print "\n".(split(/\t/, $_))[1]."\n";
				if ($gene_name eq (split(/\t/, $_))[1]) {$valid_hgnc = 1;last;}}
			close H;
			if ($valid_hgnc != 1) {print "\n$gene_name is not a valid HGNC gene coding name. Please check and resubmit\n";&error($gene_name);next;}
			
			if ($iso == 1) {$current_gene = $gene_name}			
			if ($current_gene ne $gene_name) {$iso = 1}
			print "Treating $gene_name...";
			$gene_counter++;
			#my $gene_name = $_;
			my ($i, $j, $enst, $ensp, $nm, $hgnc, $uniprot, $name) = (0, 0, 0, 0, 0, 0, 0, 0);
			open G, "$MART" or die $!; #first file biomart file for NM, ENST, ENSP, UNIPROT, HGNC
			print "Searching transcripts IDs...";
			while (my $line = <G>) {
				chomp($line);
				if ($line !~ /ENSG/o) {#look for columns
					my @header = split(/\t/, $line);
					foreach (@header) {#define headers positions
						#print "$_\n";
						if (/Gene\sname/o){$name = $i}
						elsif (/Transcript\sstable\sID/o){$enst = $i}
						elsif (/RefSeq\smRNA\sID/o){$nm = $i}
						elsif (/Protein\sstable\sID/o){$ensp = $i}
						elsif (/HGNC\sID/o){$hgnc = $i}
						elsif (/UniProtKB/o){$uniprot = $i}
						$i++
					}
					next;
				}
				if ($line =~ /\s$gene_name\s/) {
					my @content = split(/\t/, $line);
					#print "$line\n";
					if ($content[$nm] ne '' && !$transcript_hash{$content[$nm]}) {#unique NM
						my $transcript = transcript->new($content[$nm]);
						$transcript->setENST($content[$enst]);
						$transcript->setENSP($content[$ensp]);
						if ($content[$hgnc] ne '') {
							if ($content[$hgnc] =~ /HGNC:(\d+)/o) {
								$transcript->setHGNC($1);
								my @hgnc_data = split(/\t/, `grep "$content[$hgnc]" $HGNC_FILE`);
								if ($hgnc_data[8]) {
									$hgnc_data[8] =~ s/"//g;
									if ($hgnc_data[8] =~ /([\w-]+)|/o) {$hgnc_data[8] = $1}
									#print $hgnc_data[8];
									if (length $hgnc_data[8] < 21) {$transcript->setSecondName($hgnc_data[8])}
									else {$transcript->setSecondName(substr($hgnc_data[8], 0, 20))}
								}						
								#first occurence of Uniprot ID - is also checked in mart export
								if ($hgnc_data[25] =~ /([\w]{6})$/o) {$transcript->setUniprot($1)}
							}
						}					
						$transcript->setGeneName($content[$name]);
						
						#if ($content[$uniprot] && $content[$uniprot] ne '') {
						if ($content[$uniprot] && length($content[$uniprot]) == 6 || $transcript->getUniprot() ne '') {#swissprot only
							#print "$content[$uniprot]\n";
							if ($transcript->getUniprot() eq '') {$transcript->setUniprot($content[$uniprot])}
							#get info from UNIPROT and create domain objs
							#my $response = $ua->get('http://www.uniprot.org/uniprot/'.$transcript->getUniprot().'.txt');
							
							my (@domains, $data);
							if (! defined($uniprot_hash{$transcript->getUniprot()})) {
								my $uniprot_code = $transcript->getUniprot();
								#my $data = `wget http://www.uniprot.org/uniprot/$uniprot_code.txt`;
								$data = `curl https://www.uniprot.org/uniprot/$uniprot_code.txt`;
								#print "curl results:$data\n";
								##remove ### to get REST UNIPROT call
								###my $uniprot_client = REST::Client->new(timeout => 10);
								###$uniprot_client->addHeader('contact_email'=>$email);
								###my $response = $uniprot_client->GET('https://www.uniprot.org/uniprot/'.$transcript->getUniprot().'.txt');
								#my $response = $ua->get('https://www.uniprot.org/uniprot/'.$transcript->getUniprot().'.txt');
								###if ($uniprot_client->responseCode() == 200) {
									#if ($response->is_success()) {
									#print $response->decoded_content;  # or whatever
									#$data = $response->decoded_content();
									###$data = $uniprot_client->responseContent();
									###$uniprot_hash{$transcript->getUniprot()} = $data;
									
								###}
								###else {
									#print "\nUNIPROT did not respond for ".$transcript->getGeneName()."-".$transcript->getNM()." (UNIPROT ID ".$transcript->getUniprot()."): ".$response->status_line()."\n"
									###print "\nUNIPROT did not respond for ".$transcript->getGeneName()."-".$transcript->getNM()." (UNIPROT ID ".$transcript->getUniprot().", 'https://www.uniprot.org/uniprot/".$transcript->getUniprot().".txt'): ".$uniprot_client->responseCode()."\n"
								###}
									##end remove
							}
							else {$data = $uniprot_hash{$transcript->getUniprot()}}
							my ($short_prot) = (ucfirst(lc($transcript->getGeneName())));
							$short_prot =~ s/'/'\''/og;
							if ($iso > 1) {$short_prot .= " ($iso)"}
							if ($data) {
								my @uniprot = split(/\n/, $data);								
								foreach (@uniprot) {
									#print $_;
									if (/^ID\s+\w+\s+\w+\;\s+(\d+)\sAA\./o) {$transcript->setProtSize($1)}
									elsif (/^DE\s+RecName:\sFull=([\w\s,\/'-]+)[\;\(\{]/o) {
										my $prot_name = $1;
										$prot_name =~ s/'/''/og;
										if ($iso > 1) {$transcript->setProtName("$prot_name ($iso)")}
										else {$transcript->setProtName($prot_name)}
										$transcript->setShortProt($short_prot);
									}
									elsif (/^FT\s+(DOMAIN|MOTIF|TRANSMEM|SIGNAL|TOPO_DOM|REGION|COMPBIAS|REPEAT|COILED)\s+(\d+)\s+(\d+)\s+(.+)\./o) {
										my ($type, $start_aa, $end_aa, $dom_name) = ($1, $2, $3, $4);
										if ($start_aa <= $transcript->getProtSize() || $end_aa <= $transcript->getProtSize()) {
											if ($type eq 'SIGNAL') {$dom_name = lc($type);$dom_name = ucfirst($dom_name)." peptide";}
											elsif ($type eq 'TRANSMEM') {$dom_name = lc($type)."brane";$dom_name = ucfirst($dom_name);}
											elsif ($type eq 'COILED') {$dom_name = 'Coiled Coil'}
											elsif ($type eq 'TOPO_DOM' || $type eq 'REGION') {
												if ($dom_name =~ /(\w+)\.\s\{ECO:\d+\}/o) {$dom_name = $1}
											}
											$dom_name =~ s/ \(Potential\)//og;
											my $domain = domain->new($transcript->getNM(), $gene_name);
											
											$domain->setName($dom_name);
											$domain->setStartAA($start_aa);
											$domain->setEndAA($end_aa);
											push @domains, $domain;
										}
										
									}
								}
								if((!$transcript->getProtName()))  {									
									$transcript->setProtName($short_prot);
									$transcript->setShortProt($short_prot);
								}
								$domain_hash{"$gene_name-$content[$nm]"} = \@domains;
							}
							else {
								$transcript->setProtName($short_prot);
								$transcript->setShortProt($short_prot);
								$transcript->setProtSize('NULL');
							}
							$iso++;
							push @transcript , $transcript;
							print $transcript->getNM();
							$transcript_hash{"$gene_name-$content[$nm]"} = $transcript;
							$j = 1;
						}
						elsif ($opts{'n'} && $j == 0) {
							$transcript->setUniprot('NULL');
							$transcript->setProtName('NULL');
							$transcript->setShortProt('NULL');
							$transcript->setProtSize('NULL');
							push @transcript , $transcript;
							print $transcript->getNM();
							$transcript_hash{"$gene_name-$content[$nm]"} = $transcript;
	;						}
						#}
					}
				}
			}
			close G;
			
			$iso++;
			if ($j == 0) {&error($gene_name)}
			my ($ng, $np, $main, $nm_ver, $actual_nm);
			($i, $j, $ng, $nm, $np, $main, $name,) = (0, 0, 0, 0, 0, 0, 0);
			my $k = 0;
			open G, "$REFGENE" or die $!; #second file refgene file for NM version, NG, NP, main
			print "  Searching refSeq IDs...";
			while (my $line = <G>) {
				chomp($line);
				if ($line =~ /^#/o) {#look for columns
					my @header = split(/\t/, $line);
					foreach (@header) {#define headers positions
						#print "$_\n";
						if (/GeneID$/o){$name = $i}
						elsif (/RSG$/o){$ng = $i}
						elsif (/RNA$/o){$nm = $i}
						elsif (/Protein$/o){$np = $i}
						elsif (/Category$/o){$main = $i}
						$i++
					}
					next;
				}
				if ($line =~ /\s$gene_name\s/) {
					my @content = split(/\t/, $line);
					my $complete_nm = $content[$nm];
					if ($complete_nm =~ /(N[RM]_\d+)\.(\d)/o) {
						#print $complete_nm;exit;
						($actual_nm, $nm_ver) = ($1, $2);
						if ($content[$nm] ne '' && $transcript_hash{"$gene_name-$actual_nm"}) {#si on a un NM d'intérêt
							my $transcript = $transcript_hash{"$gene_name-$actual_nm"};
							$transcript->setNG($content[$ng]);
							$transcript->setNP($content[$np]);
							$transcript->setNMVersion($nm_ver);
							if ($content[$main] =~ /reference/o) {$transcript->setMain('t')}
							else {$transcript->setMain('f')}#beware, there can be several main; must be checked by hand
							$k++;
						}
					}
					else {$j++}
				}
			}
			close G;
			if ($k == 0 && $j != 0) {&error($gene_name)}
			$j = 0;
			#ucsc via togows to get chr, strand, exons
			#http://togows.org/api/ucsc/hg19/refGene/name2=actg1
			print "Searching general infos and defining exons/introns...";
			my $togows_client = REST::Client->new(timeout => 10);
			$togows_client->GET("http://togows.org/api/ucsc/$genome/refGene/name2=$gene_name");
			if ($togows_client->responseCode() == 200) {
				my ($chr, $strand, $txstart, $txend, $cdsstart, $cdsend, $exon_count, $exon_start, $exon_end, $exon_frames);
				$i = 0;
				undef($nm);
				push my @info, split(/\n/, $togows_client->responseContent());
				foreach my $line (@info)	{
					#print "$line--------\n";		
					if ($line =~ /bin/o) {
						my @header = split(/\t/, $line);
						foreach (@header) {#define headers positions
							#print "$_\n";
							if (/name$/o){$nm = $i}
							elsif (/chrom$/o){$chr = $i}
							elsif (/strand$/o){$strand = $i}
							elsif (/txStart$/o){$txstart = $i}
							elsif (/txEnd$/o){$txend = $i}
							elsif (/cdsStart$/o){$cdsstart = $i}
							elsif (/cdsEnd$/o){$cdsend = $i}
							elsif (/exonCount$/o){$exon_count = $i}
							elsif (/exonStarts$/o){$exon_start = $i}
							elsif (/exonEnds$/o){$exon_end = $i}
							elsif (/exonFrames$/o){$exon_frames = $i}
							$i++;
						}
						next;
					}
					else {
						if ($line =~ /\t$gene_name\t/) {					
							my @content = split(/\t/, $line);
							#print "$gene_name-$content[$nm]\n";
							if ($content[$nm] ne '' && $content[$chr] =~ /^chr[\dXYM]{1,2}$/ && $transcript_hash{"$gene_name-$content[$nm]"}) {#si on a un NM d'intérêt
								my $transcript = $transcript_hash{"$gene_name-$content[$nm]"};
								$transcript->setChr($content[$chr]);
								$transcript->setStrand($content[$strand]);
								$transcript->setNbExons($content[$exon_count]);
								#my ($pos1, $pos2) = ($content[$txstart], $content[$cdsstart]);
								my (@exon_start, @exon_end, @exon_frames);
								if ($transcript->getStrand() eq '+') {
									$transcript->setTss(($content[$cdsstart]-$content[$txstart]+1));
									@exon_start = split(/,/, $content[$exon_start]);
									@exon_end = split(/,/, $content[$exon_end]);
									@exon_frames = split(/,/, $content[$exon_frames]);
								}
								else {
									$transcript->setTss(($content[$txend]-$content[$cdsend]+1));
									@exon_start = reverse(split(/,/, $content[$exon_end]));
									@exon_end = reverse(split(/,/, $content[$exon_start]));
									@exon_frames = reverse(split(/,/, $content[$exon_frames]));
								}
								my ($k, $prev) = (0, 0);#$prev is meant to keep in memory the last relevant nuc position to create intronic objects
								my @segments;
								for ($k = 0;$k <= $transcript->getNbExons();$k++) {
									my ($start, $end);									
									my $strand = $transcript->getStrand();
									if ($k != $transcript->getNbExons()) {
										($start, $end) = ($exon_start[$k]+1, $exon_end[$k]);#ucsc is 0-based for start and 1-based for end
										if ($strand eq '-') {($start, $end) = ($exon_start[$k], $exon_end[$k]+1)}
									}
									my $segment = segment->new($content[$nm], $gene_name);
									if ($k == 0) {
										#UTR5-exon1
										$segment->setType('5UTR');
										$segment->setNumber('-1');
										$segment->setName('5UTR');
										if ($strand eq '+') {
											if ($genome eq 'hg19') {
												$segment->setStartG($start-2001);
												$segment->setEndG($start-1);
											}
											elsif ($genome eq 'hg38') {
												$segment->setStartG38($start-2001);
												$segment->setEndG38($start-1);
											}
										}
										else {
											if ($genome eq 'hg19') {
												$segment->setStartG($start+2001);
												$segment->setEndG($start+1);
											}
											elsif ($genome eq 'hg38') {
												$segment->setStartG38($start+2001);
												$segment->setEndG38($start+1);
											}
										}
										$segment->setSize(2000);
										push @segments, $segment;
										
										my $segment_exon1 = segment->new($content[$nm], $gene_name);
										$segment_exon1->setNumber('1');
										$segment_exon1->setType('exon');
										if ($genome eq 'hg19') {
											$segment_exon1->setStartG($start);
											$segment_exon1->setEndG($end);
										}
										elsif ($genome eq 'hg38') {
											$segment_exon1->setStartG38($start);
											$segment_exon1->setEndG38($end);
										}
										
										$segment_exon1->setExonFrame($exon_frames[$k]);
										if ($strand eq '+') {$segment_exon1->setSize($end-$start+1)}
										else {$segment_exon1->setSize($start-$end+1)}
										#if ($genome eq 'hg19') {$prev = $segment_exon1->getEndG()}
										#elsif ($genome eq 'hg38') {$prev = $segment_exon1->getEndG38()}
										#liftover $segment_exon1 if hg19
										if ($genome eq 'hg19') {
											&pyliftover($transcript->getChr(),$segment_exon1->getStartG(), $segment_exon1->getEndG(), $strand, $segment_exon1);
											$prev = $segment_exon1->getEndG();
										}
										elsif ($genome eq 'hg38') {
											&pyliftover($transcript->getChr(),$segment_exon1->getStartG38(), $segment_exon1->getEndG38(), $strand, $segment_exon1);
											$prev = $segment_exon1->getEndG38()
										}
										push @segments, $segment_exon1;
									}
									elsif ($k == $transcript->getNbExons()) {
										#UTR3
										$segment->setType('3UTR');
										$segment->setNumber($transcript->getNbExons()+1);
										$segment->setName('3UTR');
										if ($strand eq '+') {
											if ($genome eq 'hg19') {
												$segment->setStartG($prev+1);
												$segment->setEndG($prev+2001);
											}
											elsif ($genome eq 'hg38') {
												$segment->setStartG38($prev+1);
												$segment->setEndG38($prev+2001);
											}
										}
										else {
											if ($genome eq 'hg19') {
												$segment->setStartG($prev-1);
												$segment->setEndG($prev-2001);
											}
											elsif ($genome eq 'hg38') {
												$segment->setStartG38($prev-1);
												$segment->setEndG38($prev-2001);
											}
										}
										$segment->setSize(2000);
										push @segments, $segment;
									}
									else {
										#exons 2->n / introns 1->n-1
										$segment->setType('exon');
										$segment->setNumber($k+1);
										if ($genome eq 'hg19') {
											$segment->setStartG($start);
											$segment->setEndG($end);
										}
										elsif ($genome eq 'hg38') {
											$segment->setStartG38($start);
											$segment->setEndG38($end);
										}
										$segment->setExonFrame($exon_frames[$k]);
										#print $segment->getExonFrame();
										
										my $segment_intron = segment->new($content[$nm], $gene_name);
										$segment_intron->setNumber($k);
										if ($strand eq '+') {
											$segment->setSize($end-$start+1);
											if ($genome eq 'hg19') {
												$segment_intron->setStartG($prev+1);
												$segment_intron->setEndG($segment->getStartG()-1);
												$segment_intron->setSize($segment->getEndG()-$prev+1);
											}
											elsif ($genome eq 'hg38') {
												$segment_intron->setStartG38($prev+1);
												$segment_intron->setEndG38($segment->getStartG38()-1);
												$segment_intron->setSize($segment->getEndG38()-$prev+1);
											}											
										}
										else {
											$segment->setSize($start-$end+1);
											if ($genome eq 'hg19') {
												$segment_intron->setStartG($prev-1);
												$segment_intron->setEndG($segment->getStartG()+1);
												$segment_intron->setSize($prev-1-$segment->getStartG());
											}
											elsif ($genome eq 'hg38') {
												$segment_intron->setStartG38($prev-1);
												$segment_intron->setEndG38($segment->getStartG38()+1);
												$segment_intron->setSize($prev-1-$segment->getStartG38());
											}
											
										}
										#liftover $segment_intron
										if ($genome eq 'hg19') {
											&pyliftover($transcript->getChr(),$segment_intron->getStartG(), $segment_intron->getEndG(), $strand, $segment_intron);
											$prev = $segment->getEndG();
										}
										elsif ($genome eq 'hg38') {
											&pyliftover($transcript->getChr(),$segment_intron->getStartG38(), $segment_intron->getEndG38(), $strand, $segment_intron);
											$prev = $segment->getEndG38();
										}
										
										push @segments, $segment_intron;
										push @segments, $segment;
									}
									#liftover $segment
									if ($genome eq 'hg19') {&pyliftover($transcript->getChr(),$segment->getStartG(), $segment->getEndG(), $strand, $segment)}
									if ($genome eq 'hg38') {&pyliftover($transcript->getChr(),$segment->getStartG38(), $segment->getEndG38(), $strand, $segment)}
								}
								#print "$segments[0]\n";
								$segment_hash{"$gene_name-$content[$nm]"} = \@segments;
								
								#####put here toBed toSQL for segments
								#open BED, '>>'.$list."_exons_$genome.bed";
								#foreach my $seg (@segments) {
								#	if ($seg->getType() eq 'exon') {
								#		print BED $seg->toBed($transcript->getChr());
								#	}
								#}
								#close BED;
								#####
								#my $segment = segment->new($content[$nm]);
							}
						}
					}
				}
			}
			else {
				print "\n\nBEWARE: togows REST server answered ".$togows_client->responseCode()." to our request for gene $gene_name\ngene will be put in error list, you will need to rerun it\n\n";
				&error($gene_name);
			}
			print "\n";
		}
	}
	close F;
	sleep 1;
}

sub main {
	my $isoform_counter = 0;
	open LOVD, '>results/'.$filename.'_LOVD_domains.txt';
	open BED, '>results/'.$filename."_exons_$genome.bed";
	print BED "track name=$genome.$filename description=\"RefSeq exons with offset $offset of User supplied gene list\"\n";
	open INFO, '>results/'.$filename.'_info.txt';
	open SUMMARY, '>results/'.$filename.'_summary.txt';
	print SUMMARY "Gene-RefSeq\tNG AccNo\tMain Isoform\n";
	if ($opts{'s'}) {open SQL, '>results/'.$filename.'_SQL.sql'}
	#if ($genome eq 'hg19' && $opts{'s'}) {open SQL, '>results/'.$filename.'_SQL.sql'}
	#elsif ($opts{'s'}) {print "\nIgnoring SQL option -s in non hg19 context\n"}
	foreach my $key (sort keys %transcript_hash) {	
	#foreach my $obj (@transcript) {
		#$key =~ /(\w+)-(N[RM]_\d+)/o;
		#my ($gene, $nm) = ($1, $2);
		$isoform_counter++;
		my $obj_transcript = $transcript_hash{$key};
		print INFO "$key\n";
		print SUMMARY "$key.".$obj_transcript->getNMVersion()."\t".$obj_transcript->getNG()."\t".$obj_transcript->getMain()."\n";
		if ($obj_transcript->getChr()){#if false, rest client failed => the gene is in the error file and needs to be reran
			print INFO $obj_transcript->toPrint();
			#if ($genome eq 'hg19' && $opts{'s'}) {print SQL $obj_transcript->toSQL()}
			if ($opts{'s'}) {print SQL $obj_transcript->toSQL()}
			my $segment_list = $segment_hash{$key};
			if ($genome eq 'hg19') {print INFO "#Gene\tNM\tsType\tsNumber\tsName\tsSize\tsStartG\tsEndG\tsStartG38\tsEndG38\tsExonFrame\n"}
			else {print INFO "#Gene\tNM\tsType\tsNumber\tsName\tsSize\tsStartG\tsEndG\tsExonFrame\n"}
			####put here toBed toSQL for segments		
			foreach my $obj_segment (@{$segment_list}) {
				print INFO $obj_segment->toPrint();
				if ($obj_segment->getType() eq 'exon') {print BED $obj_segment->toBed($obj_transcript->getChr(), $obj_transcript->getStrand(), $offset)}
				#if ($obj_segment->getType() ne 'intron') {if ($genome eq 'hg19' && $opts{'s'}) {print SQL $obj_segment->toSQL()}}
				#if ($genome eq 'hg19' && $opts{'s'}) {print SQL $obj_segment->toSQL()}
				if ($opts{'s'}) {print SQL $obj_segment->toSQL()}
			}		
			my $domain_list = $domain_hash{$key};
			print INFO "#Gene\tNM\tdName\tdStart\tdEnd\n";
			print LOVD "#".$obj_transcript->getGeneName()."--".$obj_transcript->getProtName()."\n";
			foreach my $obj_domain (@{$domain_list}) {
				print INFO $obj_domain->toPrint();			
				print LOVD $obj_domain->toLOVD();
				#if ($genome eq 'hg19' && $opts{'s'}) {print SQL $obj_domain->toSQL($obj_transcript->getShortProt())}
				if ($opts{'s'}) {print SQL $obj_domain->toSQL($obj_transcript->getShortProt())}
			}
		}
	}
	close LOVD;
	close BED;
	close INFO;
	close SUMMARY;
	#if ($genome eq 'hg19' && $opts{'s'}) {close SQL}
	if ($opts{'s'}) {close SQL}
	#bedtools to merge intervals
	if ($BEDTOOLS) {
		my $bed = 'results/'.$filename."_exons_$genome";
		system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed; $BEDTOOLS merge -i $bed.sorted.bed -header -c 4,5,6 -o collapse,distinct,distinct > $bed.sorted.merged.bed";
		unlink "$bed.bed", "$bed.sorted.bed";
	}
	print "\n\n$gene_counter genes and $isoform_counter isoforms treated\n\n";
	#foreach (@{$test}) {print "$_\n"}
}

sub HELP_MESSAGE {
	print "\nUsage: perl getGeneInfo.pl -l path/to/annotated/hgnc_gene_list.txt -g genome_version \nSupports --help or --version\n\n
### This script retrieves various information about a list of genes provided as input
### -l txt file, gene list HGNC approved
### -g genome version, hg19/hg38
### -n considers genes without NCBI NG accession number
### -o offset supplementary region to be applied to BED file (numerical, in bp)
### -s generates SQL files for USHVaM2 (https://neuro-2.iurc.montp.inserm.fr/usher), hg19 only
### see README.md for installation instructions
### contact: david.baux\@inserm.fr\n\n"
}

sub VERSION_MESSAGE {
	print "\nVersion 1.0 01/09/2017\n"
}

sub error {
	my $gene = shift;
	open E, ">>results/".$filename."_error.txt";
	print E "$gene\n";
	#print E "No match for $gene (no UNIPROT)\n";
	close E;
}

sub liftover {
	my ($chr, $pos1, $pos2, $strand, $segment) = @_;
	#not optimized as liftover is ran (and files created) each time a segment is achieved - would be more efficient to run it once in the end, but would also be more error prone, and as the script is meant to run only sometimes...
	#print input file
	my $bed = "$chr\t$pos1\t$pos2\n";
	if ($strand eq '-') {$bed = "$chr\t$pos2\t$pos1\n";}
	open(S, '>tmp/input.bed') or die $!;
        print S $bed;
        close S;
	#liftover
	system "$LIFTOVER tmp/input.bed $LIFTOVER_CHAIN tmp/output.bed  tmp/unmapped.bed &>/dev/null";
	open(T, 'tmp/output.bed') or print "Conversion pb with $bed\n";
	while (my $line = <T>) {
		if ($line =~ /chr[0-9XY]{1,2}\t(\d+)\t(\d+)$/o) {
			if ($strand eq '+') {
				if ($genome eq 'hg19') {
					$segment->setStartG38($1);
					$segment->setEndG38($2);
				}
				elsif ($genome eq 'hg38') {
					$segment->setStartG($1);
					$segment->setEndG($2);
				}
			}
			else {
				if ($genome eq 'hg19') {
					$segment->setStartG38($2);
					$segment->setEndG38($1);
				}
				elsif ($genome eq 'hg38') {
					$segment->setStartG($2);
					$segment->setEndG($1);
				}
			}
		}	
	}
	close T;
}

sub pyliftover {
	my ($chr, $pos1, $pos2, $strand, $segment) = @_;
	$pos1 = $pos1-1;
	$pos2 = $pos2-1;
	if ($chr =~ /^([\dXYM]{1,2})$/o) {$chr = "chr$1"}
	#my $ret =  or die "hg38 gene mutalyzer gene only and $!";
	#print STDERR "/usr/local/bin/python $path/liftover$way.py chr$chr $pos";
	#my $way = 'hg38ToHg19.over.chain.gz';
	#if ($genome eq 'hg19') {$way = 'hg19ToHg38.over.chain.gz'}
	my ($s, $e, $chr_tmp);
	($chr_tmp, $s) = split(/,/, `python liftover/liftover.py $LIFTOVER_CHAIN $chr $pos1`);
	$s =~ s/\)//g;	
	$s =~ s/ //g;
	$s =~ s/'//g;
	#print "-$pos1-$s-\n";
	if ($s =~ /^\d+$/o) {$s = $s+1}
	($chr_tmp, $e) = split(/,/, `python liftover/liftover.py $LIFTOVER_CHAIN $chr $pos2`);
	$e =~ s/\)//g;	
	$e =~ s/ //g;
	$e =~ s/'//g;
	#print "-$pos2-$e-\n";
	if ($e =~ /^\d+$/o) {$e = $e+1}
	#if ($strand eq '+') {
	if ($genome eq 'hg19') {
		$segment->setStartG38($s);
		$segment->setEndG38($e);
	}
	elsif ($genome eq 'hg38') {
		$segment->setStartG($s);
		$segment->setEndG($e);
	}
	#}
	#else {
	#	if ($genome eq 'hg19') {
	#		$segment->setStartG38($e);
	#		$segment->setEndG38($s);
	#	}
	#	elsif ($genome eq 'hg38') {
	#		$segment->setStartG($e);
	#		$segment->setEndG($s);
	#	}
	#}

}

