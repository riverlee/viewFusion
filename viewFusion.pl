#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Thu 05 Jul 2012 12:16:52 PM CDT 
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
use Getopt::Long;
use FindBin;
use File::Which;
#View fusion event by circos plot
my $usage=<<USAGE;
##################################################################
#Author: Jiang Li
#email:  riverlee2008\@gmail.com
#URL:    http://github.com/riverlee/viewFusion
#Info:   Use circos to plot fusion event, please modify circos.conf.template to meet your requirement

#Usage:  perl viewFusion.pl -i/--input <inputFusionResult> -f/--format [fusionResultFormat] -s/--species [human] -b/--build [hg19]
   
     -i/--input   fusion result file from fusion detecting program(e.g, FusionHunter)
     -f/--format  the format of the fusion result, use program name as it format, default is FusionHunter
     -s/--species which species's karyotype will be use in the plot(e.g,human,mouse),default is human
     -b/--build   which version(e.g., hg19,mm8,mm9),default is hg19
     -h/--help    print out this help
     -l/--list    list current supported species' karyotypes in karyotype folder,then exists
##################################################################
USAGE
#
my($in,$format,$species,$build,$list,$help)=(undef,"FusionHunter","human",'hg19',undef,undef);
GetOptions("i|input=s"=>\$in,
		   "f|format=s"=>\$format,
		   "s|species=s"=>\$species,
		   "b|build=s"=>\$build,
		   "h|help"=>\$help,
		   "l|list"=>\$list);	   
###########################
#Start parameter checking
if($help){
	print $usage;
	exit(0);
}

if($list){
	_list();
	exit(0);
}

if(!defined($in)){
	print "[Error] input is not provided yet\n\n";
	print $usage;
	exit(1);
}

if(! -e $in){
	print "[Error] input file '$in' is not exists\n\n";
	print $usage;
	exit(1);
}

my $karyotypefile=$FindBin::Bin."/karyotype/karyotype.".lc($species).".".lc($build).".txt";
if(! -e $karyotypefile){
	print "[Error] Species='$species', build='$build' is not supported yet\n\n";
	_list();
	exit(1);
}

# check whether circos in your $PATH
unless(which('circos')){
	print <<CIRCOS;
	
circos in not in your \$PATH, make sure you have circos installed in your PC.
Visit http://circos.ca/software/install to see how to intall it.

CIRCOS
	exit (1);
}

#End parameter checking
###########################

#Main program
my %mapping=('human'=>'hs','mouse'=>'mm');

if($format eq "FusionHunter"){
	_parseFusionHunter($in,$mapping{lc($species)});
	
}else{
	print <<FORMAT;

format='$format' has not been implemented yet.
Please contact Jiang Li <riverlee2008\@gmail.com> by sending the program name
together with its published paper and example output if possible.
 
FORMAT
	exit(1);
}


_copy($karyotypefile);

system("circos -conf circos.conf");

sub _copy{
	my ($karyotypefile) = @_;
	open(IN,$FindBin::Bin."/circos.conf.template") or die $!;
	my $str = join "",<IN>;close IN;
	$str=~s/replacekaryotype/$karyotypefile/;
	open(OUT,">circos.conf") or die $!;
	print OUT $str;
	close OUT;
}

sub _list{

	opendir(DIR,$FindBin::Bin."/karyotype") or die $!;
	my @files = grep { /^karyotype\.\w+\.\w+\.txt/} readdir DIR;
    #print join "\n",@files;
	print "Current avaiable karyotypes of species\n";
	#print '=' x 40;
	#print "\n";
	#print "#Sepcies\tBuild\tFilename\n";
	printf("%-10s%-10s%-30s\n","Species","Build","Filename");
	foreach (@files){
		/^karyotype\.(\w+)\.(\w+)\.txt/;
		printf("%-10s%-10s%-30s\n",$1,$2,$_);
		#print join "\t",($1,$2,$_);
		#print "\n";
	}
	print "\n";
}





=head FusionHunter result description
All these answers are from the author
1)What's the difference between FusionHunter.fusion and FusionHunter.readthrough?
Readthrough is fusion gene linking two adjacent genes on chromosome, while 
FusionHunter.fusion reports fusion genes from distinct genomic regions.

2)How to get the correlative reads for each pair of candidate fusion/readthrough
  genes by program 'allWriteReadsToDir'?
All unmapped reads with their pair mapped onto those region, and all unmapped reads
with partial alignment (say 25bp on either end of the read) onto those region.

3) Taking the first fusion event in the FusionHunter.fusion file (The first line is shown below) as an example, 
   a) what does '++' mean? Does that mean the two genes are both from forward strand, 
      while '+-' means the first gene is from forward strand and the second form reverse strand? 
   b) What does 5 mean in the '(5)'? 
      Does that mean there are 5 paired-end reads spanning this candidate fusion? 
   c) Saying there are 14 Junction Spanning Reads in the first fusion,
      but there is a line with '-------------------', 
      does that mean there are two potential fusion junctions based on the data.
   d) What does 'known x unknown', 'known x known', etc. mean?

a) Your interpretation is right.

b) Yes

c) Due to alternative splicing, fusion gene may have multiple fusion junctions

d) Known means the fusion junction site is on annotated exon boundary,
   while unknown stands for unannotated exon boundary, which is always less reliable.
   
Example 

# Fusion: FusionHunter_temp/R28 chr17:39654604-39662248+chr12:53199081-53209724[++](5)
-> ACB034_0204:6:1201:3040:66271#TGACCA/1       CCCCACCAAAGCCACCTCCAAGGCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGACACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTT        chr17:39661553-39661579 chr12:53207596-53207668  KRT13 x KRT4	known x unknown
-> ACB034_0204:5:1202:12094:42789#TGACCA/1            CAAAGCCACCTCCAAGGCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGGCACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTTGTTCCC  chr17:39661559-39661579 chr12:53207596-53207674  KRT13 x KRT4	known x unknown
-> ACB034_0204:4:2110:12822:77315#TGACCA/1          ACCAAAGCCACCTCCAAGGCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGGCACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTTGTTC    chr17:39661557-39661579 chr12:53207596-53207672  KRT13 x KRT4	known x unknown
----------------
-> ACB034_0204:6:2304:2761:75222#TGACCA/1                 CTTCTCATGGCCAGTGAGGAGGCCGCCATCACAAGCACCAAAGTCAACAAAGCCACCAGCAAAACCCCCACCAAAGCCAC CAGTGCCAAAGCCTCCAGCA                                                                    chr17:39661488-39661567 chr12:53207596-53207615  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:6:1108:3836:25929#TGACCA/1                                                                               ACCCCCCCCCAAAGCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGGCACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTTGTTCCCCCT      chr17:39661550-39661567 chr12:53207596-53207677  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:6:1102:8875:80367#TGACCA/1             GGATCTTCTCATTGCCAGTGAGGAGGCCGCCATCACAAGCACCAAAGTCAACAAAGCCACCAGCAAAACCCCCACCAAAGCCAC CAGTGCCAAAGCCTCC                                                                        chr17:39661484-39661567 chr12:53207596-53207611  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:5:2115:7480:84701#TGACCA/1       GCATGGTGATCTTCTCATTGCCAGTGAGGAGGCCGCCATCACAAGCACCAAAGTCAACAAAGCCACCAGCAAAACCCCCACCAAAGCCAC CAGTGCCAAA                                                                              chr17:39661478-39661567 chr12:53207596-53207605  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:5:1304:16841:61176#TGACCA/1             GATCTTCTCATTGCCAGTGAGGAGGCCGCCATCACAAGCACCAAAGTCAACAAAGCCACCAGCAAAACCCCCACCAAAGCCAC CAGTGCCAAAGCCTCCA                                                                       chr17:39661485-39661567 chr12:53207596-53207612  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:4:2304:6567:39093#TGACCA/1                                                                                   CCCACCAAAGCCAC CAGTGCCGAAGCCTCCAGCACCCCCAAAGCAGACACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTTGTTCCCCCTGAGG  chr17:39661554-39661567 chr12:53207596-53207681  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:4:1311:14621:5178#TGACCA/1                                                                             AAAGCCCCCACCAAACCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGGCACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTTGTTCCCC        chr17:39661548-39661567 chr12:53207596-53207675  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:4:1305:4642:68848#TGACCA/2                                                                                  CCGCACCAAAGCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGACACCTTGTCGTGACCCAGCCACACTCATGGAGATGCTTTTGTTCCCCCTGAG   chr17:39661553-39661567 chr12:53207596-53207680  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:4:1214:15237:68940#TGACCA/1                 TTCTCATTGCCAGTGAGGAGGCCGCCATCACAAGCACCAAAGTCAACAAAGCCACCAGCAAAACCCCCACCAAAGCCAC CAGTGCCAAAGCCTCCAGCAC                                                                   chr17:39661489-39661567 chr12:53207596-53207616  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:4:1203:6862:57302#TGACCA/1                TCTTCTCATTGCCAGTGAGGAGGCCGCCATCACAAGCACCAAAGTCAACAAAGCCACCAGCAAAACCACCACCAAAGCCAC CAGTGCCAAAGCCTCCAGC                                                                     chr17:39661487-39661567 chr12:53207596-53207614  KRT13 x KRT4	unknown x unknown
-> ACB034_0204:4:1109:4615:72042#TGACCA/1                                                      CCAAAGTCAACAAAGCCACCAGCAAAACCCCCACCAAAGCCAC CAGTGCCAAAGCCTCCAGCACCCCCAAAGCAGACACCTTGTCGTGACCCAGCCACAC                               chr17:39661525-39661567 chr12:53207596-53207652  KRT13 x KRT4	unknown x unknown
# Total # of Junction Spanning Reads: 14
=cut

sub _parseFusionHunter{
	my ($in,$chrsuffix) = @_;
	open(IN,$in) or die $!;
	
	#Store fusion genes, to plot this gene names in the circos plot
	my %fusiongenes;	
	#Store fusion boundaries, to plot in circos as link
	my %fusionboundaries;
	
	my $str = "";
	while(<IN>){
		next if (/^$/);
		if(/^# Total/){ #begin to parse
			my @lines = split "\n",$str;
			#$header will look like following
			## Fusion: FusionHunter_temp/R28 chr17:39654604-39662248+chr12:53199081-53209724[++](5)
			## Readthrough: Readthrough_temp/R15 chr17:4013413-4114023+chr17:4119260-4216718[++](4)
			my $header = shift(@lines);
			my ($key1,$key2); #fusion gene, value is chr\tstart\tend
			my ($gene1chr,$gene1start,$gene1end,$gene2chr,$gene2start,$gene2end,
			    $gene1strand,$gene2strand,$supportPairReadsNum);
			   # print $header;
			if($header=~/temp\/R\d+ (\w+):(\d+)\-(\d+)\+(\w+):(\d+)\-(\d+)\[(.)(.)\]\((\d+)\)/){
				($gene1chr,$gene1start,$gene1end,$gene2chr,$gene2start,
				$gene2end,$gene1strand,$gene2strand,$supportPairReadsNum)=($1,$2,$3,$4,$5,$6,$7,$8,$9);
				
				$gene1chr=_dealwithchr($gene1chr,$chrsuffix);
				$gene2chr=_dealwithchr($gene2chr,$chrsuffix);
				$key1=join "\t",($gene1chr,$gene1start,$gene1end);
				$key2=join "\t",($gene2chr,$gene2start,$gene2end);			
				#print join "\t",($1,$2,$3,$4,$5,$6,$7,$8,$9,"\n");
			}else{
				#Maybe the regular expression is not powful
				print "[Warning]:Not a standard FunsionHunter result\n";
				$str="";
				next;
			}
			
			#Parse boundaries supported by Junction spanning reads
			foreach my $l (@lines){
				chomp($l);
				next if ($l!~/^\->/);
				my(undef,$read,$r1,$r2,$chrpos1,$chrpos2,$gene1,undef,$gene2,$isknow1,undef,$isknow2)=split /\s+/,$l;
				#print join "\t",($chrpos1,$chrpos2,$gene1,$gene2,$isknow1,$isknow2,"\n");
				$fusiongenes{$key1}=$gene1;
				$fusiongenes{$key2}=$gene2;
				
				#Determine fusion boundaries
				my($chromosome1,$start1,$end1) = split /:|\-/,$chrpos1;
				my($chromosome2,$start2,$end2) = split /:|\-/,$chrpos2;
				$chromosome1=_dealwithchr($chromosome1,$chrsuffix);
				$chromosome2=_dealwithchr($chromosome2,$chrsuffix);
				my ($boundary1,$boundary2)=($start1,$start2);
				if($gene1strand eq "+"){
					$boundary1=$end1;
				}
				if($gene2strand eq "-"){
					$boundary2=$end2;
				}
				
				my $boundary1end=$boundary1+1;
				my $boundary2end=$boundary2+1;
				my $boundarykey = $chromosome1."\t".$boundary1."\t".$boundary1end.":".
								  $chromosome2."\t".$boundary2."\t".$boundary2end;
				$fusionboundaries{$boundarykey} ++;
			}
			$str="";
		}else{
			$str.=$_;
		}
	}
	close(IN);
	#Write out fusiongenes
	open(OUT,">fusion-gene.txt") or die $!;
	foreach my $key(sort keys %fusiongenes){
		print OUT join "\t",($key,$fusiongenes{$key});
		print OUT "\n";
	}
	close OUT;
	
	#Write out fusion boundaries
	open(OUT,">fusion-link.txt") or die $!;
	my $count=1;
	foreach my $boundary (sort keys %fusionboundaries){
		my($s1,$s2)=split ":",$boundary;
		my $key="fusion".$count++;
		print OUT join "\t" ,($key,$s1,"color=".
							_getColor($fusionboundaries{$boundary}).
							",spanningReads=".$fusionboundaries{$boundary});
		print OUT "\n";
		print OUT join "\t" ,($key,$s2,"color=".
							_getColor($fusionboundaries{$boundary}).
							",spanningReads=".$fusionboundaries{$boundary});
		print OUT "\n";
	}
}

sub _dealwithchr{
	my($str,$chrsuffix) = @_;
	$str=uc($str);
	$str=~s/CHR//g;
	return($chrsuffix.$str);
}


#According to value, generate color name
#Need to impove later
sub _getColor{
	my ($n)=@_;
	if($n<5){
		return('blue');
	}elsif($n<10){
		return('orange');
	}elsif($n<15){
		return('green');
	}elsif($n<20){
		return('blue');
	}else{
		return('purple');
	}
}

