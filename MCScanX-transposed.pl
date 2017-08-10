use Getopt::Std;
%options=();
getopts("i:t:o:c:d:k:g:s:e:m:w:a:x:", \%options);
if(!exists $options{i} || !exists $options{t} || !exists $options{o} || !exists $options{c})
{
print "Usage: perl MCScanX-transposed.pl -i data_directory -t target_species -c outgroup_species(comma_delimited) -o output_directory\n";
print "#####################\nOptional:\n";
print "-x number_of_different_epoches (if specified, outgroup species must be provided in the order of divergence from the target species(most recent first), default: 1, only consider the transposed duplications that occurred after the divergence between target species and all outgroups )\n";
print "-a 1 or 0(are segmental duplicates ancestral loci or not? default: 1, yes)\n";
print "-d number_of_genes(maximum distance to call proximal, default: 10)\n";
print "#####################\nThe following are optional MCScanX parameters:\n";
print "-k match_score(cutoff score of collinear blocks for MCScanX, default: 50)\n";
print "-g gap_penalty(gap penalty for MCScanX, default: -1)\n";
print "-s match_size(number of genes required to call a collinear block for MCScanX, default: 5)\n";
print "-e e_value(alignment significance for MCScanX, default: 1e-05)\n";
print "-m max_gaps(maximum gaps allowed for MCScanX, default: 25)\n";
print "-w overlap_window(maximum distance in terms of gene number, to collapse BLAST matches for MCScanX, default: 5)\n";
exit;
}
$options{i}=~s/\/$//;
$options{o}=~s/\/$//;
$good=1;
$t_gff="$options{i}\/$options{t}\.gff";
$t_bla="$options{i}\/$options{t}\.blast";
unless(-e $t_gff)
{
$good=0;
print "Error: Cannot find the gff file for the target species at $t_gff\n";
}
unless(-e $t_bla)
{
$good=0;
print "Error: Cannot find the blast file for the target species at $t_bla\n";
}
@ogrp=split("\,",$options{c});
############################
$n_epoch=1;
if(exists $options{x})
{
if($options{x}>=1 && $options{x}<=$#ogrp+1)
{
$n_epoch=$options{x};
}
else
{
print "Eorror: 'x' must be >=1 and <= # of outgroup species\n";
}
}
###########################
for($i=0;$i<=$#ogrp;$i++)
{
$o_gff="$options{i}\/$options{t}_$ogrp[$i]\.gff";
$o_bla="$options{i}\/$options{t}_$ogrp[$i]\.blast";
unless(-e $o_gff)
{
$good=0;
print "Error: Cannot find the gff file between the target species and the outgroup species \"$ogrp[$i]\" at $o_gff\n";
}
unless(-e $o_bla)
{
$good=0;
print "Error: Cannot find the blast file between the target species and the outgroup species \"$ogrp[$i]\" at $o_bla\n";
}
}
if($good==0)
{
print "##############################################################\n";
print "!!!The execution is terminated due to incorrect input files!!!\n";
print "##############################################################\n";
exit;
}
unless(-d $options{o})
{
mkdir $options{o} or die "Error: Cannot create the output directory\n";
}
@para=("k","g","s","e","m","w");
$mcscanx_para="";
for($i=1;$i<6;$i++)
{
if(exists $options{$para[$i]})
{
$mcscanx_para=$mcscanx_para." -$para[$i] $options{$para[$i]}";
}
}
########################################################################
$in="$options{i}\/$options{t}";
$out="$options{o}\/$options{t}";
system("./MCScanX $in $out $mcscanx_para");
for($i=0;$i<=$#ogrp;$i++)
{
$in="$options{i}\/$options{t}_$ogrp[$i]";
$out="$options{o}\/$options{t}_$ogrp[$i]";
system("./MCScanX $in $out $mcscanx_para");
system("rm $out.gff.sorted");
}
########################################################################
$sorted_gff="$options{o}\/$options{t}.gff.sorted";
open(input,$sorted_gff);
%gid=();
%gch=();
%glc=();
%anc=();
%gmd=();
$i=0;
while($line=<input>)
{
chomp($line);
@a=split("\t",$line);
$gid{$a[0]}=$i;
$gch{$a[0]}=$a[1];
$glc{$a[0]}=$a[2];
$gmd{$a[0]}=0;
$anc{$a[0]}=0;
$i++;
}
close(input);
$cut=10;
if(exists $options{d})
{
$cut=$options{d};
}
#############read blastp of target genome#########################
%bla1=();
%bla2=();
%blae=();
open(input,$t_bla);
while($line=<input>)
{
@a=split("\t",$line);
if($a[0] ne $a[1] && exists $gid{$a[0]} && exists $gid{$a[1]})
{
$key="$a[0]\t$a[1]";
if(!exists $bla2{$key})
{
$bla2{$key}=$a[2];
}
else
{
if($a[2]>$bla2{$key})
{
$bla2{$key}=$a[2];
}
}
$blae{$key}=$a[10];
if($gid{$a[0]} gt $gid{$a[1]})
{
$key="$a[1]\t$a[0]";
}
$bla1{$key}=0;
}
}
#######################segmental#############################################
$in="$options{o}\/$options{t}.collinearity";
$out1=">$options{o}\/$options{t}.segmental.pairs";
$out2=">$options{o}\/$options{t}.segmental.genes";
$incl=1;
if(exists $options{a})
{
$incl=$options{a};
}
%h=();
open(input,$in);
open(output1,$out1);
open(output2,$out2);
print output1 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output2 "Duplicate\tLocation\n";
while($line=<input>)
{
chomp($line);
if($line!~/\#/ && $line ne "")
{
@a=split("\t",$line);
$t_key="$a[1]\t$a[2]";
if(exists $blae{$t_key})
{
$temp=$blae{$t_key};
}
else
{
$t_key="$a[2]\t$a[1]";
$temp=$blae{$t_key};
}
print output1 "$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$a[2]\t$gch{$a[2]}\:$glc{$a[2]}\t$temp\n";
$h{$a[1]}=0;
$h{$a[2]}=0;
if($incl==1)
{
$anc{$a[1]}=1;
$anc{$a[2]}=1;
}
$gmd{$a[1]}=1;
$gmd{$a[2]}=1;
}
}
foreach $key (sort(keys %h))
{
print output2 "$key\t$gch{$key}\:$glc{$key}\n";
}
####################tandem&proximal#####################
$out1=">$options{o}\/$options{t}.tandem.pairs";
$out2=">$options{o}\/$options{t}.tandem.genes";
$out3=">$options{o}\/$options{t}.proximal.pairs";
$out4=">$options{o}\/$options{t}.proximal.genes";
open(output1,$out1);
open(output2,$out2);
open(output3,$out3);
open(output4,$out4);
print output1 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output2 "Duplicate\tLocation\n";
print output3 "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\n";
print output4 "Duplicate\tLocation\n";
%tan_p=();
%tan_d=();
foreach $key (keys %bla1)
{
@a=split("\t",$key);
$dis= abs($gid{$a[1]}-$gid{$a[0]});
if($gch{$a[0]} eq $gch{$a[1]} && $dis<$cut)
{
for($i=0;$i<2;$i++)
{
if(! exists $tan_p{$a[$i]})
{
$tan_p{$a[$i]}=$a[1-$i];
$tan_d{$a[$i]}=$dis;
}
else
{
if($dis<$tan_d{$a[$i]})
{
$tan_p{$a[$i]}=$a[1-$i];
$tan_d{$a[$i]}=$dis;
}
elsif($dis==$tan_d{$a[$i]})
{
if($gid{$a[1-$i]}>$gid{$a[$i]})
{
$tan_p{$a[$i]}=$a[1-$i];
}
}
else
{
}
}
}
}
}
%tan_pairs=();
foreach $key (keys %tan_p)
{
$newkey=$key."\t".$tan_p{$key};
if($gid{$key} gt $gid{$tan_p{$key}})
{
$newkey=$tan_p{$key}."\t".$key;
}
$tan_pairs{$newkey}=$tan_d{$key};
}
%h=();
%k=();
foreach $key (sort(keys %tan_pairs))
{
@a=split("\t",$key);
if($tan_pairs{$key}==1)
{
$t_key="$a[0]\t$a[1]";
if(exists $blae{$t_key})
{
$temp=$blae{$t_key};
}
else
{
$t_key="$a[1]\t$a[0]";
$temp=$blae{$t_key};
}
print output1 "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$temp\n";
for($i=0;$i<2;$i++)
{
if($gmd{$a[$i]}==0)
{
$h{$a[$i]}=2;
$gmd{$a[$i]}=2;
}
}
}
else
{
$t_key="$a[0]\t$a[1]";
if(exists $blae{$t_key})
{
$temp=$blae{$t_key};
}
else
{
$t_key="$a[1]\t$a[0]";
$temp=$blae{$t_key};
}
print output3 "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$temp\n";
for($i=0;$i<2;$i++)
{
if($gmd{$a[$i]}==0 || $gmd{$a[$i]}==3)
{
$k{$a[$i]}=3;
$gmd{$a[$i]}=3;
}
}
}
}
foreach $key (sort(keys %h))
{
print output2 "$key\t$gch{$key}\:$glc{$key}\n";
}
foreach $key (sort(keys %k))
{
print output4 "$key\t$gch{$key}\:$glc{$key}\n";
}
####################transposed&dispersed################
%epoch_pair=();
for($k=$n_epoch;$k>=1;$k--)
{
for($i=$#ogrp;$i>=$k-1;$i--)
{
$in="$options{o}\/$options{t}_$ogrp[$i].collinearity";
open(input,$in);
while($line=<input>)
{
chomp($line);
if($line!~/\#/ && $line ne "")
{
@a=split("\t",$line);
for($j=1;$j<3;$j++)
{
if(exists $anc{$a[$j]})
{
$anc{$a[$j]}=1;
}
}
}
}
}
%tran_genes=();
%tran_ident=();
foreach $key (keys %bla2)
{
$good=0;
@a=split("\t",$key);
if($gmd{$a[0]}==0)
{
if($gch{$a[0]} ne $gch{$a[1]})
{
$good=1;
}
else
{
if(abs($gid{$a[0]}-$gid{$a[1]})>=$cut)
{
$good=1;
}
}
}
if($good==1)
{
if($anc{$a[0]}==0 && $anc{$a[1]}==1)
{
if(!exists $tran_genes{$a[0]})
{
$tran_genes{$a[0]}=$a[1];
$tran_ident{$a[0]}=$bla2{$key};
}
else
{
if($a[2]>$tran_ident{$a[0]})
{
$tran_genes{$a[0]}=$a[1];
$tran_ident{$a[0]}=$bla2{$key};
}
}
}
}
}
if($n_epoch>1)
{
$out1=">$options{o}\/$options{t}.transposed_after_$ogrp[$k-1].pairs";
$out2=">$options{o}\/$options{t}.transposed_after_$ogrp[$k-1].genes";
}
else
{
$out1=">$options{o}\/$options{t}.transposed.pairs";
$out2=">$options{o}\/$options{t}.transposed.genes";
}
open(output1,$out1);
open(output2,$out2);
print output1 "Transposed\tLocation\tParental\tLocation\tE-value\n";
print output2 "Duplicate\tLocation\n";
foreach $key (sort(keys %tran_genes))
{
print output2 "$key\t$gch{$key}\:$glc{$key}\n";
$temp="$key\t$tran_genes{$key}";
print output1 "$key\t$gch{$key}\:$glc{$key}\t$tran_genes{$key}\t$gch{$tran_genes{$key}}\:$glc{$tran_genes{$key}}\t$blae{$temp}\n";
$epoch_pair{$temp}=$k;
}
}
if($n_epoch>1)
{
for($k=$n_epoch;$k>1;$k--)
{
$tp1="out_epoch_pair$k";
$tp2="out_epoch_gene$k";
open($tp1,">$options{o}\/$options{t}.transposed_between_$ogrp[$k-2]_$ogrp[$k-1].pairs");
open($tp2,">$options{o}\/$options{t}.transposed_between_$ogrp[$k-2]_$ogrp[$k-1].genes");
print {$tp1} "Transposed\tLocation\tParental\tLocation\tE-value\n";
print {$tp2} "Duplicate\tLocation\n";
}
open(out_all_pair,">$options{o}\/$options{t}.transposed.pairs");
open(out_all_gene,">$options{o}\/$options{t}.transposed.genes");
print out_all_pair "Transposed\tLocation\tParental\tLocation\tE-value\n";
print out_all_gene "Duplicate\tLocation\n";
%ugene=();
foreach $key (sort(keys %epoch_pair))
{
@a=split("\t",$key);
print out_all_pair "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$blae{$key}\n";
#print out_all_gene "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\n";
$ugene{$a[0]}=1;
if($epoch_pair{$key}>1)
{
$hd1="out_epoch_pair$epoch_pair{$key}";
$hd2="out_epoch_gene$epoch_pair{$key}";
print {$hd1} "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\t$a[1]\t$gch{$a[1]}\:$glc{$a[1]}\t$blae{$key}\n";
print {$hd2} "$a[0]\t$gch{$a[0]}\:$glc{$a[0]}\n";
}
}
foreach $key (sort(keys %ugene))
{
print out_all_gene "$key\t$gch{$key}\:$glc{$key}\n";
}
}
