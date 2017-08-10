use Getopt::Std;
%options=();
getopts("i:d:o:", \%options);
if(! exists $options{i}||! exists $options{d})
{
print "Usage:\nperl detect_tra_epochs_for_a_family.pl -i gene_family_file -d output_directory_of_MCScanX-transposed\/target_species\n";
print "optional: -o output_file\n";
exit;
}
$b=`ls $options{d}.transposed_*.pairs`;
chomp($b);
if($b eq "")
{
print "Please check the path for the MCScanX-transposed output\n";
exit;
}
@f=split("\n",$b);
for($i=0;$i<=$#f;$i++)
{
$pos=index($f[$i],"\.transposed_");
$tmp=substr($f[$i],$pos+1);
$tmp=~s/\.pairs//;
$mode[$i]=$tmp;
}
%dup=();
%gene=();
%pos=();
%eva=();
for($i=0;$i<=$#f;$i++)
{
$tag=$mode[$i];
open($tag,$f[$i]);
$line=<$tag>;
while($line=<$tag>)
{
chomp($line);
@a=split("\t",$line);
$gene{$a[0]}=1;
$gene{$a[2]}=1;
$pos{$a[0]}=$a[1];
$pos{$a[2]}=$a[3];
$temp="$a[0]\t$a[2]";
$dup{$temp}=$mode[$i];
$eva{$temp}=$a[4];
}
}
#print scalar keys %dup,"\n";
open(input1,$options{i}) or die "Cannot open the gene_family_file!\n";
$write_file=0;
if(exists $options{o})
{
open(output,">$options{o}") or die "Cannot open output_file!\n";
$write_file=1;
}
############################################
$j=0;
while($line=<input1>)
{
chomp($line);
@a="";
if($line ne "")
{
@a=split(" ",$line);
for($i=0;$i<=$#a;$i++)
{
$m[$j]=$a[$i];
$j++;
}
}
}
$num=0;
for($k=0;$k<$j;$k++)
{
if(exists $gene{$m[$k]})
{
$g[$num]=$m[$k];
$num++;
}
}
print "Transposed\tLocation\tParental\tLocation\tE-value\tEpoch\n";
print output "Transposed\tLocation\tParental\tLocation\tE-value\tEpoch\n";
for($i=0;$i<$num-1;$i++)
{
for($j=$i+1;$j<$num;$j++)
{
$gid1=$g[$i]."\t".$g[$j];
$gid2=$g[$j]."\t".$g[$i];
if(exists $dup{$gid1})
{
print "$g[$i]\t$pos{$g[$i]}\t$g[$j]\t$pos{$g[$j]}\t$eva{$gid1}\t$dup{$gid1}\n";
if($write_file==1)
{
print output "$g[$i]\t$pos{$g[$i]}\t$g[$j]\t$pos{$g[$j]}\t$eva{$gid1}\t$dup{$gid1}\n";
}
}
if(exists $dup{$gid2})
{
print "$g[$j]\t$pos{$g[$j]}\t$g[$i]\t$pos{$g[$i]}\t$eva{$gid2}\t$dup{$gid2}\n";
if($write_file==1)
{
print output "$g[$j]\t$pos{$g[$j]}\t$g[$i]\t$pos{$g[$i]}\t$eva{$gid2}\t$dup{$gid2}\n";
}
}
}
}
