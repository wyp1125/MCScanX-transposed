use Getopt::Std;
%options=();
getopts("i:d:o:", \%options);
if(! exists $options{i}||! exists $options{d})
{
print "Usage:\nperl detect_dup_modes_for_a_gene.pl -i gene_name -d output_directory_of_MCScanX-transposed\/target_species\n";
print "optional: -o output_file\n";
exit;
}
$write_file=0;
if(exists $options{o})
{
open(output,">$options{o}") or die "Cannot open output_file!\n";
$write_file=1;
}
print "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\tMode\n";
if($write_file==1)
{
print output "Duplicate 1\tLocation\tDuplicate 2\tLocation\tE-value\tMode\n";
}
@mode=("segmental","tandem","proximal","transposed");
for($i=0;$i<4;$i++)
{
$dir=$options{d}."\.".$mode[$i]."\.pairs";
$tag="mode$i";
open($tag,$dir) or die "Cannot open $mode[$i]_duplication_file!\n";
$line=<$tag>;
while($line=<$tag>)
{
chomp($line);
if($line=~/$options{i}/)
{
print "$line\t$mode[$i]\n";
if($write_file==1)
{
print output "$line\t$mode[$i]\n";
}
}
}
}
######################################################
$b=`ls $options{d}.transposed_*.pairs`;
chomp($b);
if($b eq "")
{
print "Transposed duplciations are not classified into different epochs\n";
exit;
}
else
{
print "#################epochs for transposed dupliclations######\n";
if($write_file==1)
{
print output "#################epochs for transposed dupliclations######\n";
}
}
print "Transposed\tLocation\tParental\tLocation\tE-value\tEpoch\n";
if($write_file==1)
{
print output "Transposed\tLocation\tParental\tLocation\tE-value\tEpoch\n";
}
@f=split("\n",$b);
for($i=0;$i<=$#f;$i++)
{
$pos=index($f[$i],"\.transposed_");
$tmp=substr($f[$i],$pos+1);
$tmp=~s/\.pairs//;
$mode[$i]=$tmp;
}
for($i=0;$i<=$#f;$i++)
{
$tag=$mode[$i];
open($tag,$f[$i]);
$line=<$tag>;
while($line=<$tag>)
{
chomp($line);
if($line=~/$options{i}/)
{
print "$line\t$mode[$i]\n";
if($write_file==1)
{
print output "$line\t$mode[$i]\n";
}
}
}
}
