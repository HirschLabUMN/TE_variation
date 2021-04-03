#! /usr/bin/perl -w

use strict;

open (my $freq_fh, '<', "/Users/cnhirsch/Documents/Research/Nathan/Christine_TE_methods/Analysis/LD_analysis_revisions/File_S1_refB73_PAV.txt") || die;
open (my $ld_fh, '<', "/Users/cnhirsch/Documents/Research/Nathan/Christine_TE_methods/Analysis/LD_analysis_revisions/R2/plink_results_SNPs-highest-LD-TE_R2.ld") || die;
open (my $out_fh, '>', "/Users/cnhirsch/Documents/Research/Nathan/Christine_TE_methods/Analysis/LD_analysis_revisions/R2/TE_SNP_LD_with_freq.txt") || die;

my %freq;
<$freq_fh>;
while (my $line = <$freq_fh>) {
  chomp $line;
  my @line = split ("\t", $line);
  my $name = $line[0];
  my $frequency = $line[513];
  $freq{$name} = $frequency;
}


my $line1 = <$ld_fh>;
chomp $line1;
print $out_fh "$line1\tprop_prsent\n";

while (my $line = <$ld_fh>) {
  chomp $line;
  my @line = split ("\t", $line);

  if ($line[2] !~ /^S/) {
    my $TE = $line[2];
    print $out_fh "$line\t$freq{$TE}\n";
  }

  elsif ($line[5] !~ /^S/) {
    my $TE = $line[6];
    print $out_fh "$line\t$freq{$TE}\n";
  }

  else {
    print "$line\n";
  }

}

close $freq_fh;
close $ld_fh;
close $out_fh;
