#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV) != 2) {
	print "Usage: perl change_matrix_gene_id.pl <expression_matrix_file> <gene_id_conversion_file>\n";
	exit;
}

my ($expr_matrix,$gene_ids) = @ARGV;

my %list_hash;

open (my $fh, $gene_ids) || die ("\nERROR: the file ".$gene_ids." could not be found\n");
while (my $line = <$fh>) {
  chomp($line);

  my ($ensembl_id, $new_id) = split("\t",$line);
  
  if ($ensembl_id) {
    $list_hash{$ensembl_id} = $new_id;
  }
}

my $counter = 0;

open (my $fh2, $expr_matrix) || die ("\nERROR: the file ".$expr_matrix." could not be found\n");
while (my $line = <$fh2>) {
	chomp($line);
  
  if ($counter == 0) {
    print "$line\n";
    $counter++;
  }
  
	my @line_array = split("\t",$line);

  my $ref_name = shift(@line_array);

	if ($ref_name) {
	  if ($list_hash{$ref_name}) {
	    print $list_hash{$ref_name}."\t".join("\t",@line_array)."\n";
	  } else {
			print "\n\n\n!!!!!!!!!!!!!! NO ID? !!!!!!!!!!!!!!!!\n\n\n";
	  }
	}

}
