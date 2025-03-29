#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV != 6) {
    die "Input: @ARGV;  Usage: $0 file_A file_B index_0 index_1 index_2 index_3\n";
}

my ($file_A, $file_B, $index_0, $index_1, $index_2, $index_3) = @ARGV;

open(my $IN1, $file_A) or die "Cannot open A file: $file_B !";
open(my $IN2, $file_B) or die "Cannot open B file: $file_B !";

my %hash;

# Process A file and store information in the hash
while (<$IN1>) {
    chomp;
    my @A = split /\s+/;
    my $key = "$A[0],$A[1],$A[2],$A[3]";
    $hash{$key} = $_; # Store the entire line in the hash
}

# Process B file and print matching rows from both files
while (<$IN2>) {
    chomp;
    my @B = split /\s+/;
    my $key = "$B[$index_0],$B[$index_1],$B[$index_2],$B[$index_3]";
    if (exists $hash{$key}) {
	    #  print "$hash{$key}\t$_\n"; # Print both A and B lines
	print "$_\n"; # Print only A lines
    }
}

close($IN1);
close($IN2);
