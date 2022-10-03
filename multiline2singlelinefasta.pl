use strict;
use warnings;

my $fasta_path = shift @ARGV;
my $file_path = $fasta_path;
$file_path =~ s/.*\///;
my ($file_name, $file_extension) = $file_path =~ /^(.+)\.([^.]+)$/;

open my $fasta, "<:utf8", $fasta_path or die;
open my $fasta_single_line, ">:utf8", "singleline_${file_name}.${file_extension}";
while (<$fasta>) {
    chomp unless m/\A>(.+)/;
    if ($. > 1 && m/\A>(.+)/) {
        $_ =~ s/\A/\n/;
    };
    print $fasta_single_line "$_";
    print $fasta_single_line "\n" if eof;
};
close $fasta;
close $fasta_single_line;
