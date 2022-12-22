use strict;
use warnings;
use Getopt::Long;

GetOptions(
    "input|i=s" => \my $file,
    "mode|m=s" => \my $mode,
    "keywords|k=s{,}" => \my @keywords,
);

open my $input, "<:utf8", $file or die;
my (%id_sequence, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(.+)/;
    } else {
        $id_sequence{$id} = uc $_;
    };
    last if eof $input;
};
close $input;

my (@filtered_id);
if ($mode eq "in") {
    @filtered_id = grep {my $id = $_; grep {$id =~ /$_/} @keywords;} keys %id_sequence;
} elsif ($mode eq "out") {
    @filtered_id = grep {my $id = $_; !grep {$id =~ /$_/} @keywords;} keys %id_sequence;
};

$file =~ s/.*\///; 
my ($name, $extension) = $file =~ /^(.+)\.([^.]+)$/;
open my $output, ">:utf8", "filtered_${name}.${extension}" or die;
foreach (@filtered_id) {
    print $output join("\n", ">${_}", $id_sequence{$_} . "\n");
};
