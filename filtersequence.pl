use strict;
use warnings;

my $file = $ARGV[0];
my $opt = $ARGV[1];
my @kws = @ARGV[2 .. $#ARGV];

open my $input, "<:utf8", $file or die;
my (%id_sequence, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(\S+)/;
    } else {
        $id_sequence{$id} = uc $_;
    };
    last if eof $input;
};
close $input;

my (@filt_id);
if ($opt eq "in") {
    @filt_id = grep {my $id = $_; grep {$id =~ /$_/} @kws;} keys %id_sequence;
} elsif ($opt eq "out") {
    @filt_id = grep {my $id = $_; !grep {$id =~ /$_/} @kws;} keys %id_sequence;
};

$file =~ s/.*\///; 
my ($name, $extension) = $file =~ /^(.+)\.([^.]+)$/;
open my $output, ">:utf8", "filtered_${name}.${extension}" or die;
foreach (@filt_id) {
    print $output join("\n", ">${_}", $id_sequence{$_} . "\n");
};

