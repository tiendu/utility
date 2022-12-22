use strict;
use warnings;

my $info_path = shift @ARGV;
my $singlelined_fasta_path = shift @ARGV;

my $info = $info_path =~ s/.*\///r;
my $name = $1 if $info =~ /^(.+)\.([^.]+)$/;
my $fasta = $singlelined_fasta_path =~ s/.*\///r;
my $extension = $1 if $fasta =~ /\.([^.]+)$/;

open my $input, "<:utf8", $info_path or die;
my %info_hash;
while (<$input>) {
    chomp;
    my @fields = split "\t";
    my $id = $fields[0];
    my ($start, $end, $annotation) = ($fields[1], $fields[2], $fields[3]);
    $info_hash{$id} = [] unless exists $info_hash{$id};
    push @{$info_hash{$id}}, [$start, $end, $annotation];
    last if eof $input;
};
close $input;

open my $fasta_file, "<:utf8", $singlelined_fasta_path or die;
my (%id_sequence_hash, $id);
while (<$fasta_file>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(.+)/;
    } else {
        $id_sequence_hash{$id} = uc $_;
    };
    last if eof $fasta_file;
};
close $fasta_file;

open my $orf_fasta, ">:utf8", "orf_${name}.${extension}" or die;
for my $id (sort keys %id_sequence_hash) {
    for my $item (@{$info_hash{$id}}) {
        my ($start, $end, $annotation) = @$item;
        my $seq = $id_sequence_hash{$id}; 
        my $label = ">${id}_${start}_${end}_${annotation}";
        $label =~ s/\s//;
        if ($start < $end) {
            my @orf = split("", $seq);
            @orf = @orf[$start-1 .. $end-1];
            print $orf_fasta "${label}\n";
            print $orf_fasta join("", @orf, "\n");
        } elsif ($start > $end) {
            my @orf = split("", reverse $seq);
            @orf = reverse @orf[$end-1 .. $start-1];
            print $orf_fasta "${label}\n";
            print $orf_fasta join("", @orf, "\n");
        };
    };
};
close $orf_fasta;
