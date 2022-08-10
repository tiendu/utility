use strict;

my $functional_profile_path = shift @ARGV;
my $fasta_single_line_path = shift @ARGV;

my $functional_profile_path_cp = $functional_profile_path;
$functional_profile_path_cp =~ s/.*\///;
my $file_name = $1 if = $functional_profile_path_cp =~ /^(.+)\.([^.]+)$/;

open my $functional_profile, "<:utf8", $functional_profile_path or die;
my %positions_hash;
while (<$functional_profile>) {
    chomp;
    my @fields = split "\t";
    my $id = $fields[0];
    my $start = $fields[1];
    my $end = $fields[2];
    $positions_hash{$id} = [] unless exists $positions_hash{$id};
    push @{$positions_hash{$id}}, [$start, $end];
    last if eof $functional_profile;
};
close $functional_profile;

my %sequences_hash;
my $flag = 0;
my $stored_id;
open my $fasta, "<:utf8", $fasta_single_line_path or die;
while (<$fasta>) {
    chomp;
    if ($_ =~ m/\A>(.+)/) {
        my @fields = split " ";
        my $id = $fields[0];
        $id =~ s/>//;
        $sequences_hash{$id} = '' unless exists $sequences_hash{$id};
        $stored_id = $id;
        $flag = 1;
    } elsif ($flag) {
        $sequences_hash{$stored_id} = $_;
        $flag = 0;
    };
    last if eof $fasta;
};
close $fasta;

open my $orf_fasta, ">:utf8", "orf_${file_name}.fasta";
for my $id (sort keys %sequences_hash) {
    for my $item (@{$positions_hash{$id}}) {
        my ($start, $end) = @$item;
        my $seq = $sequences_hash{$id};
        if ($start < $end) {
            my @orf = split(undef, $seq);
            @orf = @orf[$start-1 .. $end-1];
            print $orf_fasta ">${id}_(${start}:${end})\n";
            print $orf_fasta join("", @orf, "\n");
        } elsif ($start > $end) {
            my @orf = split(undef, reverse $seq);
            @orf = reverse @orf[$end-1 .. $start-1];
            print $orf_fasta ">${id}_(${start}:${end})\n";
            print $orf_fasta join("", @orf, "\n");
        };
    };
};
close $orf_fasta;
