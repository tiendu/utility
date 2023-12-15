use strict;

my @files = @ARGV;

my %records;
my $tally = 0;
for my $file (@files) {
    open my $file, "<:utf8", $file or die;
    while (<$file>) {
        chomp;
        my ($column_1, $column_2) = split "\t";
        $records{$column_1}[$tally] = $column_2;
    };
    $tally++;
};

open my $merged, ">:utf8", "merged.tsv";
# open my $shared, ">:utf8", "shared.tsv";
# open my $unshared, ">:utf8", "unshared.tsv";
for my $key (sort keys %records) {
    my @values = @{$records{$key}};
    print $merged join("\t", $key, map { $_ // 0 } @values[0 .. $tally - 1], "\n");
#     unless (grep {! defined($_) } @values) {
#         print $shared join("\t", $key, @values, "\n");
#     };
#     if (grep {! defined($_) } @values) {
#         print $unshared join("\t", $key, map { $_ // 0 } @values[0 .. $tally - 1], "\n");
#     };
};
close $merged;
# close $shared;
# close $unshared;
