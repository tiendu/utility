use strict;

my $file_path = shift @ARGV;
my $file_path_cp = $file_path;
$file_path_cp =~ s/.*\///;
my ($file_name, $file_extension) = $file_path_cp =~ /^(.+)\.([^.]+)$/;

open my $file, "<:utf8", $file_path or die;
my @sequence_lengths;
my $flag = 0;
while (<$file>) {
    chomp;
    if ($_ =~ m/\A>(.+)/) {
        $flag = 1;
    } elsif ($flag) {
        push @sequence_lengths, length($_);
        $flag = 0;
    };
    last if eof $file;
};
close $file;

open my $result, ">:utf8", "${file_name}_Nx.tsv";
for (my $i = 0; $i <= 100; $i += 1) {
    print $result Nx(@sequence_lengths, $i);
};
close $result;

sub Nx {
    my @lengths = \@{$_[0]};
    my $x = $_[1];
    my $cumulative_sum;
    my $total_sum;
    for my $len (@lengths) {
        $total_sum += $len;
    };
    for my $len (sort {$b <=> $a} @lengths) { 
        $cumulative_sum += $len;
        if ($cumulative_sum >= $total_sum * $x / 100) {
            return "$x\t$len\n";
        };
    };
};
