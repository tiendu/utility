use strict;
use warnings;

my @files = @ARGV;
my %records;

for my $file (@files) {
    $records{$file} = [];
    open my $input, "<:utf8", $file or die;
    print "Reading $file\n";
    while (<$input>) {
        chomp;
        if ($_ =~ m/\A>(.+)/) {
            next;
        } else {
            push $records{$file}->@*, length($_);
        };
        last if eof;
    };
    close $input;
};
print "Calculating. Please wait...\n";
for my $record (sort keys %records) {
    my $size = format_number(scalar $records{$record}->@*);
    my @sorted = sort {$a <=> $b} $records{$record}->@*;
    my ($min_length, $max_length) = @sorted[0, -1];
    my ($sum_length, $sum_squared_length) = (0, 0);
    $sum_length += $_ foreach @sorted;
    $sum_squared_length += ($_ ** 2) foreach @sorted;
    my $auN = format_number(sprintf("%.3f", $sum_squared_length / $sum_length));
    $min_length = format_number($min_length);
    $max_length = format_number($max_length);
    $sum_length = format_number(sprintf("%.3f", $sum_length / 1_000_000));
    my $n50 = format_number(Nx(\@sorted, 50));
    my $n90 = format_number(Nx(\@sorted, 90));
    print "=" x 36 . "\n";
    my $filename = $record =~ s/\A(.+)\///r;
    print join(" | ", $filename, "Seqs: ${size}", "Total (Mb): ${sum_length}", "Minlen (bp): ${min_length}", "Maxlen (bp): ${max_length}", "auN: ${auN}", "N50: ${n50}", "N90: ${n90}", "\n");
};

# for my $record (sort keys %records) {
#     open my $result, ">:utf8", "${record}_Nx.tsv" or die;
#     for (my $i = 0; $i <= 100; $i += 1) {
#         print $result join("\t", $i, Nx($records{$record}, $i), "\n");
#     };
#     close $result;
# };

sub Nx {
    my @lengths = @{$_[0]};
    my $x = $_[1];
    my $cumulative_sum;
    my $total_sum;
    for my $len (@lengths) {
        $total_sum += $len;
    };
    for my $len (sort {$a <=> $b} @lengths) { 
        $cumulative_sum += $len;
        if ($cumulative_sum >= $total_sum * (100 - $x) / 100) {
            return $len;
        };
    };
};

sub format_number {
    my $num = $_[0];
    while ($num =~ s/^(-?\d+)(\d{3}(?:,\d{3})*(?:\.\d+)*)$/$1,$2/) {};
    return $num;
};

# sub table {
#     my ($record, $length, $sum_length, $min_length, $max_length, $auN, $n50, $n90) = @_;
#     
# #     $~ = "FMT1_TOP";
# #     write;
#     format_name STDOUT "FMT1";
#     format_top_name STDOUT "FMT1_TOP";
#     write;
# 
# format FMT1_TOP =
#     File  | No of seqs |   Total bp   |    Min  |     Max    |    auN     |    n50   |   n90  |
# ----------+------------+--------------+---------+------------+------------+----------+--------+
# .
# 
# format FMT1 =
# @<<<<<<<< | @>>>>>>>>> | @>>>>>>>>>>> | @>>>>>> | @>>>>>>>>> | @>>>>>>>>> | @>>>>>>> | @>>>>> |
# $record, $length, $sum_length, $min_length, $max_length, $auN, $n50, $n90
# .
# }
