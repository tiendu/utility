use strict;

my @files = @ARGV;
my %records;

for my $file (@files) {
    my @lengths;
    my $file_name = $1 if $file =~ /^(.+)\.([^.]+)$/;
    open my $file, "<:utf8", $file or die;
    while (<$file>) {
        chomp;
        if ($_ =~ m/\A>(.+)/) {
            next;
        } else {
            push @lengths, length($_);
        };
        $records{$file_name} = [@lengths];
        last if eof;
    };
    close $file;
};

print "File\tNo of seqs\tTotal bp\tMin\tMax\tauN\tN50\tN90\n";
for my $record (sort keys %records) {
    my ($sum_length, $sum_squared_length) = (0, 0);
    $sum_length += $_ foreach $records{$record}->@*;
    my $length = scalar $records{$record}->@* + 1;
    my $max_length = max($records{$record}->@*);
    my $min_length = min($records{$record}->@*);
    $sum_squared_length += ($_ ** 2) foreach $records{$record}->@*;
    my $auN = sprintf("%.3f", $sum_squared_length / $sum_length);
    my $n50 = Nx(\$records{$record}->@*, 50);
    my $n90 = Nx(\$records{$record}->@*, 90);
    print join("\t", $record, $length, $sum_length, $min_length, $max_length, $auN, $n50, $n90, "\n");
};

# for my $record (sort keys %records) {
#     open my $result, ">:utf8", "${record}_Nx.tsv" or die;
#     for (my $i = 0; $i <= 100; $i += 1) {
#         print $result join("\t", $i, Nx($records{$record}, $i), "\n");
#     };
#     close $result;
# };


sub max {
    my $max = shift;
    foreach (@_) {
        $max = $_ if $_ > $max;
    };
    return $max;
}

sub min {
    my $min = shift;
    foreach (@_) {
        $min = $_ if $_ < $min;
    };
    return $min;
};

sub Nx {
    my @lengths = @{$_[0]};
    my $x = $_[1];
    my $cumulative_sum;
    my $total_sum;
    for my $len (@lengths) {
        $total_sum += $len;
    };
    for my $len (sort {$b <=> $a} @lengths) { 
        $cumulative_sum += $len;
        if ($cumulative_sum >= $total_sum * $x / 100) {
            return $len;
        };
    };
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
