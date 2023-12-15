use strict;
use warnings;
use Getopt::Long;

GetOptions(
    "work|w=i" => \my $work,
    "break|b=i" => \my $break,
);

if ($work < $break || $work > 60 || $break > 60) {
    die("Not appropriate work/break timings!");
}

sub clear_terminal {
    system("clear");
}

sub display_time {
    my $x = $_[0];
    my $div = int($x / 60);
    my $mod = $x % 60;
    my $minutes = ($div < 10) ? "0$div" : $div;
    my $seconds = ($mod < 10) ? "0$mod" : $mod;

    print "$minutes:$seconds\n";
}

sub run {
    my $x = 0;

    while ($x < ($work * 60)) {
        clear_terminal();
        print "Work ($work minutes):\n";
        display_time($x);
        sleep(1);
        $x += 1;
    }
}

sub pause {
    my $break = $_[0];
    my $x = 0;

    while ($x < ($break * 60)) {
        clear_terminal();
        print "Break ($break minutes):\n";
        display_time($x);
        sleep(1);
        $x += 1;
    }
}

my $count_run = 0;

while (1) {
    run();
    $count_run += 1;
    if (($count_run > 0) && ($count_run % 4 == 0)) {
        pause($break * 3);
    } else {
        pause($break);
    }
}
