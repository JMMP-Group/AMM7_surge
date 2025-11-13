#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Config::Default';
    use_ok($class);
}

__END__
