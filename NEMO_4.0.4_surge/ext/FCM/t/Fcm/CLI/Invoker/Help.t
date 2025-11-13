#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Invoker::Help';
    use_ok($class);
}

# TODO: actual unit tests

__END__
