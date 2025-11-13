#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $module = 'Fcm::Keyword::Config';
    use_ok($module);
}

__END__
