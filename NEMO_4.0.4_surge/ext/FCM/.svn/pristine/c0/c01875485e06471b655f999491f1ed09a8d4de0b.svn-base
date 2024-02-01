#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my ($class) = 'Fcm::Keyword::Entry::Location';
    use_ok($class);
    test_constructor($class);
}

################################################################################
# Tests constructor
sub test_constructor {
    my ($class) = @_;
    my $prefix = 'constructor';
    isa_ok($class->new(), $class, "$prefix: empty");
    my $entry = $class->new({key => 'key', value => 'value'});
    isa_ok($entry, $class, "$prefix: normal");
    is($entry->get_key(), 'key', "$prefix: normal: get_key()");
    is($entry->get_value(), 'value', "$prefix: normal: get_value()");
    isa_ok($entry->get_revision_entries(), 'Fcm::Keyword::Entries');
}

__END__
