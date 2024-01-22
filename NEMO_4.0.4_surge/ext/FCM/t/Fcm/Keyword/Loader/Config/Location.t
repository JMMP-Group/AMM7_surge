#!/usr/bin/perl

use strict;
use warnings;

use Fcm::Config;
use Fcm::Keyword::Entries;
use Test::More qw{no_plan};

my %VALUE_OF = (
    foo => 'fcm-test://foo/foo',
    bar => 'fcm-test://bar/bar',
    baz => 'fcm-test://baz/baz',
);

main();

sub main {
    my $class = 'Fcm::Keyword::Loader::Config::Location';
    use_ok($class);
    test_constructor($class);
    test_load_to($class);
}

################################################################################
# Tests simple usage of the constructor
sub test_constructor {
    my ($class) = @_;
    my $prefix = "constructor";
    my $loader = $class->new();
    isa_ok($loader, $class);
    is($loader->get_source(), 'Fcm::Config', "$prefix: get_source()");
}

################################################################################
# Tests loading to an Fcm::Keyword::Entries object
sub test_load_to {
    my ($class) = @_;
    my $prefix = 'load to';
    my $config = Fcm::Config->instance();
    for my $key (keys(%VALUE_OF)) {
        $config->setting(['URL', $key], $VALUE_OF{$key});
    }
    my $loader = $class->new();
    my $entries = Fcm::Keyword::Entries->new({
        entry_class => 'Fcm::Keyword::Entry::Location',
    });
    isnt($loader->load_to($entries), 0, "$prefix: number loaded");
    for my $key (keys(%VALUE_OF)) {
        my $entry = $entries->get_entry_by_key($key);
        if ($entry) {
            is($entry->get_key(), uc($key), "$prefix: by key: $key");
            is($entry->get_value(), $VALUE_OF{$key}, "$prefix: by value: $key");
            is(
                $entries->get_entry_by_value($VALUE_OF{$key}),
                $entry,
                "$prefix: by key: $key: object",
            );
        }
        else {
            fail("$prefix: by key: $key");
        }
    }
}

# TODO: tests loading of browser mapping

__END__
