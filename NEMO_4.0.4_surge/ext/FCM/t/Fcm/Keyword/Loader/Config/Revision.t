#!/usr/bin/perl

use strict;
use warnings;

use Fcm::Config;
use Fcm::Keyword::Entries;
use Test::More qw{no_plan};

my %VALUE_OF = (
    bar => {
        'bar3'    => 3,
        'bar3.1'  => 31,
        'bar3.14' => 314,
    },
    baz => {
        'bear'   => 4,
        'bee'    => 6,
        'spider' => 8,
    },
);

main();

sub main {
    my $class = 'Fcm::Keyword::Loader::Config::Revision';
    use_ok($class);
    test_constructor($class);
    test_load_to($class);
}

################################################################################
# Tests simple usage of the constructor
sub test_constructor {
    my ($class) = @_;
    my $prefix = "constructor";
    my $loader = $class->new({namespace => 'namespace'});
    isa_ok($loader, $class);
    is($loader->get_namespace(), 'namespace', "$prefix: get_namespace()");
    is($loader->get_source(), 'Fcm::Config', "$prefix: get_source()");
}

################################################################################
# Tests loading to an Fcm::Keyword::Entries object
sub test_load_to {
    my ($class) = @_;
    my $prefix = 'load to';
    my $config = Fcm::Config->instance();
    for my $key (keys(%VALUE_OF)) {
        for my $rev_key (keys(%{$VALUE_OF{$key}})) {
            my $value = $VALUE_OF{$key}{$rev_key};
            $config->setting(['URL_REVISION', uc($key), uc($rev_key)], $value);
        }
        my $entries = Fcm::Keyword::Entries->new();
        my $loader = $class->new({namespace => $key});
        isnt($loader->load_to($entries), 0, "$prefix: number loaded");
        for my $rev_key (keys(%{$VALUE_OF{$key}})) {
            my $entry = $entries->get_entry_by_key($rev_key);
            my $value = $VALUE_OF{$key}{$rev_key};
            if ($entry) {
                is(
                    $entry->get_key(),
                    uc($rev_key),
                    "$prefix: by key: $rev_key",
                );
                is($entry->get_value(), $value, "$prefix: by value: $rev_key");
                is(
                    $entries->get_entry_by_value($value),
                    $entry,
                    "$prefix: by key: $key: object",
                );
            }
            else {
                fail("$prefix: by key: $rev_key");
            }
        }
    }
}

__END__
