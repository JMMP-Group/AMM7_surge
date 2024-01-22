#!/usr/bin/perl

use strict;
use warnings;

################################################################################
# A Fcm::Keyword::Entry sub-class for testing
{
    package TestEntry;
    use base qw{Fcm::Keyword::Entry};
}

################################################################################
# A mock loader implementing the Fcm::Keyword::Loader interface
{
    package MockLoader0;
    use Scalar::Util qw{blessed};

    ############################################################################
    # Constructor
    sub new {
        my ($class) = @_;
        return bless({number_of_calls_to_load_to => 0}, $class);
    }

    ############################################################################
    ##Returns the package name
    sub get_source {
        my ($self) = @_;
        return blessed($self);
    }

    ############################################################################
    # Returns number of times $self->load_to() has been called
    sub get_number_of_calls_to_load_to {
        my ($self) = @_;
        return $self->{number_of_calls_to_load_to};
    }

    ############################################################################
    # Loads data into $entries, and returns number of entries loaded
    sub load_to {
        my ($self, $entries) = @_;
        $self->{number_of_calls_to_load_to}++;
        return $self->load_to_impl($entries);
    }

    ############################################################################
    # Returns 0
    sub load_to_impl {
        my ($self, $entries) = @_;
        return 0;
    }
}

################################################################################
# A mock loader implementing the Fcm::Keyword::Loader interface
{
    package MockLoader1;
    our @ISA = qw{MockLoader0};

    my %VALUE_OF = (foo => 'foo1', bar => 'bar2', baz => 'baz3');

    ############################################################################
    # Returns a reference to the mock data
    sub get_data {
        my ($class) = @_;
        return \%VALUE_OF;
    }

    ############################################################################
    ##Writes mock data to the $entries object
    sub load_to_impl {
        my ($self, $entries) = @_;
        my $counter = 0;
        for my $key (keys(%{$self->get_data()})) {
            $entries->add_entry($key, $self->get_data()->{$key});
            $counter++;
        }
        return $counter;
    }
}

################################################################################
# A mock loader implementing the Fcm::Keyword::Loader interface
{
    package MockLoader2;
    our @ISA = qw{MockLoader1};

    my %VALUE_OF = (sausages => 'pig', eggs => 'hen', chips => 'potato');

    ############################################################################
    # Returns a reference to the mock data
    sub get_data {
        my ($class) = @_;
        return \%VALUE_OF;
    }
}

package main;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::Keyword::Entries';
    use_ok($class);
    test_empty_constructor($class);
    test_constructor($class);
    test_add_entry($class);
    test_loaders($class);
}

################################################################################
# Tests empty constructor
sub test_empty_constructor {
    my ($class) = @_;
    my $prefix = 'empty constructor';
    my $entries = $class->new();
    isa_ok($entries, $class);
    is($entries->get_entry_class(), 'Fcm::Keyword::Entry',
        "$prefix: default entry class");
    is_deeply([$entries->get_loaders()], [], "$prefix: empty list of loaders");
    is_deeply([$entries->get_all_entries()], [],
        "$prefix: empty list of entries");
    for my $arg ('foo', undef) {
        is($entries->get_entry_by_key($arg), undef,
            "$prefix: entry by key: undef");
        is($entries->get_entry_by_value($arg), undef,
            "$prefix: entry by value: undef");
    }
}

################################################################################
# Tests other constructor usages
sub test_constructor {
    my ($class) = @_;
    my $prefix = 'constructor';
    my @loaders = (MockLoader1->new(), MockLoader2->new());
    my $entries = $class->new({
        entry_class => 'not-a-class',
        loaders     => \@loaders,
    });
    isa_ok($entries, $class);
    is($entries->get_entry_class(), 'not-a-class', "$prefix: entry class");
    is_deeply([$entries->get_loaders()], \@loaders, "$prefix: list of loaders");
    eval {
        $entries->add_entry('key', 'value');
    };
    isnt($@, undef, "$prefix: invalid entry class");
}

################################################################################
# Tests adding entries
sub test_add_entry {
    my ($class) = @_;
    my $prefix = 'add entry';
    my %VALUE_OF = (key1 => 'value1', egg => 'white and yolk', 'xyz.abc' => '');
    for my $entry_class ('Fcm::Keyword::Entry', 'TestEntry') {
        my $entries = $class->new({entry_class => $entry_class});
        my $number_of_entries = 0;
        for my $key (keys(%VALUE_OF)) {
            my $entry = $entries->add_entry($key, $VALUE_OF{$key});
            isa_ok($entry, $entry_class);
            is(scalar(@{$entries->get_all_entries()}), ++$number_of_entries,
                "$prefix: number of entries: $number_of_entries");
        }
        for my $key (keys(%VALUE_OF)) {
            my $entry = $entries->get_entry_by_key($key);
            isa_ok($entry, $entry_class);
            is($entry->get_key(), uc($key), "$prefix: get by key: $key");
            is($entry->get_value(), $VALUE_OF{$key},
                "$prefix: get by key: $key: value");
        }
        for my $key (keys(%VALUE_OF)) {
            my $entry = $entries->get_entry_by_value($VALUE_OF{$key});
            isa_ok($entry, $entry_class);
            is($entry->get_key(), uc($key), "$prefix: get by value: $key");
            is($entry->get_value(), $VALUE_OF{$key},
                "$prefix: get by value: $key: value");
        }
        is($entries->get_entry_by_key('no-such-key'), undef,
            "$prefix: get by key: no-such-key");
        is($entries->get_entry_by_value('no-such-value'), undef,
            "$prefix: get by value: no-such-value");
    }
}

################################################################################
# Tests usage of loaders
sub test_loaders {
    my ($class) = @_;
    my $prefix = "loader";
    my @loaders = (MockLoader0->new(), MockLoader1->new(), MockLoader2->new());
    my $entries = $class->new({loaders => \@loaders});
    for my $loader (@loaders) {
        is($loader->get_number_of_calls_to_load_to(), 0, "$prefix: not loaded");
    }
    for my $key (keys(%{$loaders[1]->get_data()})) {
        my $value = $loaders[1]->get_data()->{$key};
        my $entry = $entries->get_entry_by_key($key);
        is($entry->get_key(), uc($key), "$prefix: by key: $key: key");
        is($entries->get_entry_by_value($value), $entry,
            "$prefix: by value: $key: object");
    }
    is($loaders[0]->get_number_of_calls_to_load_to(), 1,
        "$prefix: loaded once: 0");
    is($loaders[1]->get_number_of_calls_to_load_to(), 1,
        "$prefix: loaded once: 1");
    is($loaders[2]->get_number_of_calls_to_load_to(), 0,
        "$prefix: not loaded: 2");
    for my $key (keys(%{$loaders[2]->get_data()})) {
        my $value = $loaders[2]->get_data()->{$key};
        my $entry = $entries->get_entry_by_key($key);
        is($entry->get_key(), uc($key), "$prefix: by key: $key: key");
        is($entries->get_entry_by_value($value), $entry,
            "$prefix: by value: $key: object");
    }
    is($loaders[0]->get_number_of_calls_to_load_to(), 2,
        "$prefix: loaded once: 0");
    is($loaders[1]->get_number_of_calls_to_load_to(), 1,
        "$prefix: loaded once: 1");
    is($loaders[2]->get_number_of_calls_to_load_to(), 1,
        "$prefix: not loaded: 2");
}

__END__
