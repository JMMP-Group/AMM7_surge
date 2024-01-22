#!/usr/bin/perl

use strict;
use warnings;

use Fcm::CLI::Option;
use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Subcommand';
    use_ok($class);
    test_constructor($class);
    test_has_a_name($class);
    test_as_string($class);
}

################################################################################
# Tests the constructor
sub test_constructor {
    my ($class) = @_;
    my $prefix = 'constructor';
    my %OPTION_OF = (
        description    => 'description value',
        invoker_class  => 'invoker_class value',
        invoker_config => 'invoker_config value',
        is_vc          => 'is_vc value',
        names          => 'names value',
        options        => 'options value',
        synopsis       => 'synopsis value',
        usage          => 'usage value',
    );
    my $subcommand = Fcm::CLI::Subcommand->new(\%OPTION_OF);
    isa_ok($subcommand, $class, $prefix);
    for my $key (keys(%OPTION_OF)) {
        my $getter = index($key, 'is') == 0 ? $key : "get_$key";
        is($subcommand->$getter(), $OPTION_OF{$key}, "$prefix: $getter");
    }
}

################################################################################
# Tests match a string name to a subcommand
sub test_has_a_name {
    my ($class) = @_;
    my $prefix = 'has a name';
    my @NAMES = ('foo', 'bar', 'baz');
    my $subcommand = $class->new({names => \@NAMES});
    for my $name (@NAMES) {
        ok($subcommand->has_a_name($name), "$prefix: $name");
    }
    for my $name (qw{egg ham mayo}) {
        ok(!$subcommand->has_a_name($name), "$prefix: $name");
    }
}

################################################################################
# Tests string representation of a subcommand
sub test_as_string {
    my ($class) = @_;
    my $prefix = 'as string';
    my %OPTION_OF = (
        'foo (bar, baz)' => ['foo', 'bar', 'baz'],
        'foo (bar)'      => ['foo', 'bar'],
        'foo'            => ['foo'],
        q{}              => [],
    );
    for my $key (keys(%OPTION_OF)) {
        my $subcommand = $class->new({names => $OPTION_OF{$key}});
        is($subcommand->as_string(), $key, "$prefix: $key");
    }
}

__END__
