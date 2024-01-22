#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Option';
    use_ok($class);
    test_simplest($class);
    test_simplest_scalar_arg($class);
    test_simplest_array_arg($class);
    test_simplest_hash_arg($class);
    test_simple($class);
    test_simple_scalar_arg($class);
    test_simple_array_arg($class);
    test_simple_hash_arg($class);
    test_long_letter($class);
}

################################################################################
# Tests simplest usage
sub test_simplest {
    my ($class) = @_;
    my $prefix = 'simplest';
    my $option = $class->new({
        delimiter   => 'delimiter-value',
        description => 'description value',
        name        => 'name-value',
    });
    isa_ok($option, $class);
    is($option->get_delimiter(), 'delimiter-value', "$prefix: delimiter");
    is($option->get_description(), 'description value', "$prefix: description");
    is($option->get_name(), 'name-value', "$prefix: name");
    is($option->get_letter(), undef, "$prefix: letter");
    ok(!$option->has_arg(), "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value', "$prefix: has arg");
}

################################################################################
# Tests simplest usage with a scalar argument
sub test_simplest_scalar_arg {
    my ($class) = @_;
    my $prefix = 'simplest scalar arg';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        has_arg     => $class->SCALAR_ARG,
    });
    isa_ok($option, $class);
    is($option->has_arg(), $class->SCALAR_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value=s', "$prefix: has arg");
}

################################################################################
# Tests simplest usage with array argument
sub test_simplest_array_arg {
    my ($class) = @_;
    my $prefix = 'simplest array arg';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        has_arg     => $class->ARRAY_ARG,
    });
    isa_ok($option, $class);
    is($option->has_arg(), $class->ARRAY_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value=s@', "$prefix: has arg");
}

################################################################################
# Tests simplest usage with hash argument
sub test_simplest_hash_arg {
    my ($class) = @_;
    my $prefix = 'simplest hash arg';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        has_arg     => $class->HASH_ARG,
    });
    isa_ok($option, $class);
    is($option->has_arg(), $class->HASH_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value=s%', "$prefix: has arg");
}

################################################################################
# Tests simple usage
sub test_simple {
    my ($class) = @_;
    my $prefix = 'simple';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        letter      => 'n',
    });
    isa_ok($option, $class);
    is($option->get_description(), 'description value', "$prefix: description");
    is($option->get_name(), 'name-value', "$prefix: name");
    is($option->get_letter(), 'n', "$prefix: letter");
    is($option->has_arg(), $class->NO_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value|n', "$prefix: has arg");
}

################################################################################
# Tests simple usage with a scalar argument
sub test_simple_scalar_arg {
    my ($class) = @_;
    my $prefix = 'simple scalar arg';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        letter      => 'n',
        has_arg     => $class->SCALAR_ARG,
    });
    isa_ok($option, $class);
    is($option->has_arg(), $class->SCALAR_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value|n=s', "$prefix: has arg");
}

################################################################################
# Tests simplest usage with array argument
sub test_simple_array_arg {
    my ($class) = @_;
    my $prefix = 'simple array arg';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        letter      => 'n',
        has_arg     => $class->ARRAY_ARG,
    });
    isa_ok($option, $class);
    is($option->has_arg(), $class->ARRAY_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value|n=s@', "$prefix: has arg");
}

################################################################################
# Tests simplest usage with hash argument
sub test_simple_hash_arg {
    my ($class) = @_;
    my $prefix = 'simple hash arg';
    my $option = $class->new({
        description => 'description value',
        name        => 'name-value',
        letter      => 'n',
        has_arg     => $class->HASH_ARG,
    });
    isa_ok($option, $class);
    is($option->has_arg(), $class->HASH_ARG, "$prefix: has arg");
    is($option->get_arg_for_getopt_long(), 'name-value|n=s%', "$prefix: has arg");
}

################################################################################
# Tests longer than 1 letter
sub test_long_letter {
    my ($class) = @_;
    my $prefix = 'long letter';
    my $option = $class->new({
        name   => 'name-value',
        letter => 'name',
    });
    isa_ok($option, $class);
    is($option->get_letter(), 'n', "$prefix: letter");
    is($option->get_arg_for_getopt_long(), 'name-value|n', "$prefix: has arg");
}

__END__
