#!/usr/bin/perl
# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

use strict;
use warnings;

use Fcm::CfgLine;
use Fcm::Config;
use Scalar::Util qw{reftype};
use Test::More (tests => 90);

BEGIN: {
    use_ok('Fcm::ConfigSystem');
}

my $CONFIG = undef;

# ------------------------------------------------------------------------------
if (!caller()) {
    main(@ARGV);
}

# ------------------------------------------------------------------------------
sub main {
    local @ARGV = @_;
    test_compare_setting_in_config();
}

# ------------------------------------------------------------------------------
# Tests "compare_setting_in_config".
sub test_compare_setting_in_config {
    my $PREFIX = 'TEST';
    my %S = (egg => [qw{boiled poached}], ham => 'roasted', bacon => 'fried');
    my %S_MOD = (ham => 'boiled');
    my %S_MOD_ARRAY = (egg => [qw{scrambled omelette}]);
    my %S_ADD = (mushroom => 'sauteed');
    my %S_DEL = (bacon => undef);

    my @ITEMS = (
        {
            name     => 'empty',
            original => {},
            added    => {},
            removed  => {},
            modified => {},
        },
        {
            name     => 'add keys to empty',
            original => {},
            added    => {%S},
            removed  => {},
            modified => {%S},
        },
        {
            name     => 'remove all',
            original => {%S},
            added    => {},
            removed  => {},
            modified => {map {($_, undef)} keys(%S)},
        },
        {
            name     => 'no change',
            original => {%S},
            added    => {%S},
            removed  => {},
            modified => {},
        },
        {
            name     => 'modify key',
            original => {%S},
            added    => {%S, %S_MOD},
            removed  => {},
            modified => {%S_MOD},
        },
        {
            name     => 'modify an array key',
            original => {%S},
            added    => {%S, %S_MOD_ARRAY},
            removed  => {},
            modified => {%S_MOD_ARRAY},
        },
        {
            name     => 'add a key',
            original => {%S},
            added    => {%S, %S_ADD},
            removed  => {},
            modified => {%S_ADD},
        },
        {
            name     => 'delete a key',
            original => {%S},
            added    => {%S},
            removed  => {%S_DEL},
            modified => {%S_DEL},
        },
        {
            name     => 'modify a key and delete a key',
            original => {%S},
            added    => {%S, %S_MOD},
            removed  => {%S_DEL},
            modified => {%S_MOD, %S_DEL},
        },
        {
            name     => 'add a key and delete a key',
            original => {%S},
            added    => {%S, %S_ADD},
            removed  => {%S_DEL},
            modified => {%S_ADD, %S_DEL},
        },
    );

    # A naive function to serialise an array reference
    my $flatten = sub {
        if (ref($_[0]) && reftype($_[0]) eq 'ARRAY') {
            join(q{ }, sort(@{$_[0]}))
        }
        else {
            $_[0];
        }
    };

    my $CONFIG = Fcm::Config->instance();
    for my $item (@ITEMS) {
        # New settings
        $CONFIG->{setting}{$PREFIX} = {%{$item->{added}}};
        for my $key (keys(%{$item->{removed}})) {
            delete($CONFIG->{setting}{$PREFIX}{$key});
        }

        # Old lines
        my @old_lines = map {
            Fcm::CfgLine->new(
                LABEL => $PREFIX . $Fcm::Config::DELIMITER . $_,
                VALUE => $flatten->($item->{original}{$_}),
            )
        } keys(%{$item->{original}});

        # Invokes the method
        my $system = Fcm::ConfigSystem->new();
        my ($changed_hash_ref, $new_cfg_lines_ref)
            = $system->compare_setting_in_config($PREFIX, \@old_lines);

        # Tests the return values
        my $T = $item->{name};
        is_deeply(
            $changed_hash_ref, $item->{modified},
            "$T: \$changed_hash_ref content",
        );
        is(
            scalar(@{$new_cfg_lines_ref}),
            scalar(keys(%{$item->{added}})) - scalar(keys(%{$item->{removed}})),
            "$T: \$new_cfg_lines_ref length",
        );
        for my $line (@{$new_cfg_lines_ref}) {
            my $key = $line->label_from_field(1);
            ok(exists($item->{added}{$key}), "$T: expected label $key");
            ok(!exists($item->{removed}{$key}), "$T: unexpected label $key");
            is(
                $line->value(), $flatten->($item->{added}{$key}),
                "$T: line content $key",
            );
        }
    }
}

__END__
