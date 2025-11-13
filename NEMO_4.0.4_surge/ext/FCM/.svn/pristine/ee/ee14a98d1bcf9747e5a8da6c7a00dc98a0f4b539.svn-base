#!/usr/bin/perl
# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../../../lib";

use Test::More (tests => 3);

if (!caller()) {
    main(@ARGV);
}

sub main {
    my $CLASS = 'Fcm::Build::Fortran';
    use_ok($CLASS);
    my $util = $CLASS->new();
    isa_ok($util, $CLASS);
    test_extract_interface($util);
}

sub test_extract_interface {
    my ($util) = @_;
    my $root = ($0 =~ qr{\A(.+)\.t\z}msx)[0];
    my $f90 = $root . '-extract-interface-source.f90';
    my $f90_interface = $root . '-extract-interface-result.f90';
    open(my($handle_for_source), '<', $f90) || die("$f90: $!");
    my @actual_lines = $util->extract_interface($handle_for_source);
    close($handle_for_source);
    open(my($handle_for_result), '<', $f90_interface)
        || die("$f90_interface: $!");
    my @expected_lines = readline($handle_for_result);
    close($handle_for_result);
    is_deeply(\@actual_lines, \@expected_lines, 'extract_interface');
}

__END__
