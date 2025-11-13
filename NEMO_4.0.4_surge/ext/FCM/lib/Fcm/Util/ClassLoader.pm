# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Util::ClassLoader;
use base qw{Exporter};

our @EXPORT_OK = qw{load};

use Carp qw{croak};
use Fcm::Exception;

sub load {
    my ($class, $test_method) = @_;
    if (!$test_method) {
        $test_method = 'new';
    }
    if (!UNIVERSAL::can($class, $test_method)) {
        eval('require ' . $class);
        if ($@) {
            croak(Fcm::Exception->new({message => sprintf(
                "%s: class loading failed: %s", $class, $@,
            )}));
        }
    }
    return $class;
}

1;
__END__

=head1 NAME

Fcm::ClassLoader

=head1 SYNOPSIS

    use Fcm::Util::ClassLoader;
    $load_ok = Fcm::Util::ClassLoader::load($class);

=head1 DESCRIPTION

A wrapper for loading a class dynamically.

=head1 FUNCTIONS

=over 4

=item load($class,$test_method)

If $class can call $test_method, returns $class. Otherwise, attempts to
require() $class and returns it. If this fails, croak() with a
L<Fcm::Exception|Fcm::Exception>.

=item load($class)

Shorthand for C<load($class, 'new')>.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::Exception|Fcm::Exception>

The load($class,$test_method) function croak() with this exception if it fails
to load the specified class.

=back

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
