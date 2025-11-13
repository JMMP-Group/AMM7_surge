# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Exception;
use overload (q{""} => \&as_string);

use Scalar::Util qw{blessed};

# ------------------------------------------------------------------------------
# Returns true if $e is a blessed instance of this class.
sub caught {
    my ($class, $e) = @_;
    return (blessed($e) && $e->isa($class));
}

# ------------------------------------------------------------------------------
# Constructor
sub new {
    my ($class, $args_ref) = @_;
    return bless(
        {message => q{unknown problem}, ($args_ref ? %{$args_ref} : ())},
        $class,
    );
}

# ------------------------------------------------------------------------------
# Returns a string representation of this exception
sub as_string {
    my ($self) = @_;
    return sprintf("%s: %s\n", blessed($self), $self->get_message());
}

# ------------------------------------------------------------------------------
# Returns the message of this exception
sub get_message {
    my ($self) = @_;
    return $self->{message};
}

1;
__END__

=head1 NAME

Fcm::Exception

=head1 SYNOPSIS

    use Fcm::Exception;
    eval {
        croak(Fcm::Exception->new({message => $message}));
    };
    if ($@) {
        if (Fcm::Exception->caught($@)) {
            print({STDERR} $@);
        }
    }

=head1 DESCRIPTION

This exception is raised when there is a generic problem in FCM.

=head1 METHODS

=over 4

=item $class->caught($e)

Returns true if $e is a blessed instance of this class.

=item $class->new({message=E<gt>$message})

Returns a new instance of this exception. Its first argument must be a
reference to a hash containing the detailed I<message> of the exception.

=item $e->as_string()

Returns a string representation of this exception.

=item $e->get_message()

Returns the detailed message of this exception.

=back

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
