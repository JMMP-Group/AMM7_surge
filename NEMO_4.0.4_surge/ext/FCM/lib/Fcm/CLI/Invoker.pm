# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker;

use Carp qw{croak};
use Fcm::CLI::Exception;

################################################################################
# Constructor
sub new {
    my ($class, $args_ref) = @_;
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Returns the name of the (sub)command as given by the user
sub get_command {
    my ($self) = @_;
    return $self->{command};
}

################################################################################
# Returns a reference to a hash containing the options
sub get_options {
    my ($self) = @_;
    return (wantarray() ? %{$self->{options}} : $self->{options});
}

################################################################################
# Returns a reference to an array containing the arguments
sub get_arguments {
    my ($self) = @_;
    return (wantarray() ? @{$self->{arguments}} : $self->{arguments});
}

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my $message = "command not implemented\n";
    $message .= sprintf("opts:");
    for my $key (sort keys(%{$self->get_options()})) {
        my $value = $self->get_options()->{$key};
        $message .= sprintf(
            " [%s=%s]",
            $key,
            ($value && ref($value) eq 'ARRAY' ? join(q{, }, @{$value}) : $value)
        );
    }
    $message .= sprintf("\n");
    $message .= sprintf("args: [%s]\n", join(q{] [}, $self->get_arguments()));
    croak(Fcm::CLI::Exception->new({message => $message}));
}

1;
__END__

=head1 NAME

Fcm::CLI::Invoker

=head1 SYNOPSIS

    use Fcm::CLI::Invoker;
    $invoker = Fcm::CLI::Invoker->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke();

=head1 DESCRIPTION

This is the base class for an invoker of a FCM sub-system from the CLI.
Sub-classes should override the invoke() method.

=head1 METHODS

=over 4

=item new($args_ref)

Constructor. It accepts a hash reference as an argument. The element I<command>
should be set to the actual (sub)command as specified by the user. The element
I<options> should be a reference to a hash containing the specified command line
options. The element I<arguments> should be a reference to an array containing
the remaining command line arguments.

=item get_command()

Returns the actual (sub)command as specified by the user.

=item get_options()

Returns a hash containing the specified command line options. In scalar context,
returns a reference to the hash.

=item get_arguments()

Returns an array containing the (remaining) command line arguments. In scalar
context, returns a reference to the array.

=item invoke()

Sub-classes should override this method. Calling the method in this base
class causes the system to croak() with a
L<Fcm::CLI::Exception|Fcm::CLI::Exception>.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::CLI::Exception|Fcm::CLI::Exception>

The C<invoke()> croak() with this exception.

=back

=head1 SEE ALSO

L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
