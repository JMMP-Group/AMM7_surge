# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Interactive::InputGetter;

use Carp qw{croak};

################################################################################
# Constructor
sub new {
    my ($class, $args_ref) = @_;
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Methods: get_*
for my $key (
    ############################################################################
    # Returns the title of the prompt
    'title',
    ############################################################################
    # Returns the message of the prompt
    'message',
    ############################################################################
    # Returns the of the prompt
    'type',
    ############################################################################
    # Returns the default return value
    'default',
) {
    no strict qw{refs};
    my $getter = "get_$key";
    *$getter = sub {
        my ($self) = @_;
        return $self->{$key};
    }
}

################################################################################
# Invokes the getter
sub invoke {
    my ($self) = @_;
    croak("Fcm::Interactive::InputGetter->invoke() not implemented.");
}

1;
__END__

=head1 NAME

Fcm::Interactive::TxtInputGetter

=head1 SYNOPSIS

    use Fcm::Interactive::TxtInputGetter;
    $answer = Fcm::Interactive::get_input(
        title   => 'My title',
        message => 'Would you like to ...?',
        type    => 'yn',
        default => 'n',
    );

=head1 DESCRIPTION

An object of this abstract class is used by
L<Fcm::Interactive|Fcm::Interactive> to get a user reply.

=head1 METHODS

=over 4

=item new($args_ref)

Constructor, normally invoked via L<Fcm::Interactive|Fcm::Interactive>.

Input options are: I<title>, for a short title of the prompt, I<message>, for
the message prompt, I<type> for the prompt type, and I<default> for the default
value of the return value.

Prompt type can be YN (yes or no), YNA (yes, no or all) or input (for an input
string).

=item get_title()

Returns the title of the prompt.

=item get_message()

Returns the message of the prompt.

=item get_type()

Returns the type of the prompt, can be YN (yes or no), YNA (yes, no or all) or
input (for an input string).

=item get_default()

Returns the default return value of invoke().

=item invoke()

Gets an input string from the user, and returns it. Sub-classes must override
this method.

=back

=head1 SEE ALSO

L<Fcm::Interactive|Fcm::Interactive>,
L<Fcm::Interactive::TxtInputGetter|Fcm::Interactive::TxtInputGetter>,
L<Fcm::Interactive::GUIInputGetter|Fcm::Interactive::GUIInputGetter>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
