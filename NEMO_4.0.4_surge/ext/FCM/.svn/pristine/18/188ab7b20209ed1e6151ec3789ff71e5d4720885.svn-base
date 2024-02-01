# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Option;

use constant NO_ARG     => 0;
use constant SCALAR_ARG => 1;
use constant ARRAY_ARG  => 2;
use constant HASH_ARG   => 3;
use constant ARG_STRING_SUFFIX_FOR => (q{}, q{=s}, q{=s@}, q{=s%});

################################################################################
# Constructor
sub new {
    my ($class, $args_ref) = @_;
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Methods: get_*
for my $key (
    # Returns the delimiter of this option, if it is an array
    'delimiter',
    # Returns the description of this option
    'description',
    # Returns the (long) name of this option
    'name',
) {
    no strict qw{refs};
    my $getter = "get_$key";
    *$getter = sub {
        my ($self) = @_;
        return $self->{$key};
    }
}

################################################################################
# Returns the letter of this option
sub get_letter {
    my ($self) = @_;
    if (defined($self->{letter})) {
        return substr($self->{letter}, 0, 1);
    }
    else {
        return;
    }
}

################################################################################
# Returns whether the current option has no, scalar, array or hash arguments
sub has_arg {
    my ($self) = @_;
    return (defined($self->{has_arg}) ? $self->{has_arg} : $self->NO_ARG);
}

################################################################################
# Returns true if this option is associated with help
sub is_help {
    my ($self) = @_;
    return $self->{is_help};
}

################################################################################
# Returns an option string/reference pair for Getopt::Long::GetOptions
sub get_arg_for_getopt_long {
    my ($self) = @_;
    my $option_string
        = $self->get_name()
          . ($self->get_letter() ? q{|} . $self->get_letter() : q{})
          . (ARG_STRING_SUFFIX_FOR)[$self->has_arg()]
          ;
    return $option_string;
}

1;
__END__

=head1 NAME

Fcm::CLI::Option

=head1 SYNOPSIS

    use Fcm::CLI::Option;
    $option = Fcm::CLI::Option->new({
        name        => 'name',
        letter      => 'n',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
        is_help     => 1,
        description => 'an example option',
    });

    # time passes ...
    use Getopt::Long qw{GetOptions};
    $success = GetOptions(
        \%hash,
        $option->get_arg_for_getopt_long(), # ('name|n=s')
        # and other options ...
    );
    $option_value = $option->get_value();

=head1 DESCRIPTION

An object of this class represents a CLI option.

=head1 METHODS

=over 4

=item new($args_ref)

Constructor.

=item get_arg_for_getopt_long()

Returns an option string for this option that is suitable for use as arguments
to L<Getopt::Long|Getopt::Long>.

=item get_description()

Returns a description of this option.

=item get_delimiter()

Returns the delimiter of this option. This is only relevant if has_arg() is
equal to C<ARRAY_ARG>. If set, the argument for this option should be re-grouped
using this delimiter.

=item get_name()

Returns the (long) name of this option.

=item get_letter()

Returns the option letter of this option.

=item has_arg()

Returns whether this option has no, scalar, array or hash arguments. See
L</CONSTANTS> for detail.

=item is_help()

Returns true if this option is associated with help.

=back

=head1 CONSTANTS

=over 4

=item NO_ARG

An option has no argument. (Default)

=item SCALAR_ARG

An option has a single scalar argument.

=item ARRAY_ARG

An option has multiple arguments, which can be placed in an array.

=item HASH_ARG

An option has multiple arguments, which can be placed in an hash.

=back

=head1 SEE ALSO

L<Getopt::Long|Getopt::Long>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
