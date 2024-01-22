# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Subcommand;

use Carp qw{croak};
use Fcm::CLI::Exception;
use Fcm::Util::ClassLoader;

################################################################################
# Constructor
sub new {
    my ($class, $args_ref) = @_;
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Methods: get_*
for my $key (
    # Returns the long description of this subcommand
    'description',
    # Returns the class of the invoker of this subcommand
    'invoker_class',
    # Returns the extra config to be given to the invoker of this subcommand
    'invoker_config',
    # Returns the names of this subcommand
    'names',
    # Returns the options of this subcommand
    'options',
    # Returns the synopsis of this subcommand
    'synopsis',
    # Returns the usage of this subcommand
    'usage',
) {
    no strict qw{refs};
    my $getter = "get_$key";
    *$getter = sub {
        my ($self) = @_;
        if (defined($self->{$key})) {
            if (ref($self->{$key}) eq 'ARRAY') {
                return (wantarray() ? @{$self->{$key}} : $self->{$key});
            }
            else {
                return $self->{$key};
            }
        }
        else {
            return;
        }
    }
}

################################################################################
# Returns true if this subcommand represents a command in the CM sub-system
sub is_vc {
    my ($self) = @_;
    return $self->{is_vc};
}

################################################################################
# Returns true if $string matches a name of this subcommand
sub has_a_name {
    my ($self, $string) = @_;
    if ($self->get_names() && ref($self->get_names()) eq 'ARRAY') {
        my %name_of = map {$_, 1} @{$self->get_names()};
        return exists($name_of{$string});
    }
    else {
        return;
    }
}

################################################################################
# Invokes this subcommand based on current @ARGV
sub get_invoker {
    my ($self, $command) = @_;
    my %options = ();
    if (($self->get_options())) {
        my $problem = q{};
        local($SIG{__WARN__}) = sub {
            ($problem) = @_;
        };
        my $success = GetOptions(
            \%options,
            (map {$_->get_arg_for_getopt_long()} ($self->get_options())),
        );
        if (!$success) {
            croak(Fcm::CLI::Exception->new({message => sprintf(
                "%s: option parse failed: %s", $command, $problem,
            )}));
        }
    }
    my $invoker_class
        = $self->get_invoker_class() ? $self->get_invoker_class()
        :                              'Fcm::CLI::Invoker'
        ;
    Fcm::Util::ClassLoader::load($invoker_class);
    my $invoker = $invoker_class->new({
        command   => $command,
        options   => \%options,
        arguments => [@ARGV],
    });
    return $invoker;
}

################################################################################
# Returns a simple string representation of this subcommand
sub as_string {
    my ($self) = @_;
    # FIXME: can do with using Text::Template or Perl6::Form
    if (
           $self->get_names()
        && ref($self->get_names()) eq 'ARRAY'
        && @{$self->get_names()}
    ) {
        my @names = $self->get_names();
        my $return = $names[0];
        for my $i (1 .. $#names) {
            if ($names[$i]) {
                $return
                    .= $i == 1       ? q{ (} . $names[$i]
                    :                  q{, } . $names[$i]
                    ;
            }
            if ($i == $#names) {
                $return .= q{)};
            }
        }
        return $return;
    }
    else {
        return q{};
    }
}

1;
__END__

=head1 NAME

Fcm::CLI::Subcommand

=head1 SYNOPSIS

    use Fcm::CLI::Subcommand;
    $subcommand = Fcm::CLI::Subcommand->new({
        names          => ['build', 'bld'],
        options        => [
            Fcm::CLI::Option->new(
                # ... some arguments ...
            ),
            # ... more options
        ],
        synopsis       => 'invokes the build system',
        description    => $description,
        usage          => '[OPTIONS] [CONFIG]',
        invoker_class  => $invoker_class,
        invoker_config => {
            option1 => $option1,
            # ... more options
        },
    });
    $boolean = $subcommand->has_a_name($string);
    $invoker_class = $subcommand->get_invoker_class();

=head1 DESCRIPTION

An object of this class is used to store the configuration of a subcommand of
the FCM CLI.

=head1 METHODS

=over 4

=item new($args_ref)

Constructor.

=item get_description()

Returns the long description of this subcommand.

=item get_invoker_class()

Returns the invoker class of this subcommand, which should be a sub-class of
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>.

=item get_invoker_cconfig()

Returns a reference to a hash containing the extra configuration to be given to
the constructor of the invoker of this subcommand.

=item get_names()

Returns an array containing the names of this subcommand.

=item get_options()

Returns an array containing the options of this subcommand. Each element of
the array should be a L<Fcm::CLI::Option|Fcm::CLI::Option> object.

=item get_synopsis()

Returns a short synopsis of this subcommand.

=item get_usage()

Returns a short usage statement of this subcommand.

=item is_vc()

Returns true if this subcommand represents commands in the underlying VC system.

=item has_a_name($string)

Returns true if a name in C<$self-E<gt>get_names()> matches $string.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::CLI::Exception|Fcm::CLI::Exception>

The invoke() method may croak() with this exception.

=back

=head1 SEE ALSO

L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Option|Fcm::CLI::Option>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
