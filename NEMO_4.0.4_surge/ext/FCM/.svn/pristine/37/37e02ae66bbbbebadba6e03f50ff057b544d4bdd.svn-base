# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Config;

use Carp;
use Fcm::Keyword::Entries;
use Fcm::Keyword::Exception;
use Fcm::Util::ClassLoader;

our %CONFIG_OF = (
    LOCATION_ENTRIES => {
        entry_class => 'Fcm::Keyword::Entry::Location',
        loaders => [
            {
                class => 'Fcm::Keyword::Loader::Config::Location',
            },
        ],
    },
    REVISION_ENTRIES => {
        entry_class => 'Fcm::Keyword::Entry',
        loaders => [
            {
                class   => 'Fcm::Keyword::Loader::Config::Revision',
                options => [{key => 'namespace', valuekey => 'key'}],
            },
            {
                class   => 'Fcm::Keyword::Loader::VC::Revision',
                options => [{key => 'source', valuekey => 'value'}],
            },
        ],
    },
);

################################################################################
# Returns a Fcm::Keyword::Entries object for given configuration
sub get_entries {
    my ($context, $args_ref) = @_;
    if (!exists($CONFIG_OF{$context})) {
        croak(Fcm::Keyword::Exception->new({message => sprintf(
            "%s: keyword configuration not found", $context,
        )}));
    }
    my $config_ref = $CONFIG_OF{$context};
    my @loaders;
    if (exists($config_ref->{loaders})) {
        for my $loader_config (@{$config_ref->{loaders}}) {
            my $class = $loader_config->{class};
            Fcm::Util::ClassLoader::load($class);
            my %options;
            if (exists($loader_config->{options})) {
                for my $option_ref (@{$loader_config->{options}}) {
                    my $key = $option_ref->{key};
                    my $value;
                    if (exists($option_ref->{value})) {
                        $value = $option_ref->{value};
                    }
                    elsif (
                           exists($option_ref->{valuekey})
                        && $args_ref
                        && ref($args_ref) eq 'HASH'
                        && exists($args_ref->{$option_ref->{valuekey}})
                    ) {
                        $value = $args_ref->{$option_ref->{valuekey}};
                    }
                    $options{$key} = $value;
                }
            }
            push @loaders, $class->new(\%options);
        }
    }
    my %entries_options = (
        (@loaders ? (loaders => \@loaders) : ()),
        (
            exists($config_ref->{entry_class})
            ? (entry_class => $config_ref->{entry_class})
            : ()
        ),
    );
    return Fcm::Keyword::Entries->new(\%entries_options);
}

1;
__END__

=head1 NAME

Fcm::Keyword::Config

=head1 SYNOPSIS

    use Fcm::Keyword::Config;

=head1 DESCRIPTION

This module stores the default configuration used by modules in the
L<Fcm::Keyword> family.

=head1 FUNCTIONS

=over 4

=item get_entries($context,$args_ref)

Returns a L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object for a given
$context. If there is no matching $context in the configuration, croak() with a
L<Fcm::Keyword::Exception|Fcm::Keyword::Exception>. $args_ref is an optional
argument, which should be a reference to a hash containing a I<key> and a
I<value> element. It can be used by this function to set up the constructor
options in the loaders of the returned
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object.

=back

=head1 DIAGNOSTICS

=head1 TO DO

Allow configuration to be changed in runtime.

Convert this module to OO?

Separate configuration from logic if this module becomes any bigger.

Unit tests.

=head1 SEE ALSO

L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Exception|Fcm::Keyword::Exception>,
L<Fcm::Keyword::Entry::Location|Fcm::Keyword::Entry::Location>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
