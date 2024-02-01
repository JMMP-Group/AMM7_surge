# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Entry::Location;
use base qw{Fcm::Keyword::Entry};

use Fcm::Keyword::Config;

sub new {
    my ($class, $args_ref) = @_;
    if (!$args_ref) {
        $args_ref = {};
    }
    $args_ref = {
        browser_rev_template       => undef,
        browser_url_template       => undef,
        implied_entry_list         => [],
        is_implied                 => 0,
        location_component_pattern => undef,
        revision_entries           => Fcm::Keyword::Config::get_entries(
            'REVISION_ENTRIES', $args_ref,
        ),
        %{$args_ref},
    },
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Methods: get_*
for my $key (
    # Returns a template for constructing the browser URL
    'browser_url_template',
    # Returns a template for constructing the revision part in the browser URL
    'browser_rev_template',
    # Returns a list of entries implied this entry
    'implied_entry_list',
    # Returns the component pattern for a location matching this entry
    'location_component_pattern',
    # Returns the entries for revision keywords
    'revision_entries',
) {
    no strict qw{refs};
    my $getter = "get_$key";
    *$getter = sub {
        my ($self) = @_;
        return $self->{$key};
    }
}

################################################################################
# Returns true if this is an implied entry
sub is_implied {
    my ($self) = @_;
    return $self->{is_implied};
}

1;
__END__

=head1 NAME

Fcm::Keyword::Entry::Location

=head1 SYNOPSIS

    use Fcm::Keyword::Entry::Location;

    $entry = Fcm::Keyword::Entry::Location->new({
        key => $key, value => $value, # ...
    });

    $key = $entry->get_key();
    $value = $entry->get_value();
    $revision_entries = $entry->get_revision_entries();

=head1 DESCRIPTION

This is a sub-class of L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>. An object of
this class represents a FCM location keyword entry.

=head1 METHODS

See L<Fcm::Keyword::Entry|Fcm::Keyword::Entry> for inherited methods.

=over 4

=item new($args_ref)

Constructor.

=item get_browser_url_template()

Returns the template string for constructing the browser URL. The string {1},
{2}, {3}, etc in the template string will be substituted by the components
captured by the location component pattern and the revision template. See
C<get_url_component_pattern()> and C<get_browser_rev_template()>.

=item get_browser_rev_template()

Returns the template string for constructing the revision part of the browser
URL. The string {1} in the template string will be substituted by the revision.
See C<get_browser_url_template()>.

=item get_implied_entry_list()

Returns a list of entries implied by this entry.

=item get_location_component_pattern()

Returns a regular expression, when matched against the scheme-specific-part in
the actual URI of a location in the namespace of this keyword entry, will
capture a list of components, which can then be used to replace the numbered
fields in the browser URL template. See C<get_browser_url_template()>.

=item get_revision_entries()

Returns a L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object containing the
revision keyword entries of this location.

=item is_implied()

Returns true if this is an implied entry.

=back

=head1 TO DO

Introduce a Fcm::Keyword::Config module to store entries constructor setting.

=head1 SEE ALSO

L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Config|Fcm::Keyword::Config>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
