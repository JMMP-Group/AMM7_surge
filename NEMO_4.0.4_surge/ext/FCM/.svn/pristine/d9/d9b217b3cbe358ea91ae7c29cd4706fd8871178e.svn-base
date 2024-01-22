# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Loader::Config::Location;

use Fcm::Config;

my %IMPLIED_NAMESPACE_SUFFIX = (tr => 'trunk', br => 'branches', tg => 'tags');

sub new {
    my ($class) = @_;
    return bless(\do{my $annon_scalar}, $class);
}

################################################################################
# Returns 'Fcm::Config'
sub get_source {
    my ($self) = @_;
    return 'Fcm::Config';
}

################################################################################
# Loads location keywords from Fcm::Config to $entries
sub load_to {
    my ($self, $entries) = @_;
    my $config = $self->get_source()->instance();
    my $load_counter = 0;
    for my $key (keys(%{$config->setting('URL')})) {
        my $value = $config->setting('URL', $key);
        my $location_component_pattern = $config->setting(
            'URL_BROWSER_MAPPING', $key, 'LOCATION_COMPONENT_PATTERN');
        my $browser_url_template = $config->setting(
            'URL_BROWSER_MAPPING', $key, 'BROWSER_URL_TEMPLATE');
        my $browser_rev_template = $config->setting(
            'URL_BROWSER_MAPPING', $key, 'BROWSER_REV_TEMPLATE');
        my $entry = $entries->add_entry(
            $key,
            $value,
            {
                location_component_pattern => $location_component_pattern,
                browser_url_template       => $browser_url_template,
                browser_rev_template       => $browser_rev_template,
            },
        );
        $load_counter++;

        # Set up implied keywords
        for my $suffix (keys(%IMPLIED_NAMESPACE_SUFFIX)) {
            my $value_suf = $value . '/' . $IMPLIED_NAMESPACE_SUFFIX{$suffix};
            for my $join (q{_}, q{-}) {
                my $implied_entry = $entries->add_entry(
                    uc($key . $join . $suffix),
                    $value_suf,
                    {is_implied => 1},
                );
                push(@{$entry->get_implied_entry_list()}, $implied_entry);
                $load_counter++;
            }
        }
    }
    return ($config->is_initialising() ? 0 : defined($load_counter));
}

1;
__END__

=head1 NAME

Fcm::Keyword::Loader::Config::Location

=head1 SYNOPSIS

    $loader = Fcm::Keyword::Loader::Config::Location->new();
    $loader->load_to($entries);

=head1 DESCRIPTION

This class implements the L<Fcm::Keyword::Loader|Fcm::Keyword::Loader>
interface.

Loads location keywords from L<Fcm::Config|Fcm::Config> into a
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object containing
L<Fcm::Keyword::Entry::Location|Fcm::Keyword::Entry::Location> objects.

=head1 METHODS

=over 4

=item new()

Constructor.

=item get_source()

Returns the string "L<Fcm::Config|Fcm::Config>".

=item load_to($entries)

Loads location keywords and implied keywords from L<Fcm::Config|Fcm::Config> to
$entries. It also loads settings for mapping location to browser URL. Returns
true on success. (However, if L<Fcm::Config|Fcm::Config> is initialising,
returns false to force a reload next time.)

=back

=head1 TO DO

Need a more flexible system for implied keywords.

=head1 SEE ALSO

L<Fcm::Config|Fcm::Config>,
L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Loader|Fcm::Keyword::Loader>,
L<Fcm::Keyword::Loader::Config::Revision|Fcm::Keyword::Loader::Config::Revision>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
