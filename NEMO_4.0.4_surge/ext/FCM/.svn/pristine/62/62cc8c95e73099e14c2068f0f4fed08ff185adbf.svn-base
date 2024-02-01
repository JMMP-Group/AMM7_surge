# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Loader::Config::Revision;

use Fcm::Config;

sub new {
    my ($class, $args_ref) = @_;
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Returns the namespace where the revision keywords belong
sub get_namespace {
    my ($self) = @_;
    return $self->{namespace};
}

################################################################################
# Returns 'Fcm::Config'
sub get_source {
    my ($self) = @_;
    return 'Fcm::Config';
}

################################################################################
# Loads revision keywords from Fcm::Config to $entries
sub load_to {
    my ($self, $entries) = @_;
    my $load_counter = 0;
    my $config = $self->get_source()->instance();
    my $rev_keyword_ref = $config->setting(
        qw/URL_REVISION/,
        uc($self->get_namespace()),
    );
    if ($rev_keyword_ref) {
        for my $key (keys(%{$rev_keyword_ref})) {
            $entries->add_entry($key, $rev_keyword_ref->{$key});
            $load_counter++;
        }
    }
    return ($config->is_initialising() ? 0 : defined($load_counter));
}

1;
__END__

=head1 NAME

Fcm::Keyword::Loader::Config::Revision

=head1 SYNOPSIS

    $loader = Fcm::Keyword::Loader::Config::Revision->new({namespace => $name});
    $loader->load_to($entries);

=head1 DESCRIPTION

Loads revision keywords from L<Fcm::Config|Fcm::Config> into a
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object containing
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry> objects.

=head1 METHODS

=over 4

=item C<new({namespace =E<gt> $namespace})>

Constructor. The argument $namespace is the namespace where the revision
keywords belong.

=item get_namespace()

Returns the namespace where the revision keywords belong.

=item get_source()

Returns the string "L<Fcm::Config|Fcm::Config>".

=item load_to($entries)

Loads revision keywords in the namespace given by C<$self-E<gt>get_namespace()>
from L<Fcm::Config|Fcm::Config> to $entries. Returns true on success. (However,
if L<Fcm::Config|Fcm::Config> is initialising, returns false to force a reload
next time.)

=back

=head1 SEE ALSO

L<Fcm::Config|Fcm::Config>,
L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Entry::Location|Fcm::Keyword::Entry::Location>,
L<Fcm::Keyword::Loader|Fcm::Keyword::Loader>,
L<Fcm::Keyword::Loader::Config::Location|Fcm::Keyword::Loader::Config::Location>
L<Fcm::Keyword::Loader::VC::Revision|Fcm::Keyword::Loader::VC::Revision>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
