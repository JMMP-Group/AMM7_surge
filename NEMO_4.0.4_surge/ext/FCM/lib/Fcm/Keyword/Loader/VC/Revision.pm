# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Loader::VC::Revision;

use Fcm::Util qw{run_command};

sub new {
    my ($class, $args_ref) = @_;
    return bless({%{$args_ref}}, $class);
}

################################################################################
# Returns the VC location where revision keywords will be loaded from
sub get_source {
    my ($self) = @_;
    return $self->{source};
}

################################################################################
# Loads revision keywords from $self->get_source() to $entries
sub load_to {
    my ($self, $entries) = @_;
    my @lines = run_command(
        [qw{svn pg fcm:revision}, $self->get_source()],
        DEVNULL => 1,
        ERROR   => 'ignore',
        METHOD  => 'qx',
    );
    my $load_counter = 0;
    for my $line (@lines) {
        chomp($line);
        my ($key, $value) = split(qr{\s+ = \s+}xms, $line);
        if ($key && $value) {
            $entries->add_entry($key, $value);
            $load_counter++;
        }
    }
    return defined($load_counter);
}

1;
__END__

=head1 NAME

Fcm::Keyword::Loader::VC::Revision

=head1 SYNOPSIS

    $loader = Fcm::Keyword::Loader::VC::Revision->new({source => $source});
    $loader->load_to($entries);

=head1 DESCRIPTION

Loads revision keywords from a VC location into a
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object containing
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry> objects.

=head1 METHODS

=over 4

=item C<new({source =E<gt> $source})>

Constructor. The argument $source is the VC location from which revision
keywords will be loaded from.

=item get_source()

Returns the source VC location from which revision keywords will be loaded
from.

=item load_to($entries)

Loads revision keywords from C<$self-E<gt>get_source()> to $entries.

=back

=head1 TO DO

Abstract away the call to the VC system, which assumes the Subversion shell
client at the moment.

=head1 SEE ALSO

L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Entry::Location|Fcm::Keyword::Entry::Location>,
L<Fcm::Keyword::Loader|Fcm::Keyword::Loader>,
L<Fcm::Keyword::Loader::Config::Revision|Fcm::Keyword::Loader::Config::Revision>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
