# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Exception;
use base qw{Fcm::Exception};

1;
__END__

=head1 NAME

Fcm::Keyword::Exception

=head1 SYNOPSIS

    use Carp qw{croak};
    use Fcm::Keyword::Exception;
    croak(Fcm::Keyword::Exception->new({message => 'something is wrong'}));

=head1 DESCRIPTION

This class extends L<Fcm::Exception|Fcm::Exception>. This exception is thrown
on errors associated with the command line interface.

=head1 METHODS

See L<Fcm::Exception|Fcm::Exception> for a list of methods.

=head1 SEE ALSO

L<Fcm::Exception|Fcm::Exception>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
