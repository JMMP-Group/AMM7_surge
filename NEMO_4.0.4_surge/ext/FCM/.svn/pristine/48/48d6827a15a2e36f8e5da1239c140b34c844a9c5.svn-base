# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword;

use Carp qw{croak};
use Fcm::Config;
use Fcm::Exception;
use Fcm::Keyword::Config;
use Fcm::Keyword::Exception;
use URI;

my $ENTRIES;

my $PREFIX_OF_LOCATION_KEYWORD = 'fcm';
my $PATTERN_OF_RESERVED_REVISION_KEYWORDS
    = qr{\A (?:\d+|HEAD|BASE|COMMITTED|PREV|\{[^\}]+\}) \z}ixms;

################################################################################
# Returns the Fcm::Keyword::Entries object for storing the location entries
sub get_entries {
    my ($reset) = @_;
    if ($reset || !$ENTRIES) {
        $ENTRIES = Fcm::Keyword::Config::get_entries('LOCATION_ENTRIES');
    }
    return $ENTRIES;
}

################################################################################
# Returns a list of Fcm::Keyword::Entry::Location objects matching $in_loc
sub get_location_entries_for {
    my ($in_loc) = @_;
    my ($loc, $rev, @entry_trail_refs) = _parse_loc($in_loc);
    return (map {$_->[0]} @entry_trail_refs);
}

################################################################################
# Returns the prefix of location keyword (with or without the delimiter).
sub get_prefix_of_location_keyword {
    my ($with_delimiter) = @_;
    return $PREFIX_OF_LOCATION_KEYWORD . ($with_delimiter ? ':' : '');
}

################################################################################
# Expands (the keywords in) the specfied location (and REV), and returns them
sub expand {
    my ($in_loc, $in_rev) = @_;
    my ($loc, $rev) = _expand($in_loc, $in_rev);
    return _unparse_loc($loc, $rev, $in_rev);
}

################################################################################
# Returns the corresponding browser URL for the input VC location
sub get_browser_url {
    my ($in_loc, $in_rev) = @_;

    my ($loc, $rev, @entry_trail_refs) = _expand($in_loc, $in_rev);
    if (!@entry_trail_refs) {
        croak(Fcm::Keyword::Exception->new({message => sprintf(
            "%s: cannot be mapped to a browser URL", $in_loc,
        )}));
    }

    my @entries = map {$_->[0]} @entry_trail_refs;
    my $location_component_pattern
        = _get_browser_url_setting(\@entries, 'location_component_pattern');
    my $browser_url_template
        = _get_browser_url_setting(\@entries, 'browser_url_template');
    my $browser_rev_template
        = _get_browser_url_setting(\@entries, 'browser_rev_template');

    if (
           $location_component_pattern
        && $browser_url_template
        && $browser_rev_template
    ) {
        my $uri = URI->new($loc);
        my $sps = $uri->opaque();
        my @matches = $sps =~ $location_component_pattern;
        if (@matches) {
            my $result = $browser_url_template;
            for my $field_number (1 .. @matches) {
                my $match = $matches[$field_number - 1];
                $result =~ s/\{ $field_number \}/$match/xms;
            }
            my $rev_field = scalar(@matches) + 1;
            if ($rev) {
                my $rev_string = $browser_rev_template;
                $rev_string =~ s/\{1\}/$rev/xms;
                $result =~ s/\{ $rev_field \}/$rev_string/xms;
            }
            else {
                $result =~ s/\{ $rev_field \}//xms;
            }
            return $result;
        }
    }
    else {
        croak(Fcm::Keyword::Exception->new({message => sprintf(
            "%s: mapping templates not defined correctly", $in_loc,
        )}));
    }
}

################################################################################
# Returns a browser URL setting, helper function for get_browser_url()
sub _get_browser_url_setting {
    my ($entries_ref, $setting) = @_;
    my $getter = "get_$setting";
    for my $entry (@{$entries_ref}) {
        my $setting = $entry->$getter();
        if ($setting) {
            return $setting;
        }
    }
    my $config = Fcm::Config->instance();
    return $config->setting('URL_BROWSER_MAPPING_DEFAULT', uc($setting));
}

################################################################################
# Un-expands the specfied location (and REV) to keywords, and returns them
sub unexpand {
    my ($in_loc, $in_rev) = @_;
    my ($loc, $rev, @entry_trail_refs) = _parse_loc($in_loc, $in_rev);
    if (@entry_trail_refs) {
        my ($entry, $trail) = @{$entry_trail_refs[0]};
        if ($rev) {
            GET_REV_KEY:
            for my $entry_trail_ref (@entry_trail_refs) {
                my ($e, $t) = @{$entry_trail_ref};
                my $rev_key
                    = $e->get_revision_entries()->get_entry_by_value($rev);
                if ($rev_key) {
                    $rev = $rev_key->get_key();
                    last GET_REV_KEY;
                }
            }
        }
        $loc = get_prefix_of_location_keyword(1) . $entry->get_key() . $trail;
        return _unparse_loc($loc, $rev, $in_rev);
    }
    return _unparse_loc($in_loc, $in_rev, $in_rev);
}

################################################################################
# Expands (the keywords in) the specfied location (and REV), and returns them
sub _expand {
    my ($in_loc, $in_rev) = @_;
    my ($loc, $rev, @entry_trail_refs) = _parse_loc($in_loc, $in_rev);
    if (@entry_trail_refs) {
        my ($entry, $trail) = @{$entry_trail_refs[0]};
        $loc = $entry->get_value() . $trail;
        if ($rev && $rev !~ $PATTERN_OF_RESERVED_REVISION_KEYWORDS) {
            my $r;
            GET_REV:
            for my $entry_trail_ref (@entry_trail_refs) {
                my ($e, $t) = @{$entry_trail_ref};
                $r = $e->get_revision_entries()->get_entry_by_key($rev);
                if ($r) {
                    $rev = $r->get_value();
                    last GET_REV;
                }
            }
            if (!$r) {
                croak(Fcm::Keyword::Exception->new({message => sprintf(
                    "%s: %s: unknown revision keyword",
                    $loc, $rev,
                )}));
            }
        }
    }
    return ($loc, $rev, @entry_trail_refs);
}

################################################################################
# Parses $in_loc (and $in_rev)
sub _parse_loc {
    my ($in_loc, $in_rev) = @_;
    if (!$in_loc) {
        croak(Fcm::Exception->new({
            message => 'internal error: $in_loc not defined',
        }));
    }
    if ($in_loc) {
        if (!defined($in_rev)) {
            my ($loc, $rev) = $in_loc =~ qr{\A (.+) \@ ([^/\@]+) \z}xms;
            if ($loc && $rev) {
                return ($loc, $rev, _get_loc_entry($loc));
            }
            else {
                return ($in_loc, $in_rev, _get_loc_entry($in_loc));
            }
        }
        return ($in_loc, $in_rev, _get_loc_entry($in_loc));
    }
    return;
}

################################################################################
# Returns a list of keyword entries/trailing path pairs for the input location
sub _get_loc_entry {
    my ($loc) = @_;
    if ($loc) {
        my $uri = URI->new($loc);
        if (
               $uri->scheme()
            && $uri->scheme() eq get_prefix_of_location_keyword()
        ) {
            my ($key, $trail) = $uri->opaque() =~ qr{\A ([^/\@]+) (.*) \z}xms;
            my $entry = get_entries()->get_entry_by_key($key);
            if (!$entry || !$entry->get_value()) {
                die(Fcm::Keyword::Exception->new({message => sprintf(
                    "%s: unknown FCM location keyword", $loc,
                )}));
            }
            $loc = $entry->get_value() . ($trail ? $trail : q{});
        }
        my @entry_trail_pairs = ();
        my $lead = $loc;
        GET_ENTRY:
        while ($lead) {
            my $entry = get_entries()->get_entry_by_value($lead);
            if ($entry) {
                my $trail = substr($loc, length($lead));
                push @entry_trail_pairs, [$entry, $trail];
            }
            if (!($lead =~ s{/+ [^/]* \z}{}xms)) {
                last GET_ENTRY;
            }
        }
        if (@entry_trail_pairs) {
            return @entry_trail_pairs;
        }
        else {
            return;
        }
    }
    return;
}

################################################################################
# If $in_rev, returns (LOC, REV). Otherwise, returns LOC@REV
sub _unparse_loc {
    my ($loc, $rev, $in_rev) = @_;
    if (!$loc) {
        return;
    }
    return ($in_rev ? ($loc, $rev) : join(q{@}, $loc, ($rev ? $rev : ())));
}

1;
__END__

=head1 NAME

Fcm::Keyword

=head1 SYNOPSIS

    use Fcm::Keyword;

    $loc = Fcm::Keyword::expand('fcm:namespace/path@rev-keyword');
    $loc = Fcm::Keyword::unexpand('svn://host/namespace/path@1234');

    ($loc, $rev) = Fcm::Keyword::expand('fcm:namespace/path', 'rev-keyword');
    ($loc, $rev) = Fcm::Keyword::unexpand('svn://host/namespace/path', 1234);

    $loc = Fcm::Keyword::get_browser_url('fcm:namespace/path');
    $loc = Fcm::Keyword::get_browser_url('svn://host/namespace/path');

    $loc = Fcm::Keyword::get_browser_url('fcm:namespace/path@1234');
    $loc = Fcm::Keyword::get_browser_url('svn://host/namespace/path@1234');

    $loc = Fcm::Keyword::get_browser_url('fcm:namespace/path', 1234);
    $loc = Fcm::Keyword::get_browser_url('svn://host/namespace/path', 1234);

    $entries = Fcm::Keyword::get_entries();

=head1 DESCRIPTION

This module contains utilities to expand and unexpand FCM location and revision
keywords.

=head1 FUNCTIONS

=over 4

=item expand($loc)

Expands FCM keywords in $loc and returns the result.

If $loc is a I<fcm> scheme URI, the leading part (before any "/" or "@"
characters) of the URI opaque is the namespace of a FCM location keyword. This
is expanded into the actual value. Optionally, $loc can be suffixed with a peg
revision (an "@" followed by any characters). If a peg revision is a FCM
revision keyword, it is expanded into the actual revision.

=item expand($loc,$rev)

Same as C<expand($loc)>, but $loc should not contain a peg revision. Returns a
list containing the expanded version of $loc and $rev.

=item get_browser_url($loc)

Given a repository $loc in a known keyword namespace, returns the corresponding
URL for the code browser.

Optionally, $loc can be suffixed with a peg revision (an "@" followed by any
characters).

=item get_browser_url($loc,$rev)

Same as get_browser_url($loc), but the revision should be specified using $rev
but not pegged with $loc.

=item get_entries([$reset])

Returns the L<Fcm::Keyword::Entries|Fcm::Keyword::Entries> object for storing
location keyword entries. If $reset if true, reloads the entries.

=item get_location_entries_for($loc)

Returns a list of L<Fcm::Keyword::Entry::Location|Fcm::Keyword::Entry::Location>
objects matching $loc.

=item get_prefix_of_location_keyword($with_delimiter)

Returns the prefix of a FCM location keyword, (currently "fcm"). If
$with_delimiter is specified and is true, returns the prefix with the delimiter,
(currently "fcm:").

=item unexpand($loc)

Does the opposite of expand($loc). Returns the FCM location keyword equivalence
of $loc. If the $loc can be mapped using 2 or more namespaces, the namespace
that results in the longest substitution is used. Optionally, $loc can be
suffixed with a peg revision (an "@" followed by any characters). If a peg
revision is a known revision, it is turned into its corresponding revision
keyword.

=item unexpand($loc,$rev)

Same as unexpand($loc), but $loc should not contain a peg revision. Returns a
list containing the unexpanded version of $loc and $rev 

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::Keyword::Exception|Fcm::Keyword::Exception>

Functions in this module may die() with this exception when it fails to expand
a keyword.

=back

=head1 SEE ALSO

L<Fcm::Keyword::Config|Fcm::Keyword::Config>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Entry::Location|Fcm::Keyword::Entry::Location>,
L<Fcm::Keyword::Exception|Fcm::Keyword::Exception>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
