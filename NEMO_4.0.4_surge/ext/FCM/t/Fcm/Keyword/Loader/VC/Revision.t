#!/usr/bin/perl

use strict;
use warnings;

use Carp qw{croak};
use Fcm::Keyword::Entries;
use File::Basename qw{dirname};
use File::Spec;
use File::Temp qw{tempdir};
use IO::File;
use IO::Pipe;
use POSIX qw{WIFEXITED};
use Test::More (tests => 17);

my %VALUE_OF = (
    bar => {
        'bar3'    => 3,
        'bar3.1'  => 31,
        'bar3.14' => 314,
    },
    baz => {
        'bear'   => 4,
        'bee'    => 6,
        'spider' => 8,
    },
);

main();

sub main {
    my $class = 'Fcm::Keyword::Loader::VC::Revision';
    use_ok($class);
    test_constructor($class);
    test_load_to($class);
}

################################################################################
# Tests simple usage of the constructor
sub test_constructor {
    my ($class) = @_;
    my $prefix = "constructor";
    my $loader = $class->new({source => 'foo'});
    isa_ok($loader, $class);
    is($loader->get_source(), 'foo', "$prefix: get_source()");
    ok($loader->load_to(), "$prefix: load_to"); # FIXME: should fail?
}

################################################################################
# Tests loading to an Fcm::Keyword::Entries object
sub test_load_to {
    my ($class) = @_;
    my $prefix = 'load to';
    my $temp_dir = tempdir(CLEANUP => 1);
    my $repos = File::Spec->catfile($temp_dir, 'repos');
    WIFEXITED(system(qw{svnadmin create}, $repos))
        || croak("$repos: cannot create: $?");
    my $dump_file = File::Spec->catfile(dirname($0), 'Revision.dump');
    my $handle = IO::File->new($dump_file, 'r');
    if (!$handle) {
        croak("$dump_file: cannot load: $!");
    }
    my $dump = do{local $/; $handle->getline()};
    $handle->close();
    my $pipe = IO::Pipe->new();
    $pipe->writer(qw{svnadmin load -q}, $repos);
    print($pipe $dump);
    $pipe->close();
    if ($?) {
        croak("$dump_file: cannot load: $?");
    }
    my $repos_url = "file://$repos";
    my $loader = $class->new({source => $repos_url});
    my $entries = Fcm::Keyword::Entries->new();
    ok($loader->load_to($entries), "$prefix: nothing to load");
    for my $key (keys(%VALUE_OF)) {
        my $url = "$repos_url/$key";
        my $loader = $class->new({source => $url});
        $loader->load_to($entries);
        for my $rev_key (keys(%{$VALUE_OF{$key}})) {
            my $entry = $entries->get_entry_by_key($rev_key);
            if ($entry) {
                is(
                    $entry->get_key(),
                    uc($rev_key),
                    "$prefix: by key: $rev_key",
                );
                is(
                    $entries->get_entry_by_value($VALUE_OF{$key}{$rev_key}),
                    $entry,
                    "$prefix: by value: $rev_key: object",
                );
            }
            else {
                fail("$prefix: by key: $rev_key");
            }
        }
    }
}

__END__
