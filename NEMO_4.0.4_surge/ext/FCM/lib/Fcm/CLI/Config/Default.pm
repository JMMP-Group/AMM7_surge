# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Config::Default;

use Fcm::CLI::Option;
use Fcm::CLI::Subcommand;

my %DESCRIPTION_OF = (
    # --------------------------------------------------------------------------
    BROWSER => <<'END_DESCRIPTION',
If TARGET is specified, it must be a FCM URL keyword, a Subversion URL or the
path to a local working copy. If not specified, the current working directory
is assumed to be a working copy. If the --browser option is specified, the
specified web browser command is used to launch the repository browser.
Otherwise, it attempts to use the default browser from the configuration
setting.
END_DESCRIPTION
    # --------------------------------------------------------------------------
    BUILD => <<'END_DESCRIPTION',
The path to a CFGFILE may be provided. Otherwise, the build system searches the
default locations for a bld cfg file.

If no option is specified, the options "-s 5 -t all -j 1 -v 1" are assumed.

If the option for full build is specified, the sub-directories created by
previous builds will be removed, so that the current build can start cleanly.

The -s option can be used to limit the actions performed by the build system up
to a named stage. The stages are:
    "1", "s" or "setup"                - stage 1, setup
    "2", "pp" or "pre_process"         - stage 2, pre-process
    "3", "gd" or "generate_dependency" - stage 3, generate dependency
    "4", "gi" or "generate_interface"  - stage 4, generate Fortran 9X interface
    "5", "m", "make"                   - stage 5, make

If a colon separated list of targets is specified using the -t option, the
default targets specified in the configuration file will not be used.

If archive mode is switched on, build sub-directories that are only used in the
build process will be archived to TAR files. The default is off.

If specified, the verbose level must be an integer greater than 0. Verbose
level 0 is the quiet mode. Increasing the verbose level will increase the
amount of diagnostic output.

When a build is invoked, it sets up a lock file in the build root directory.
The lock is normally removed at the end of the build. While the lock file is in
place, the build commands invoked in the same root directory will fail.  If
you need to bypass this check for whatever reason, you can invoke the build
system with the --ignore-lock option.
END_DESCRIPTION
    # --------------------------------------------------------------------------
    CFG_PRINTER => <<'END_DESCRIPTION',
If no option is specified, the output will be sent to standard output.
END_DESCRIPTION
    # --------------------------------------------------------------------------
    EXTRACT => <<'END_DESCRIPTION',
The path to a CFG file may be provided. Otherwise, the extract system searches
the default locations for an ext cfg file.

If no option is specified, the system will attempt an incremental extract where
appropriate.

If specified, the verbose level must be an integer greater than 0. Verbose
level 0 is the quiet mode. Increasing the verbose level will increase the
amount of diagnostic output.

When an extract is invoked, it sets up a lock file in the extract destination
root directory. The lock is normally removed at the end of the extract. While
the lock file is in place, other extract commands invoked in the same
destination root directory will fail. If you need to bypass this check for
whatever reason, you can invoke the extract system with the --ignore-lock
option.
END_DESCRIPTION
    # --------------------------------------------------------------------------
    EXTRACT_CONFIG_COMPARATOR => <<'END_DESCRIPTION',
Compares the extract configurations of two similar extract configuration files
CFGFILE1 and CFGFILE2.

In normal mode with verbosity level 2 or above, displays the change log of each
revision.

In wiki mode, print revision tables in wiki format. The argument to the --wiki
option must be the Subversion URL or FCM URL keyword of a FCM project
associated with the intended Trac system. The --verbose option has no effect
in wiki mode.
END_DESCRIPTION
    # --------------------------------------------------------------------------
    GUI => <<'END_DESCRIPTION',
The optional argument PATH modifies the initial working directory of the GUI.
END_DESCRIPTION
    # --------------------------------------------------------------------------
    KEYWORD => <<'END_DESCRIPTION',
If no argument is specified, prints registered location keywords. Otherwise,
prints the implied location keywords and revision keywords for the specified
target.
END_DESCRIPTION
);

my %OPTION_OF = (
    ARCHIVE => Fcm::CLI::Option->new({
        name        => 'archive',
        letter      => 'a',
        description => 'archives sub-directories on success',
    }),

    BROWSER => Fcm::CLI::Option->new({
        name        => 'browser',
        letter      => 'b',
        description => 'specifies the web browser command',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
    }),

    CLEAN => Fcm::CLI::Option->new({
        name        => 'clean',
        description => 'cleans the destination',
    }),

    FULL => Fcm::CLI::Option->new({
        name        => 'full',
        letter      => 'f',
        description => 'runs in full mode',
    }),

    HELP => Fcm::CLI::Option->new({
        name        => 'help',
        letter      => 'h',
        description => 'prints help',
        is_help     => 1,
    }),

    IGNORE_LOCK => Fcm::CLI::Option->new({
        name        => 'ignore-lock',
        description => 'ignores lock file',
    }),

    JOBS => Fcm::CLI::Option->new({
        name        => 'jobs',
        letter      => 'j',
        description => 'number of parallel jobs',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
    }),

    OUTPUT => Fcm::CLI::Option->new({
        name        => 'output',
        letter      => 'o',
        description => 'sends output to the specified file',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
    }),

    STAGE => Fcm::CLI::Option->new({
        name        => 'stage',
        letter      => 's',
        description => 'runs command up to a named stage',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
    }),

    TARGETS => Fcm::CLI::Option->new({
        name        => 'targets',
        letter      => 't',
        delimiter   => ':',
        description => 'list of build targets, delimited by (:)',
        has_arg     => Fcm::CLI::Option->ARRAY_ARG,
    }),

    VERBOSITY => Fcm::CLI::Option->new({
        name        => 'verbose',
        letter      => 'v',
        description => 'verbose level',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
    }),

    WIKI => Fcm::CLI::Option->new({
        name        => 'wiki',
        letter      => 'w',
        description => 'print revision tables in wiki format',
        has_arg     => Fcm::CLI::Option->SCALAR_ARG,
    }),
);

my %SUBCOMMAND_OF = (
    BRANCH => Fcm::CLI::Subcommand->new({
        names         => ['branch', 'br'],
        synopsis      => 'branch utilities',
        invoker_class => 'Fcm::CLI::Invoker::CM',
        is_vc         => 1,
    }),

    BROWSER => Fcm::CLI::Subcommand->new({
        names         => ['trac', 'www'],
        synopsis      => 'invokes the browser for a version controlled target',
        usage         => '[OPTIONS...] [TARGET]',
        description   => $DESCRIPTION_OF{BROWSER},
        invoker_class => 'Fcm::CLI::Invoker::Browser',
        options       => [
            $OPTION_OF{BROWSER},
            $OPTION_OF{HELP},
        ],
    }),

    BUILD => Fcm::CLI::Subcommand->new({
        names          => ['build', 'bld'],
        synopsis       => 'invokes the build system',
        usage          => '[OPTIONS...] [CFGFILE]',
        description    => $DESCRIPTION_OF{BUILD},
        invoker_class  => 'Fcm::CLI::Invoker::ConfigSystem',
        invoker_config => {
            impl_class         => 'Fcm::Build',
            cli2invoke_key_map => {
                'archive'     => 'ARCHIVE',
                'clean'       => 'CLEAN',
                'full'        => 'FULL',
                'ignore-lock' => 'IGNORE_LOCK',
                'jobs'        => 'JOBS',
                'stage'       => 'STAGE',
                'targets'     => 'TARGETS',
            },
        },
        options        => [
            $OPTION_OF{ARCHIVE},
            $OPTION_OF{CLEAN},
            $OPTION_OF{FULL},
            $OPTION_OF{HELP},
            $OPTION_OF{IGNORE_LOCK},
            $OPTION_OF{JOBS},
            $OPTION_OF{STAGE},
            $OPTION_OF{TARGETS},
            $OPTION_OF{VERBOSITY},
        ],
    }),

    CFG_PRINTER => Fcm::CLI::Subcommand->new({
        names         => ['cfg'],
        synopsis      => 'invokes the CFG file pretty printer',
        usage         => '[OPTIONS...] [CFGFILE]',
        description   => $DESCRIPTION_OF{CFG_PRINTER},
        invoker_class => 'Fcm::CLI::Invoker::CfgPrinter',
        options       => [
            $OPTION_OF{HELP},
            $OPTION_OF{OUTPUT},
        ],
    }),

    CM => Fcm::CLI::Subcommand->new({
        names => [qw{
            add
            blame     praise annotate ann
            cat
            checkout  co
            cleanup
            commit    ci
            copy      cp
            delete    del    remove   rm
            diff      di
            export
            import
            info
            list      ls
            lock
            log
            merge
            mkdir
            move      mv     rename   ren
            propdel   pdel   pd
            propedit  pedit  pe
            propget   pget   pg
            proplist  plist  pl
            propset   pset   ps
            resolved
            revert
            status    stat   st
            switch    sw
            unlock
            update    up
        }],
        invoker_class => 'Fcm::CLI::Invoker::CM',
        is_vc         => 1,
    }),

    CONFLICTS => Fcm::CLI::Subcommand->new({
        names         => ['conflicts', 'cf'],
        synopsis      => 'resolves conflicts in your working copy',
        usage         => '[PATH]',
        invoker_class => 'Fcm::CLI::Invoker::CM',
        is_vc         => 1,
    }),

    EXTRACT => Fcm::CLI::Subcommand->new({
        names          => ['extract', 'ext'],
        synopsis       => 'invokes the extract system',
        usage          => '[OPTIONS...] [CFGFILE]',
        description    => $DESCRIPTION_OF{EXTRACT},
        invoker_class  => 'Fcm::CLI::Invoker::ConfigSystem',
        invoker_config => {
            impl_class         => 'Fcm::Extract',
            cli2invoke_key_map => {
                'clean'       => 'CLEAN',
                'full'        => 'FULL',
                'ignore-lock' => 'IGNORE_LOCK',
            },
        },
        options        => [
            $OPTION_OF{CLEAN},
            $OPTION_OF{FULL},
            $OPTION_OF{HELP},
            $OPTION_OF{IGNORE_LOCK},
            $OPTION_OF{VERBOSITY},
        ],
    }),

    EXTRACT_CONFIG_COMPARATOR => Fcm::CLI::Subcommand->new({
        names         => ['cmp-ext-cfg'],
        synopsis      => 'invokes the extract configuration files comparator',
        usage         => '[OPTIONS...] CFGFILE1 CFGFILE2',
        description   => $DESCRIPTION_OF{EXTRACT_CONFIG_COMPARATOR},
        invoker_class => 'Fcm::CLI::Invoker::ExtractConfigComparator',
        options       => [
            $OPTION_OF{HELP},
            $OPTION_OF{VERBOSITY},
            $OPTION_OF{WIKI},
        ],
    }),

    GUI => Fcm::CLI::Subcommand->new({
        names         => ['gui'],
        synopsis      => 'invokes the GUI wrapper for code management commands',
        usage         => '[PATH]',
        description   => $DESCRIPTION_OF{GUI},
        invoker_class => 'Fcm::CLI::Invoker::GUI',
    }),

    HELP => Fcm::CLI::Subcommand->new({
        names         => ['help', q{?}, q{}],
        synopsis      => 'displays the usage of this program or its subcommands',
        usage         => '[SUBCOMMAND]',
        description   => q{},
        invoker_class => 'Fcm::CLI::Invoker::Help',
        options       => [$OPTION_OF{HELP}],
    }),

    KEYWORD => Fcm::CLI::Subcommand->new({
        names         => ['keyword-print', 'kp'],
        synopsis      => 'prints registered location and/or revision keywords',
        usage         => '[TARGET]',
        description   => $DESCRIPTION_OF{KEYWORD},
        invoker_class => 'Fcm::CLI::Invoker::KeywordPrinter',
        options       => [$OPTION_OF{HELP}],
    }),

    MKPATCH => Fcm::CLI::Subcommand->new({
        names         => ['mkpatch'],
        synopsis      => 'creates patches from specified revisions of a URL',
        usage         => '[OPTIONS] URL [OUTDIR]',
        invoker_class => 'Fcm::CLI::Invoker::CM',
        is_vc         => 1,
    }),
);

our @CORE_SUBCOMMANDS = (
    $SUBCOMMAND_OF{HELP},
    $SUBCOMMAND_OF{BUILD},
    $SUBCOMMAND_OF{CFG_PRINTER},
);

our @VC_SUBCOMMANDS = (
    $SUBCOMMAND_OF{BRANCH},
    $SUBCOMMAND_OF{BROWSER},
    $SUBCOMMAND_OF{CONFLICTS},
    $SUBCOMMAND_OF{EXTRACT},
    $SUBCOMMAND_OF{EXTRACT_CONFIG_COMPARATOR},
    $SUBCOMMAND_OF{GUI},
    $SUBCOMMAND_OF{KEYWORD},
    $SUBCOMMAND_OF{MKPATCH},
    $SUBCOMMAND_OF{CM},
);

1;
__END__

=head1 NAME

Fcm::CLI::Config::Default

=head1 SYNOPSIS

    use Fcm::CLI::Config::Default;
    @core_subcommands = @Fcm::CLI::Config::Default::CORE_SUBCOMMANDS;
    @vc_subcommands = @Fcm::CLI::Config::Default::VC_SUBCOMMANDS;

=head1 DESCRIPTION

This module stores the default configuration of the FCM command line interface.
It should only be used by L<Fcm::CLI::Config|Fcm::CLI::Config>.

=head1 SEE ALSO

L<Fcm::CLI::Config|Fcm::CLI::Config>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>,
L<Fcm::CLI::Option|Fcm::CLI::Option>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
