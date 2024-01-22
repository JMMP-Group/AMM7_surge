# ------------------------------------------------------------------------------
# NAME
#   Fcm::Cm
#
# DESCRIPTION
#   This module contains the FCM code management functionalities and wrappers
#   to Subversion commands.
#
# COPYRIGHT
#   (C) Crown copyright Met Office. All rights reserved.
#   For further details please refer to the file COPYRIGHT.txt
#   which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Cm;
use base qw{Exporter};

our @EXPORT_OK = qw(cli cm_check_missing cm_check_unknown cm_switch cm_update);

use Cwd            qw{cwd};
use Getopt::Long   qw{GetOptions :config bundling};
use Fcm::CLI::Exception;
use Fcm::Config;
use Fcm::CmBranch;
use Fcm::CmUrl;
use Fcm::Keyword;
use Fcm::Util      qw{
    get_url_of_wc
    get_url_peg_of_wc
    get_wct
    is_url
    is_wc
    run_command
    tidy_url
};
use File::Basename qw{basename dirname};
use File::Path     qw{mkpath rmtree};
use File::Spec;
use File::Temp     qw{tempfile};
use Pod::Usage     qw{pod2usage};

# ------------------------------------------------------------------------------

# CLI message handler
our $CLI_MESSAGE = \&_cli_message;

# List of CLI messages
our %CLI_MESSAGE_FOR = (
    q{}         => "%s",
    BRANCH_LIST => "%s at %s: %d branch(es) found for %s.\n",
    CHDIR_WCT   => "%s: working directory changed to top of working copy.\n",
    CF          => "Conflicts in: %s\n",
    MERGE       => "Performing merge ...\n",
    MERGE_CF    => "About to merge in changes from %s compared with %s\n",
    MERGE_CI    => "The following is added to the commit message file:\n%s",
    MERGE_DRY   => "This merge will result in the following change:\n",
    MERGE_REVS  => "Merge(s) available from %s: %s\n",
    OUT_DIR     => "Output directory: %s\n",
    PATCH_DONE  => "%s: patch generated.\n",
    PATCH_REV   => "Patch created for changeset %s\n",
    SEPARATOR   => q{-} x 80 . "\n",
    STATUS      => "Status of the target working copy(ies):\n%s",
);

# CLI abort and error messages
our %CLI_MESSAGE_FOR_ABORT = (
    FAIL => "%s: command failed.\n",
    NULL => "%s: command will result in no change.\n",
    USER => "%s: abort by user.\n",
);

# CLI abort and error messages
our %CLI_MESSAGE_FOR_ERROR = (
    CHDIR               => "%s: cannot change to directory.\n",
    CLI                 => "%s",
    CLI_HELP            => "Type 'fcm help %s' for usage.\n",
    CLI_MERGE_ARG1      => "Arg 1 must be the source in auto/custom mode.\n",
    CLI_MERGE_ARG2      => "Arg 2 must be the source in custom mode"
                           . " if --revision not set.\n",
    CLI_OPT_ARG         => "--%s: invalid argument [%s].\n",
    CLI_OPT_WITH_OPT    => "--%s: must be specified with --%s.\n",
    CLI_USAGE           => "incorrect usage",
    DIFF_PROJECTS       => "%s (target) and %s (source) are not related.\n",
    INVALID_BRANCH      => "%s: not a valid URL of a standard FCM branch.\n",
    INVALID_PROJECT     => "%s: not a valid URL of a standard FCM project.\n",
    INVALID_TARGET      => "%s: not a valid working copy or URL.\n",
    INVALID_URL         => "%s: not a valid URL.\n",
    INVALID_WC          => "%s: not a valid working copy.\n",
    MERGE_REV_INVALID   => "%s: not a revision in the available merge list.\n",
    MERGE_SELF          => "%s: cannot be merged to its own working copy: %s.\n",
    MERGE_UNRELATED     => "%s: target and %s: source not directly related.\n",
    MERGE_UNSAFE        => "%s: source contains changes outside the target"
                           . " sub-directory. Please merge with a full tree.\n",
    MKPATH              => "%s: cannot create directory.\n",
    NOT_EXIST           => "%s: does not exist.\n",
    PARENT_NOT_EXIST    => "%s: parent %s no longer exists.\n",
    RMTREE              => "%s: cannot remove.\n",
    ST_CONFLICT         => "File(s) in conflicts:\n%s",
    ST_MISSING          => "File(s) missing:\n%s",
    ST_OUT_OF_DATE      => "File(s) out of date:\n%s",
    SWITCH_UNSAFE       => "%s: merge template exists."
                           . " Please remove before retrying.\n",
    WC_EXIST            => "%s: working copy already exists.\n",
    WC_INVALID_BRANCH   => "%s: not a working copy of a standard FCM branch.\n",
    WC_URL_NOT_EXIST    => "%s: working copy URL does not exists at HEAD.\n",
);

# List of CLI prompt messages
our %CLI_MESSAGE_FOR_PROMPT = (
    CF_OVERWRITE      => qq{%s: existing changes will be overwritten.\n}
                         . qq{ Do you wish to continue?},
    CI                => qq{Would you like to commit this change?},
    CI_BRANCH_SHARED  => qq{\n}
                         . qq{*** WARNING: YOU ARE COMMITTING TO A %s BRANCH.\n}
                         . qq{*** Please ensure that you have the}
                         . qq{ owner's permission.\n\n}
                         . qq{Would you like to commit this change?},
    CI_BRANCH_USER    => qq{\n}
                         . qq{*** WARNING: YOU ARE COMMITTING TO A BRANCH}
                         . qq{ NOT OWNED BY YOU.\n}
                         . qq{*** Please ensure that you have the}
                         . qq{ owner's permission.\n\n}
                         . qq{Would you like to commit this change?},
    CI_TRUNK          => qq{\n}
                         . qq{*** WARNING: YOU ARE COMMITTING TO THE TRUNK.\n}
                         . qq{*** Please ensure that your change conforms to}
                         . qq{ your project's working practices.\n\n}
                         . qq{Would you like to commit this change?},
    CONTINUE          => qq{Are you sure you want to continue?},
    MERGE             => qq{Would you like to go ahead with the merge?},
    MERGE_REV         => qq{Please enter the revision you wish to merge from},
    MKPATCH_OVERWRITE => qq{%s: output location exists. OK to overwrite?},
    RUN_SVN_COMMAND   => qq{Would you like to run "svn %s"?},
);

# List of CLI warning messages
our %CLI_MESSAGE_FOR_WARNING = (
    BRANCH_SUBDIR   => "%s: is a sub-directory of a branch in a FCM project.\n",
    CF_BINARY       => "%s: ignoring binary file, please resolve manually.\n",
    INVALID_BRANCH  => $CLI_MESSAGE_FOR_ERROR{INVALID_BRANCH},
    ST_IN_TRAC_DIFF => "%s: local changes cannot be displayed in Trac.\n"
);

# CLI prompt handler and title prefix
our $CLI_PROMPT = \&_cli_prompt;
our $CLI_PROMPT_PREFIX = q{fcm };

# List of exception handlers [$class, CODE->($function, $e)]
our @CLI_EXCEPTION_HANDLERS = (
    ['Fcm::CLI::Exception', \&_cli_e_handler_of_cli_exception],
    ['Fcm::Cm::Exception' , \&_cli_e_handler_of_cm_exception],
    ['Fcm::Cm::Abort'     , \&_cli_e_handler_of_cm_abort],
);

# Event handlers
our %CLI_HANDLER_OF = (
    'WC_STATUS'      => \&_cli_handler_of_wc_status,
    'WC_STATUS_PATH' => \&_cli_handler_of_wc_status_path,
);

# Handlers of sub-commands
our %CLI_IMPL_OF = (
    'add'       => \&_cli_command_add,
    'branch'    => \&cm_branch,
    'commit'    => \&cm_commit,
    'conflicts' => \&cm_conflicts,
    'checkout'  => \&_cli_command_checkout,
    'delete'    => \&_cli_command_delete,
    'diff'      => \&cm_diff,
    'merge'     => \&cm_merge,
    'mkpatch'   => \&cm_mkpatch,
    'switch'    => \&_cli_command_switch,
    'update'    => \&_cli_command_update,
);

# List of overridden subcommands that need to display "svn help"
our %CLI_MORE_HELP_FOR = map {($_, 1)} qw{add diff delete switch update};

# The preferred name of subcommand aliases
our %CLI_PREFERRED_NAME_OF = (
    'ann'      => 'blame',
    'annotate' => 'blame',
    'br'       => 'branch',
    'ci'       => 'commit',
    'cf'       => 'conflicts',
    'co'       => 'checkout',
    'cp'       => 'copy',
    'del'      => 'delete',
    'di'       => 'diff',
    'ls'       => 'list',
    'mv'       => 'move',
    'pd'       => 'propdel',
    'pdel'     => 'propdel',
    'pe'       => 'propedit',
    'pedit'    => 'propedit',
    'pg'       => 'propget',
    'pget'     => 'propget',
    'pl'       => 'proplist',
    'plist'    => 'proplist',
    'praise'   => 'blame',
    'ps'       => 'propset',
    'pset'     => 'propset',
    'remove'   => 'delete',
    'ren'      => 'move',
    'rename'   => 'move',
    'rm'       => 'delete',
    'sw'       => 'switch',
    'up'       => 'update',
);

# List of subcommands that accept URL inputs
our %CLI_SUBCOMMAND_URL = map {($_, 1)} qw{
    blame
    branch
    cat
    checkout
    copy
    delete
    diff
    export
    import
    info
    list
    lock
    log
    merge
    mkdir
    mkpatch
    move
    propdel
    propedit
    propget
    proplist
    propset
    switch
    unlock
};

# List of subcommands that accept revision inputs
our %CLI_SUBCOMMAND_REV = map {($_, 1)} qw{
    blame
    branch
    cat
    checkout
    copy
    diff
    export
    info
    list
    log
    merge
    mkpatch
    move
    propdel
    propedit
    propget
    proplist
    propset
    switch
};

# Common patterns
our %PATTERN_OF = (
    # A CLI option
    CLI_OPT => qr{
        \A            (?# beginning)
        (--\w[\w-]*=) (?# capture 1, a long option label)
        (.*)          (?# capture 2, the value of the option)
        \z            (?# end)
    }xms,
    # A CLI revision option
    CLI_OPT_REV => qr{
        \A                      (?# beginning)
        (--revision(?:=|\z)|-r) (?# capture 1, --revision, --revision= or -r)
        (.*)                    (?# capture 2, trailing value)
        \z                      (?# end)
    }xms,
    # A CLI revision option range
    CLI_OPT_REV_RANGE => qr{
        \A                  (?# beginning)
        (                   (?# capture 1, begin)
            (?:\{[^\}]+\}+) (?# a date in curly braces)
            |               (?# or)
            [^:]+           (?# anything but a colon)
        )                   (?# capture 1, end)
        (?::(.*))?          (?# colon, and capture 2 til the end)
        \z                  (?# end)
    }xms,
    # A FCM branch path look-alike, should be configurable in the future
    FCM_BRANCH_PATH => qr{
        \A                            (?# beginning)
        /*                            (?# some slashes)
        (?:                           (?# group 1, begin)
            (?:trunk/*(?:@\d+)?\z)    (?# trunk at a revision)
            |                         (?# or)
            (?:trunk|branches|tags)/+ (?# trunk, branch or tags)
        )                             (?# group 1, end)
    }xms,
    # Last line of output from "svn status -u"
    ST_AGAINST_REV => qr{
        \A                           (?# beginning)
        Status\sagainst\srevision:.* (?# output of svn status -u)
        \z                           (?# end)
    }xms,
    # Extract path from "svn status"
    ST_PATH => qr{
        \A   (?# beginning)
        .{6} (?# 6 columns)
        \s+  (?# spaces)
        (.+) (?# capture 1, target path)
        \z   (?# end)
    }xms,
    # A legitimate "svn" revision
    SVN_REV => qr{
        \A                                      (?# beginning)
        (?:\d+|HEAD|BASE|COMMITTED|PREV|\{.+\}) (?# digit, reserved words, date)
        \z                                      (?# end)
    }ixms,
);

# Status matchers
our %ST_MATCHER_FOR = (
    MISSING     => sub {substr($_[0], 0, 1) eq '!'},
    MODIFIED    => sub {substr($_[0], 0, 6) =~ qr{\S}xms},
    OUT_OF_DATE => sub {substr($_[0], 7, 1) eq '*'},
    UNKNOWN     => sub {substr($_[0], 0, 1) eq '?'},
);

# ------------------------------------------------------------------------------
# Entry function for the FCM code management CLI. Calls the relevant FCM code
# management function or SVN command based on $function.
sub cli {
    my ($function, @args) = @_;
    if (exists($CLI_PREFERRED_NAME_OF{$function})) {
        $function = $CLI_PREFERRED_NAME_OF{$function};
    }
    if (grep {$_ eq '-h' || $_ eq '--help'} @args) {
        return _cli_help($function, 'NOEXIT');
    }
    if (exists($CLI_SUBCOMMAND_URL{$function})) {
        _cli_keyword_expand_url(\@args);
    }
    if (exists($CLI_SUBCOMMAND_REV{$function})) {
        _cli_keyword_expand_rev(\@args);
    }
    if (exists($CLI_IMPL_OF{$function})) {
        eval {
            local(@ARGV) = @args;
            return $CLI_IMPL_OF{$function}->(@args);
        };
        if ($@) {
            my $e = $@;
            for (@CLI_EXCEPTION_HANDLERS) {
                my ($class, $handler) = @{$_};
                if ($class->caught($e)) {
                    return $handler->($function, $e);
                }
            }
            die($e);
        }
    }
    else {
        return _svn($function, @args);
    }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   &Fcm::Cm::cm_branch ();
#
# DESCRIPTION
#   This is a FCM command to check information, create or delete a branch in
#   a Subversion repository.
# ------------------------------------------------------------------------------

sub cm_branch {
  # Process command line options
  # ----------------------------------------------------------------------------
  my (
    $info,
    $delete,
    $create,
    $list,
    $branch_of_branch,
    $name,
    $non_interactive,
    $password,
    $rev,
    $rev_flag,
    $show_all,
    $show_children,
    $show_other,
    $show_siblings,
    $svn_non_interactive,
    @tickets,
    $type,
    @userlist,
    $verbose,
  );
  my $rc = GetOptions(
    'info|i'              => \$info,
    'delete|d'            => \$delete,
    'create|c'            => \$create,
    'list|l'              => \$list,
    'branch-of-branch'    => \$branch_of_branch,
    'name|n=s'            => \$name,
    'non-interactive'     => \$non_interactive,
    'password=s'          => \$password,
    'revision|r=s'        => \$rev,
    'rev-flag=s'          => \$rev_flag,
    'show-all|a'          => \$show_all,
    'show-children'       => \$show_children,
    'show-other'          => \$show_other,
    'show-siblings'       => \$show_siblings,
    'svn-non-interactive' => \$svn_non_interactive,
    'ticket|k=s'          => \@tickets,
    'type|t=s'            => \$type,
    'user|u=s'            => \@userlist,
    'verbose|v'           => \$verbose,
  );
  if (!$rc) {
    _cli_err();
  }

  my $num_options = 0;
  $num_options++ if defined $info;
  $num_options++ if defined $delete;
  $num_options++ if defined $create;
  $num_options++ if defined $list;
  if ($num_options > 1) {
    _cli_err();
  }

  # Get URL of repository or branch
  # ----------------------------------------------------------------------------
  my $url;
  if ($ARGV[0]) {
    $url = Fcm::CmUrl->new (URL => $ARGV[0]);

    if (not $url->is_url) {
      # An argument is specified and is not a URL
      # Assume that it is a path with a working copy
      if (&is_wc ($ARGV[0])) {
        $url = Fcm::CmUrl->new (URL => &get_url_of_wc ($ARGV[0]));

      } else {
        return _cm_err(Fcm::Cm::Exception->INVALID_WC, $ARGV[0]);
      }
    }

  } else {
    # An argument is not specified
    # Assume that the current directory is a working copy
    if (&is_wc ()) {
      $url = Fcm::CmUrl->new (URL => &get_url_of_wc ());

    } else {
      return _cm_err(Fcm::Cm::Exception->INVALID_TARGET, '.');
    }
  }

  # Ensure $url->url_peg is a URL of a standard FCM project
  if (!$url->project_url()) {
    return _cm_err(Fcm::Cm::Exception->INVALID_PROJECT, $url->url_peg());
  }

  if ($create) {
    # The --create option is specified, create a branch
    # --------------------------------------------------------------------------

    # Check branch type flags
    if ($type) {
      $type = uc ($type);

      if ($type =~ /^(USER|SHARE)$/) {
        $type = 'DEV' . $Fcm::Config::DELIMITER . $1;

      } elsif ($type =~ /^(CONFIG|REL)$/) {
        $type = 'PKG' . $Fcm::Config::DELIMITER . $1;

      } elsif ($type =~ /^(DEV|TEST|PKG)$/) {
        $type = $1 . $Fcm::Config::DELIMITER . 'USER';

      } elsif ($type !~ /^(?:DEV|TEST|PKG)$Fcm::Config::DELIMITER(?:USER|SHARE)$/
               and $type !~ /^PKG$Fcm::Config::DELIMITER(?:CONFIG|REL)/) {
        _cli_err('CLI_OPT_ARG', 'type', $type);
      }

    } else {
      $type = 'DEV' . $Fcm::Config::DELIMITER . 'USER';
    }

    # Check branch name
    if (!$name) {
      _cli_err('CLI_OPT_WITH_OPT', 'name', 'create');
    }

    if ($name !~ qr{\A[\w.-]+\z}xms) {
      _cli_err('CLI_OPT_ARG', 'name', $name);
    }

    # Check revision flag is valid
    if ($rev_flag) {
      $rev_flag = uc ($rev_flag);
      if ($rev_flag !~ qr{\A (?:NORMAL|NUMBER|NONE) \z}xms) {
        _cli_err('CLI_OPT_ARG', 'rev-flag', $rev_flag);
      }

    } else {
      $rev_flag = 'NORMAL';
    }

    # Handle multiple tickets
    @tickets = split (
      /$Fcm::Config::DELIMITER_LIST/,
      join ($Fcm::Config::DELIMITER_LIST, @tickets)
    );
    s/^#// for (@tickets);
    @tickets = sort {$a <=> $b} @tickets;

    # Determine whether to create a branch of a branch
    $url->branch ('trunk') unless $branch_of_branch;

    # Create the branch
    my $branch = Fcm::CmBranch->new;
    $branch->create (
      SRC                 => $url,
      TYPE                => $type,
      NAME                => $name,
      PASSWORD            => $password,
      REV_FLAG            => $rev_flag,
      TICKET              => \@tickets,
      REV                 => $rev,
      NON_INTERACTIVE     => $non_interactive,
      SVN_NON_INTERACTIVE => $svn_non_interactive,
    );

  } elsif ($list) {
    # The option --list is specified
    # List branches owned by current or specified users
    # --------------------------------------------------------------------------
    # Get URL of the project "branches/" sub-directory
    $url->subdir ('');
    $url->branch ('');

    my @branches = $url->branch_list($rev);
    if (!$show_all) {
      @userlist = split(qr{:}xms, join(q{:}, @userlist));
      if (!@userlist) {
        @userlist = (Fcm::Config->instance()->user_id());
      }
      my %filter = map {($_, 1)} @userlist;
      @branches = grep {
        $filter{Fcm::CmBranch->new(URL => $_)->branch_owner()}
      } @branches
    }

    # Output, number of branches found
    $CLI_MESSAGE->(
      'BRANCH_LIST',
      $url->project_url_peg(),
      $rev ? "r$rev" : 'HEAD',
      scalar(@branches),
      ($show_all ? '[--show-all]' : join(q{, }, sort(@userlist))),
    );

    if (@branches) {
      # Output the URL of each branch
      if (not $verbose) {
        my $project = $url->project_url;
        @branches = map {Fcm::Keyword::unexpand($_)} @branches;
      }
      @branches = map {$_ . "\n"} sort @branches;
      $CLI_MESSAGE->(q{}, join(q{}, @branches));

    } else {
      # No branch found, exit with an error code
      return;
    }

  } else {
    # The option --info or --delete is specified
    # Report branch information (and/or delete a branch)
    # --------------------------------------------------------------------------
    # Set verbose level
    Fcm::Config->instance()->verbose ($verbose ? 1 : 0);

    # Set up the branch, report any error
    my $branch = Fcm::CmBranch->new (URL => $url->url_peg);
    if (!$branch->branch()) {
      return _cm_err(Fcm::Cm::Exception->INVALID_BRANCH, $branch->url_peg());
    }
    if (!$branch->url_exists()) {
      return _cm_err(Fcm::Cm::Exception->NOT_EXIST, $branch->url_peg());
    }

    # Remove the sub-directory part of the URL
    $branch->subdir ('');

    # Report branch info
    $branch->display_info (
      SHOW_CHILDREN => ($show_all || $show_children),
      SHOW_OTHER    => ($show_all || $show_other   ),
      SHOW_SIBLINGS => ($show_all || $show_siblings),
    );

    # Delete branch if --delete is specified
    $branch->del (
      PASSWORD            => $password,
      NON_INTERACTIVE     => $non_interactive,
      SVN_NON_INTERACTIVE => $svn_non_interactive,
    ) if $delete;
  }

}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   &Fcm::Cm::cm_commit ();
#
# DESCRIPTION
#   This is a FCM wrapper to the "svn commit" command.
# ------------------------------------------------------------------------------

sub cm_commit {
  my ($dry_run, $svn_non_interactive, $password);
  my $rc = GetOptions(
    'dry-run'             => \$dry_run,
    'svn-non-interactive' => \$svn_non_interactive,
    'password=s'          => \$password,
  );
  if (!$rc) {
    _cli_err();
  }

  # The remaining argument is the path to a working copy
  my ($path) = @ARGV;

  if ($path) {
    if (!-e $path) {
      return _cm_err(Fcm::Cm::Exception->NOT_EXIST, $path);
    }

  } else {
    # No argument specified, use current working directory
    $path = cwd ();
  }

  # Make sure we are in a working copy
  if (!is_wc($path)) {
    return _cm_err(Fcm::Cm::Exception->INVALID_WC, $path);
  }

  # Make sure we are at the top level of the working copy
  # (otherwise we might miss any template commit message)
  my $dir = &get_wct ($path);

  if ($dir ne cwd ()) {
    chdir($dir) || return _cm_err(Fcm::Cm::Exception->CHDIR, $dir);
    $CLI_MESSAGE->('CHDIR_WCT', $dir);
  }

  # Get update status of working copy
  # Check working copy files are not in conflict, missing, or out of date
  my @status = _svn_status_get([], 1);
  unless (defined $dry_run) {
    my (@conflict, @missing, @outdate);

    for (@status) {
      if (/^C/) {
        push @conflict, $_;
        next;
      }

      if (/^!/) {
        push @missing, $_;
        next;
      }

      if (/^.{7}\*/) {
        push @outdate, $_;
        next;
      }

      # Check that all files which have been added have the svn:executable
      # property set correctly (in case the developer adds a script before they
      # remember to set the execute bit)
      next unless /^A.{7} *\d+ +(.*)/;
      my $file = $1;

      next unless -f $file;
      my ($command, @arguments)
        = (-x $file && !-l $file) ? ('propset', '*') : ('propdel');
      run_command(['svn', $command, qw{-q svn:executable}, @arguments, $file]);
    }

    # Abort commit if files are in conflict, missing, or out of date
    if (@conflict or @missing or @outdate) {
      for (
        ['ST_CONFLICT'   , \@conflict],
        ['ST_MISSING'    , \@missing ],
        ['ST_OUT_OF_DATE', \@outdate ],
      ) {
        my ($key, $array_ref) = @{$_};
        if (@{$array_ref}) {
          $CLI_MESSAGE->($key, join(q{}, @{$array_ref}));
        }
      }
      return _cm_abort(Fcm::Cm::Abort->FAIL);
    }
  }

  # Read in any existing message
  my $ci_mesg = Fcm::CmCommitMessage->new;
  $ci_mesg->read_file;

  # Execute "svn status" for a list of changed items
  @status = grep !/^\?/, _svn_status_get();

  # Abort if there is no change in the working copy
  if (!@status) {
    return _cm_abort(Fcm::Cm::Abort->NULL);
  }

  # Get associated URL of current working copy
  my $url = Fcm::CmUrl->new (URL => &get_url_of_wc ());

  # Include URL, or project, branch and sub-directory info in @status
  unshift @status, "\n";

  if ($url->project and $url->branch) {
    unshift @status, (
      '[Project: ' . $url->project                           . ']' . "\n",
      '[Branch : ' . $url->branch                            . ']' . "\n",
      '[Sub-dir: ' . ($url->subdir ? $url->subdir : '<top>') . ']' . "\n",
    );

  } else {
    unshift @status, '[URL: ' . $url->url . ']' . "\n";
  }

  # Use a temporary file to store the final commit log message
  $ci_mesg->ignore_mesg (\@status);
  my $logfile = $ci_mesg->edit_file (TEMP => 1);

  # Check with the user to see if he/she wants to go ahead
  my $reply = 'n';
  if (!defined($dry_run)) {
    # Add extra warning for trunk commit
    my @prompt_args;
    my $user = Fcm::Config->instance()->user_id();

    if ($url->is_trunk()) {
      @prompt_args = ('CI_TRUNK');
    }
    elsif ($user && $url->is_branch() && $url->branch_owner() ne $user) {
      if (exists $Fcm::CmUrl::owner_keywords{$url->branch_owner}) {
        @prompt_args = (
          'CI_BRANCH_SHARED',
          uc($Fcm::CmUrl::owner_keywords{$url->branch_owner()}),
        );
      }
      else {
        @prompt_args = ('CI_BRANCH_USER');
      }
    }
    else {
      @prompt_args = ('CI');
    }
    $reply = $CLI_PROMPT->('commit', @prompt_args);
  }

  if ($reply eq 'y') {
    # Commit the change if user replies "y" for "yes"
    my @command = (
      qw/svn commit -F/, $logfile,
      ($svn_non_interactive  ? '--non-interactive'       : ()),
      (defined $password     ? ('--password', $password) : ()),
    );
    my $rc;
    &run_command (\@command, RC => \$rc, ERROR => 'warn');

    if ($rc) {
      # Commit failed
      # Write temporary commit log content to commit log message file
      $ci_mesg->write_file;

      # Fail the command
      return _cm_abort(Fcm::Cm::Abort->FAIL);
    }

    # Remove commit message file
    unlink $ci_mesg->file;

    # Update the working copy
    $CLI_MESSAGE->(q{}, join(q{}, _svn_update()));

  } else {
    $ci_mesg->write_file;
    if (!$dry_run) {
      return _cm_abort();
    }
  }

  return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   &Fcm::Cm::cm_conflicts ();
#
# DESCRIPTION
#   This is a FCM command for resolving conflicts within working copy using a
#   graphical merge tool.
# ------------------------------------------------------------------------------

sub cm_conflicts {
  # Path to the working copy
  my $path = $ARGV[0];
  $path    = cwd () if not $path;

  # Check for any files with conflicts
  my @status = grep /^C.{4} *(.*)/, &run_command (
    [qw/svn st/, ($path eq cwd () ? () : $path)], METHOD => 'qx',
  );
  my @files  = map {m/^C.{4} *(.*)/; $1} @status;

  # Save current working directory
  my $topdir = cwd ();

  # Set up environment for graphical merge
  # Use environment variable if set, otherwise use default setting
  local(%ENV) = %ENV;
  $ENV{FCM_GRAPHIC_MERGE}
    ||= Fcm::Config->instance()->setting (qw/TOOL GRAPHIC_MERGE/);

  FILE:
  for my $file (@files) {
    # Print name of file in conflicts
    $CLI_MESSAGE->('CF', $file);

    # Determine directory and base name of file in conflicts
    my $base = basename $file;
    my $dir  = dirname $file;

    # Change to container directory of file in conflicts
    chdir(File::Spec->catfile($topdir, $dir))
      || return _cm_err(Fcm::Cm::Exception->CHDIR, $dir);

    # Use "svn info" to determine conflict marker files
    my @info = &run_command ([qw/svn info/, $base], METHOD => 'qx');

    # Ignore if $base is a binary file
    if (-B $base) {
      $CLI_MESSAGE->('CF_BINARY', $base);
      next FILE;
    }

    # Get conflicts markers files
    my ($older, $mine, $yours);

    for (@info) {
      $older = $1 if (/^Conflict Previous Base File: (.*)/);
      $mine  = $1 if (/^Conflict Previous Working File: (.*)/);
      $yours = $1 if (/^Conflict Current Base File: (.*)/);
    }

    if (-f $base and (stat $base)[9] > (stat $mine)[9] + 1) {
      # If $base is newer (by more than a second), it may contain saved changes
      if ($CLI_PROMPT->('conflicts', 'CF_OVERWRITE', $base) ne 'y') {
        next FILE;
      }
    }

    # Launch graphic merge tool
    my $rc;
    my $command = [qw/fcm_graphic_merge/, $base, $mine, $older, $yours];
    # $rc == 0: all conflicts resovled
    # $rc == 1: some conflicts not resolved
    # $rc == 2: trouble
    eval {
      run_command($command, RC => \$rc);
    };
    if ($@) {
      if (!defined($rc) || $rc > 1) {
        die($@);
      }
    }
    next FILE if $rc;

    # Prompt user to run "svn resolved" on the file
    if ($CLI_PROMPT->('conflicts', 'RUN_SVN_COMMAND', 'resolved') eq 'y') {
      run_command([qw{svn resolved}, $base]);
    }
  }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   &Fcm::Cm::cm_diff ();
#
# DESCRIPTION
#   This is a wrapper to "svn diff". It adds two extra functionalities. The
#   first one allows the command to show differences relative to the base of
#   the branch. The second one allows differences to be displayed via a
#   graphical tool.
# ------------------------------------------------------------------------------

sub cm_diff {
  # Set up environment for graphical diff
  # Use environment variable if set, otherwise use default setting
  local(%ENV) = %ENV;
  $ENV{FCM_GRAPHIC_DIFF}
    ||= Fcm::Config->instance()->setting(qw/TOOL GRAPHIC_DIFF/);

  # Check for the --branch options
  # ----------------------------------------------------------------------------
  my $branch = grep {$_ eq '-b' or $_ eq '--branch'} @ARGV;

  if (not $branch) {
    # The --branch option not specified, just call "svn diff"
    # Convert the --graphical to qw/--diff-cmd fcm_graphical_diff/
    # Convert the --summarise to --summarize
    @ARGV = map {
      my @return;
      if ($_ eq '-g' or $_ eq '--graphical') {
        @return = (qw/--diff-cmd fcm_graphic_diff/)

      } elsif ($_ eq '--summarise') {
        @return = ('--summarize');

      } else {
        @return = ($_);
      }
      @return;
    } @ARGV;

    # Execute the command
    return _svn('diff', @ARGV);
  }

  # The --branch option is specified
  # ----------------------------------------------------------------------------

  # Determine whether the --graphical option is specified,
  # if so set the appropriate command
  # ----------------------------------------------------------------------------
  my ($diff_cmd, $extensions, $graphical, $summarise, $trac, $wiki);
  my $rc = GetOptions (
    'b|branch'            => \$branch,
    'diff-cmd=s'          => \$diff_cmd,
    'x|extensions=s'      => \$extensions,
    'g|graphical'         => \$graphical,
    'summarise|summarize' => \$summarise,
    't|trac'              => \$trac,
    'wiki'                => \$wiki,
  );
  if (!$rc) {
    _cli_err();
  }

  my @diff_cmd = ();

  if ($graphical) {
    @diff_cmd = (qw/--diff-cmd fcm_graphic_diff/);

  } elsif ($diff_cmd) {
    @diff_cmd = ('--diff-cmd', $diff_cmd);

    push @diff_cmd, '--extensions', split (/\s+/, $extensions) if $extensions;
  }

  # The remaining argument should either be a URL or a PATH
  my ($url_arg, $path_arg);

  if (@ARGV) {
    my $arg = Fcm::CmUrl->new (URL => $ARGV[0]);

    if ($arg->is_url) {
      $url_arg = $ARGV[0];

    } else {
      $path_arg = $ARGV[0];
    }
  }

  # Get repository and branch information
  # ----------------------------------------------------------------------------
  my ($url, $path);
  if (defined $url_arg) {
    # If a URL is specified, get repository and branch information from it
    $url = Fcm::CmBranch->new (URL => $url_arg);

  } else {
    # Get repository and branch information from the specified path or the
    # current directory if it is a working copy
    $path = $path_arg ? $path_arg : cwd ();
    if (!is_wc($path)) {
      return _cm_err(Fcm::Cm::Exception->INVALID_WC, $path);
    }

    $url  = Fcm::CmBranch->new (URL => &get_url_peg_of_wc ($path));
  }

  # Check that URL is a standard FCM branch
  if (!$url->is_branch()) {
    return _cm_err(Fcm::Cm::Exception->INVALID_BRANCH, $url->url_peg());
  }

  # Save and remove sub-directory part of the URL
  my $subdir = $url->subdir ();
  $url->subdir ('');

  # Check that $url exists
  if (!$url->url_exists()) {
    return _cm_err(Fcm::Cm::Exception->INVALID_URL, $url->url_peg());
  }

  # Compare current branch with its parent
  # ----------------------------------------------------------------------------
  my $parent = Fcm::CmBranch->new (URL => $url->parent->url);
  $parent->pegrev ($url->pegrev) if $url->pegrev;

  if (!$parent->url_exists()) {
    return _cm_err(
      Fcm::Cm::Exception->PARENT_NOT_EXIST, $url->url_peg(), $parent->url(),
    );
  }

  my $base = $parent->base_of_merge_from ($url);

  # Ensure the correct diff (syntax) is displayed
  # ----------------------------------------------------------------------------
  # Reinstate the sub-tree part into the URL
  $url->subdir ($subdir);
  $base->subdir ($subdir);

  # Ensure the branch URL has a peg revision
  $url->pegrev ($url->svninfo (FLAG => 'Last Changed Rev')) if not $url->pegrev;

  if ($trac or $wiki) {
    # Trac/wiki
    # --------------------------------------------------------------------------
    if (!$url_arg && _svn_status_get([$path_arg ? $path_arg : q{.}])) {
      $CLI_MESSAGE->('ST_IN_TRAC_DIFF', ($path_arg ? $path_arg : q{.}));
    }

    # Trac wiki syntax
    my $wiki_syntax = 'diff:' . $base->path_peg . '//' . $url->path_peg;

    if ($wiki) {
      # Print Trac wiki syntax only
      $CLI_MESSAGE->(q{}, "$wiki_syntax\n");

    } else { # if $trac
      # Use Trac to view "diff"
      my $browser = Fcm::Config->instance()->setting(qw/WEB_BROWSER/);
      $browser ||= 'firefox';

      my $trac_url = Fcm::Keyword::get_browser_url($url->project_url());
      $trac_url =~ s{/intertrac/.*$}{/intertrac/$wiki_syntax}xms;

      &run_command ([$browser, $trac_url], METHOD => 'exec', PRINT => 1);
    }

  } else {
    # Execute the "diff" command
    # --------------------------------------------------------------------------
    my @command = (
      qw/svn diff/, @diff_cmd,
      ($summarise ? ('--summarize') : ()),
      '--old', $base->url_peg,
      '--new', ($url_arg ? $url->url_peg : ($path_arg ? $path_arg : '.')),
    );
    &run_command (\@command, PRINT => 1);
  }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   &Fcm::Cm::cm_merge ();
#
# DESCRIPTION
#   This is a wrapper to "svn merge".
# ------------------------------------------------------------------------------

sub cm_merge {
  # Options
  # ----------------------------------------------------------------------------
  my ($custom, $dry_run, $non_interactive, $reverse, $rev, $verbose);
  my $rc = GetOptions(
    'custom'          => \$custom,
    'dry-run'         => \$dry_run,
    'non-interactive' => \$non_interactive,
    'reverse'         => \$reverse,
    'revision|r=s'    => \$rev,
    'verbose|v'       => \$verbose,
  );
  if (!$rc) {
    _cli_err();
  }

  # Find out the URL of the working copy
  # ----------------------------------------------------------------------------
  my ($target, $wct);
  if (&is_wc ()) {
    $wct = &get_wct ();

    if ($wct ne cwd ()) {
      chdir($wct) || return _cm_err(Fcm::Cm::Exception->CHDIR, $wct);
      $CLI_MESSAGE->('CHDIR_WCT', $wct);
    }

    $target = Fcm::CmBranch->new (URL => &get_url_of_wc ($wct));

  } else {
    return _cm_err(Fcm::Cm::Exception->INVALID_WC, '.');
  }

  if (!$target->url_exists()) {
    return _cm_err(Fcm::Cm::Exception->WC_URL_NOT_EXIST, '.');
  }

  # The target must be at the top of a branch
  # $subdir will be used later to determine whether the merge is allowed or not
  my $subdir = $target->subdir;
  $target->subdir ('') if $subdir;

  # Check for any local modifications
  # ----------------------------------------------------------------------------
  if (!$dry_run && !$non_interactive) {
    _svn_status_checker('merge', 'MODIFIED', $CLI_HANDLER_OF{WC_STATUS})->();
  }

  # Determine the SOURCE URL
  # ----------------------------------------------------------------------------
  my $source;

  if ($reverse) {
    # Reverse merge, the SOURCE is the the working copy URL
    $source = Fcm::CmBranch->new (URL => $target->url);

  } else {
    # Automatic/custom merge, argument 1 is the SOURCE of the merge
    my $source_url = shift (@ARGV);
    if (!$source_url) {
      _cli_err('CLI_MERGE_ARG1');
    }

    $source = _cm_get_source($source_url, $target);
  }

  # Parse the revision option
  # ----------------------------------------------------------------------------
  if ($reverse && !$rev) {
    _cli_err('CLI_OPT_WITH_OPT', 'revision', 'reverse');
  }
  my @revs = (($reverse || $custom) && $rev ? split(qr{:}xms, $rev) : ());

  # Determine the merge delta and the commit log message
  # ----------------------------------------------------------------------------
  my (@delta, $mesg);
  my $separator = '-' x 80 . "\n";

  if ($reverse) {
    # Reverse merge
    # --------------------------------------------------------------------------
    if (@revs == 1) {
      $revs[1] = ($revs[0] - 1);

    } else {
      @revs = sort {$b <=> $a} @revs;
    }

    $source->pegrev ($source->svninfo (FLAG => 'Last Changed Rev'))
      unless $source->pegrev;
    $source->subdir ($subdir);

    # "Delta" of the "svn merge" command
    @delta = ('-r' . $revs[0] . ':' . $revs[1], $source->url_peg);

    # Template message
    $mesg = 'Reversed r' . $revs[0] .
            (($revs[1] < $revs[0] - 1) ? ':' . $revs[1] : '') . ' of ' .
            $source->path . "\n";

  } elsif ($custom) {
    # Custom merge
    # --------------------------------------------------------------------------
    if (@revs) {
      # Revision specified
      # ------------------------------------------------------------------------
      # Only one revision N specified, use (N - 1):N as the delta
      unshift @revs, ($revs[0] - 1) if @revs == 1;

      $source->pegrev ($source->svninfo (FLAG => 'Last Changed Rev'))
        unless $source->pegrev;
      $source->subdir ($subdir);
      $target->subdir ($subdir);

      # "Delta" of the "svn merge" command
      @delta = ('-r' . $revs[0] . ':' . $revs[1], $source->url_peg);

      # Template message
      $mesg = 'Custom merge into ' . $target->path . ': r' . $revs[1] .
              ' cf. r' . $revs[0] . ' of ' . $source->path_peg . "\n";

    } else {
      # Revision not specified
      # ------------------------------------------------------------------------
      # Get second source URL
      my $source2_url = shift (@ARGV);
      if (!$source2_url) {
        _cli_err('CLI_MERGE_ARG2');
      }

      my $source2 = _cm_get_source($source2_url, $target);

      $source->pegrev  ($source->svninfo  (FLAG => 'Last Changed Rev'))
        unless $source->pegrev;
      $source2->pegrev ($source2->svninfo (FLAG => 'Last Changed Rev'))
        unless $source2->pegrev;
      $source->subdir  ($subdir);
      $source2->subdir ($subdir);
      $target->subdir  ($subdir);

      # "Delta" of the "svn merge" command
      @delta = ($source->url_peg, $source2->url_peg);

      # Template message
      $mesg = 'Custom merge into ' . $target->path . ': ' . $source->path_peg .
              ' cf. ' . $source2->path_peg . "\n";
    }

  } else {
    # Automatic merge
    # --------------------------------------------------------------------------
    # Check to ensure source branch is not the same as the target branch
    if (!$target->branch()) {
      return _cm_err(Fcm::Cm::Exception->WC_INVALID_BRANCH, $wct);
    }
    if ($source->branch() eq $target->branch()) {
      return _cm_err(Fcm::Cm::Exception->MERGE_SELF, $target->url_peg(), $wct);
    }

    # Only allow the merge if the source and target are "directly related"
    # --------------------------------------------------------------------------
    my $anc = $target->ancestor ($source);
    return _cm_err(
      Fcm::Cm::Exception->MERGE_UNRELATED, $target->url_peg(), $source->url_peg
    ) unless
      ($anc->url eq $target->url and $anc->url_peg eq $source->parent->url_peg)
      or
      ($anc->url eq $source->url and $anc->url_peg eq $target->parent->url_peg)
      or
      ($anc->url eq $source->parent->url and $anc->url eq $target->parent->url);

    # Check for available merges from the source
    # --------------------------------------------------------------------------
    my @revs = $target->avail_merge_from ($source, 1);

    if (@revs) {
      if ($verbose) {
        # Verbose mode, print log messages of available merges
        $CLI_MESSAGE->('MERGE_REVS', $source->path_peg(), q{});
        for (@revs) {
          $CLI_MESSAGE->('SEPARATOR');
          $CLI_MESSAGE->(q{}, $source->display_svnlog($_));
        }
        $CLI_MESSAGE->('SEPARATOR');
      }
      else {
        # Normal mode, list revisions of available merges
        $CLI_MESSAGE->('MERGE_REVS', $source->path_peg(), join(q{ }, @revs));
      }

    } else {
      return _cm_abort(Fcm::Cm::Abort->NULL);
    }

    # If more than one merge available, prompt user to enter a revision number
    # to merge from, default to $revs [0]
    # --------------------------------------------------------------------------
    if ($non_interactive || @revs == 1) {
      $source->pegrev($revs[0]);
    }
    else {
      my $reply = $CLI_PROMPT->(
        {type => q{}, default => $revs[0]}, 'merge', 'MERGE_REV',
      );
      if (!defined($reply)) {
        return _cm_abort();
      }
      # Expand revision keyword if necessary
      if ($reply) {
        $reply = (Fcm::Keyword::expand($target->project_url(), $reply))[1];
      }
      # Check that the reply is a number in the available merges list
      if (!grep {$_ eq $reply} @revs) {
        return _cm_err(Fcm::Cm::Exception->MERGE_REV_INVALID, $reply)
      }
      $source->pegrev($reply);
    }

    # If the working copy top is pointing to a sub-directory of a branch,
    # we need to check whether the merge will result in losing changes made in
    # other sub-directories of the source.
    if ($subdir and not $target->allow_subdir_merge_from ($source, $subdir)) {
      return _cm_err(Fcm::Cm::Exception->MERGE_UNSAFE, $source->url_peg());
    }

    # Calculate the base of the merge
    my $base = $target->base_of_merge_from ($source);

    # $source and $base must take into account the sub-directory
    my $s = Fcm::CmBranch->new (URL => $source->url_peg);
    my $b = Fcm::CmBranch->new (URL => $base->url_peg);

    $s->subdir ($subdir) if $subdir;
    $b->subdir ($subdir) if $subdir;

    # Diagnostic
    $CLI_MESSAGE->('MERGE_CF', $s->path_peg(), $b->path_peg());

    # Delta of the "svn merge" command
    @delta = ($b->url_peg, $s->url_peg);

    # Template message
    $mesg = 'Merged into ' . $target->path . ': ' . $source->path_peg .
            ' cf. ' . $base->path_peg . "\n";
  }

  # Run "svn merge" in "--dry-run" mode to see the result
  # ----------------------------------------------------------------------------
  my @out   = &run_command (
    [qw/svn merge --dry-run/, @delta],
    METHOD => 'qx', PRINT => ($dry_run and $verbose),
  );

  # Abort merge if it will result in no change
  if (not @out) {
    return _cm_abort(Fcm::Cm::Abort->NULL);
  }

  # Report result of "svn merge --dry-run"
  if ($dry_run || !$non_interactive) {
    $CLI_MESSAGE->('MERGE_DRY');
    $CLI_MESSAGE->('SEPARATOR');
    $CLI_MESSAGE->(q{}, join(q{}, @out));
    $CLI_MESSAGE->('SEPARATOR');
  }

  return if $dry_run;

  # Prompt the user to see if (s)he would like to go ahead
  # ----------------------------------------------------------------------------
  # Go ahead with merge only if user replies "y"
  if (!$non_interactive && $CLI_PROMPT->('merge', 'MERGE') ne 'y') {
    return _cm_abort();
  }
  $CLI_MESSAGE->('MERGE');
  run_command([qw/svn merge/, @delta], PRINT => $verbose);

  # Prepare the commit log
  # ----------------------------------------------------------------------------
  # Read in any existing message
  my $ci_mesg = Fcm::CmCommitMessage->new;
  $ci_mesg->read_file;
  $ci_mesg->auto_mesg ([$mesg, @{ $ci_mesg->auto_mesg }]);
  $ci_mesg->write_file;

  if ($verbose) {
    $CLI_MESSAGE->('SEPARATOR');
    $CLI_MESSAGE->('MERGE_CI', $mesg);
  }

  return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   &Fcm::Cm::cm_mkpatch ();
#
# DESCRIPTION
#   This is a FCM command to create a patching script from particular revisions
#   of a URL.
# ------------------------------------------------------------------------------

sub cm_mkpatch {
  # Process command line options and arguments
  # ----------------------------------------------------------------------------
  my (@exclude, $organisation, $revision);
  my $rc = GetOptions(
    'exclude=s'      => \@exclude,
    'organisation=s' => \$organisation,
    'r|revision=s'   => \$revision,
  );
  if (!$rc) {
    _cli_err();
  }

  # Excluded paths, convert glob into regular patterns
  @exclude = split (/:/, join (':', @exclude));
  for (@exclude) {
    s#\*#[^/]*#; # match any number of non-slash character
    s#\?#[^/]#;  # match a non-slash character
    s#/*$##;     # remove trailing slash
  }

  # Organisation prefix
  $organisation = $organisation ? $organisation : 'original';

  # Make sure revision option is set correctly
  my @revs = $revision ? split (/:/, $revision) : ();
  @revs    = @revs [0, 1] if @revs > 2;

  # Arguments
  my ($u, $outdir) = @ARGV;

  if (!$u) {
    _cli_err();
  }

  my $url = Fcm::CmUrl->new (URL => $u);
  if (!$url->is_url()) {
    return _cm_err(Fcm::Cm::Exception->INVALID_URL, $u);
  }
  if (!$url->url_exists()) {
    return _cm_err(Fcm::Cm::Exception->NOT_EXIST, $u);
  }
  if (!$url->branch()) {
    $CLI_MESSAGE->('INVALID_BRANCH', $u);
  }
  elsif ($url->subdir()) {
    $CLI_MESSAGE->('BRANCH_SUBDIR', $u);
  }

  if (@revs) {
    # If HEAD revision is given, convert it into a number
    # --------------------------------------------------------------------------
    for my $rev (@revs) {
      $rev = $url->svninfo (FLAG => 'Revision') if uc ($rev) eq 'HEAD';
    }

  } else {
    # If no revision is given, use the HEAD
    # --------------------------------------------------------------------------
    $revs[0] = $url->svninfo (FLAG => 'Revision');
  }

  $revs[1] = $revs[0] if @revs == 1;

  # Check that output directory is set
  # ----------------------------------------------------------------------------
  $outdir = File::Spec->catfile (cwd (), 'fcm-mkpatch-out') if not $outdir;

  if (-e $outdir) {
    # Ask user to confirm removal of old output directory if it exists
    if ($CLI_PROMPT->('mkpatch', 'MKPATCH_OVERWRITE') ne 'y') {
      return _cm_abort();
    }

    rmtree($outdir) || return _cm_err(Fcm::Cm::Exception->RMTREE, $outdir);
  }

  # (Re-)create output directory
  mkpath($outdir) || return _cm_err(Fcm::Cm::Exception->MKPATH, $outdir);
  $CLI_MESSAGE->('OUT_DIR', $outdir);

  # Get and process log of URL
  # ----------------------------------------------------------------------------
  my @script   = (); # main output script
  my %log      = $url->svnlog (REV => \@revs);
  my $url_path = $url->path;

  for my $rev (sort {$a <=> $b} keys %log) {
    # Look at the changed paths for each revision
    my $use_patch = 1;  # OK to use a patch file?
    my @paths;
    PATH: for my $path (sort keys %{ $log{$rev}{paths} }) {
      my $file = $path;

      # Skip paths outside of the branch
      next PATH unless $file =~ s#^$url_path/*##;

      # Skip excluded paths
      for my $exclude (@exclude) {
        if ($file =~ m#^$exclude(?:/*|$)#) {
          # Can't use a patch file if any files have been excluded
          $use_patch = 0;
          next PATH;
        }
      }

      # Can't use a patch file if any files have been added or replaced
      $use_patch = 0 if $log{$rev}{paths}{$path}{action} eq 'A' or
                        $log{$rev}{paths}{$path}{action} eq 'R';

      push @paths, $path;
    }

    # If a patch is being used, make sure it isn't just property changes
    if ($use_patch) {
      my @changedpaths;
      for my $path (@paths) {
        (my $file = $path) =~ s#^$url_path/*##;
        if ($log{$rev}{paths}{$path}{action} eq 'M') {
          my ($diff) = &run_command (
                         [qw/svn diff --no-diff-deleted --summarize -c/,
                          $rev, $url->url . '/' . $file. '@' . $rev],
                         METHOD => 'qx');
          next unless $diff =~ /^[A-Z]/;
        }
        push @changedpaths, $path;
      }
      @paths = @changedpaths;
    }

    next unless @paths;

    # Create the patch using "svn diff"
    my @patch = ();
    if ($use_patch) {
      @patch = &run_command ([qw/svn diff --no-diff-deleted -c/, $rev,
                              $url->url], METHOD => 'qx');
      if (@patch) {
        # Don't use the patch if it may contain subversion keywords
        for (@patch) {
          $use_patch = 0 if /\$[a-zA-Z:]+ *\$/;
        }
      } else {
        $use_patch = 0;
      }
    }

    # Create a directory for this revision in the output directory
    my $outdir_rev = File::Spec->catfile ($outdir, $rev);
    mkpath($outdir_rev)
      || return _cm_err(Fcm::Cm::Exception->MKPATH, $outdir_rev);

    # Parse commit log message
    my @msg = split /\n/, $log{$rev}{msg};
    for (@msg) {
      # Re-instate line break
      $_ .= "\n";

      # Remove line if it matches a merge template
      $_ = '' if /^Reversed r\d+(?::\d+)? of \S+$/;
      $_ = '' if /^Custom merge into \S+:.+$/;
      $_ = '' if /^Merged into \S+: \S+ cf\. \S+$/;

      # Modify Trac ticket link
      s/(?:#|ticket:)(\d+)/${organisation}_ticket:$1/g;

      # Modify Trac changeset link
      s/(?:r|changeset:)(\d+)/${organisation}_changeset:$1/g;
      s/\[(\d+)\]/${organisation}_changeset:$1/g;
    }

    push @msg, '(' . $organisation . '_changeset:' . $rev . ')' . "\n";

    # Write commit log message in a file
    my $f_revlog = File::Spec->catfile ($outdir_rev, 'log-message');
    open FILE, '>', $f_revlog or die $f_revlog, ': cannot open (', $!, ')';
    print FILE @msg;
    close FILE or die $f_revlog, ': cannot close (', $!, ')';

    # Handle each changed path
    my $export_file   = 1;  # name for next exported file (gets incremented)
    my $patch_needed  = 0;  # is a patch file required?
    my @before_script = (); # patch script to run before patch applied
    my @after_script  = (); # patch script to run after patch applied
    my @copied_dirs   = (); # copied directories
    CHANGED: for my $path (@paths) {
      (my $file = $path) =~ s#^$url_path/*##;
      my $url_file = $url->url . '/' . $file . '@' . $rev;

      # Skip paths within copied directories
      for my $copied_dir (@copied_dirs) {
        next CHANGED if $file =~ m#^$copied_dir(?:/*|$)#;
      }

      if ($log{$rev}{paths}{$path}{action} eq 'D') {
        # Script to delete file
        push @after_script, 'svn delete ' . $file;

      } else {
        my $export_required = 0;
        my $recursive_add   = 0;
        my $is_newfile      = 0;

        # Skip property changes
        if ($log{$rev}{paths}{$path}{action} eq 'M') {
          my ($diff) = &run_command (
                         [qw/svn diff --no-diff-deleted --summarize -c/,
                          $rev, $url->url . '/' . $file. '@' . $rev],
                         METHOD => 'qx');
          next CHANGED unless $diff =~ /^[A-Z]/;
        }

        # Determine if the file is a directory
        my $is_dir = 0;
        if ($log{$rev}{paths}{$path}{action} ne 'M') {
          my @info = &run_command ([qw/svn info/, $url_file], METHOD => 'qx');
          for (@info) {
            if (/^Node Kind: (\w+)/) {
              $is_dir = 1 if $1 eq 'directory';
              last;
            }
          }
        }

        # Decide how to treat added files
        if ($log{$rev}{paths}{$path}{action} eq 'A') {
          # Determine if the file is copied
          if (exists $log{$rev}{paths}{$path}{'copyfrom-path'}) {
            if ($is_dir) {
              # A copied directory needs to be treated as a new file, exported
              # and added recursively
              $is_newfile      = 1;
              $export_required = 1;
              $recursive_add   = 1;
              push @copied_dirs, $file;
            } else {
              # History exists for this file
              my $copyfrom_path = $log{$rev}{paths}{$path}{'copyfrom-path'};
              my $copyfrom_rev  = $log{$rev}{paths}{$path}{'copyfrom-rev'};
              my $cp_url = Fcm::CmUrl->new (
                URL => $url->root . $copyfrom_path . '@' . $copyfrom_rev,
              );

              if ($copyfrom_path =~ s#^$url_path/*##) {
                # File is copied from a file under the specified URL
                # Check source exists
                $is_newfile = 1 unless $cp_url->url_exists ($rev - 1);
              } else {
                # File copied from outside of the specified URL
                $is_newfile = 1;

                # Check branches can be determined
                if ($url->branch and $cp_url->branch) {

                  # Follow its history, stop on copy
                  my %cp_log = $cp_url->svnlog (STOP_ON_COPY => 1);

                  # "First" revision of the copied file
                  my $cp_rev = (sort {$a <=> $b} keys %cp_log) [0];
                  my %attrib = %{ $cp_log{$cp_rev}{paths}{$cp_url->path} }
                    if $cp_log{$cp_rev}{paths}{$cp_url->path};

                  # Check whether the "first" revision is copied from elsewhere.
                  if (exists $attrib{'copyfrom-path'}) {
                    # If source exists in the specified URL, set up the copy
                    my $cp_cp_url = Fcm::CmUrl->new (
                      URL => $url->root . $attrib{'copyfrom-path'} . '@' .
                             $attrib{'copyfrom-rev'},
                    );
                    $cp_cp_url->branch ($url->branch);
                    if ($cp_cp_url->url_exists ($rev - 1)) {
                      ($copyfrom_path = $cp_cp_url->path) =~ s#^$url_path/*##;
                      # Check path is defined - if not it probably means the
                      # branch doesn't follow the FCM naming convention
                      $is_newfile = 0 if $copyfrom_path;
                    }
                  }

                  # Note: The logic above does not cover all cases. However, it
                  # should do the right thing for the most common case. Even
                  # where it gets it wrong the file contents should always be
                  # correct even if the file history is not.
                }
              }

              # Check whether file is copied from an excluded path
              if (not $is_newfile) {
                for my $exclude (@exclude) {
                  if ($copyfrom_path =~ m#^$exclude(?:/*|$)#) {
                    $is_newfile = 1;
                    last;
                  }
                }
              }

              # Script to copy file, if required
              push @before_script, 'svn copy ' . $copyfrom_path .  ' ' . $file
                if not $is_newfile;
            }

          } else {
            # History does not exist, must be a new file
            $is_newfile = 1;
            # If it's a directory then create it (in case patch doesn't)
            push @before_script, 'mkdir ' . $file if $is_dir;
          }
        }

        if ($log{$rev}{paths}{$path}{action} eq 'R') {
          # Script to delete file
          push @before_script, 'svn delete ' . $file;

          # Now treat as new file
          $is_newfile = 1;
        }

        # Script to add the file, if required
        if ($is_newfile) {
          if ($recursive_add) {
            push @after_script, 'svn add ' . $file;
          } else {
            push @after_script, 'svn add --non-recursive ' . $file;
          }
        }

        # Decide whether the file needs to be exported
        if (not $is_dir) {
          if (not $use_patch) {
            $export_required = 1;
          } else {
            # Export the file if it is binary
            my @mime_type = &run_command
             ([qw/svn propget svn:mime-type/, $url_file], METHOD => 'qx');
            for (@mime_type) {
              $export_required = 1 if not /^text\//;
            }
            # Only create a patch file if necessary
            $patch_needed = 1 if not $export_required;
          }
        }

        if ($export_required) {
          # Download the file using "svn export"
          my $export = File::Spec->catfile ($outdir_rev, $export_file);
          &run_command ([qw/svn export -q -r/, $rev, $url_file, $export]);

          # Copy the exported file into the file
          push @before_script,
               'cp -r ${fcm_patch_dir}/' . $export_file . ' ' . $file;
          $export_file++;
        }
      }
    }

    # Write the patch file
    if ($patch_needed) {
      my $patchfile = File::Spec->catfile ($outdir_rev, 'patchfile');
      open FILE, '>', $patchfile
        or die $patchfile, ': cannot open (', $!, ')';
      print FILE @patch;
      close FILE or die $patchfile, ': cannot close (', $!, ')';
    }

    # Add line break to each line in @before_script and @after_script
    @before_script = map {($_ ? $_ . ' || exit 1' . "\n" : "\n")}
                     @before_script if (@before_script);
    @after_script  = map {($_ ? $_ . ' || exit 1' . "\n" : "\n")}
                     @after_script if (@after_script);

    # Write patch script to output
    my $out = File::Spec->catfile ($outdir_rev, 'apply-patch');
    open FILE, '>', $out or die $out, ': cannot open (', $!, ')';

    # Script header
    my $shell = Fcm::Config->instance()->setting(qw/TOOL SHELL/);
    print FILE <<EOF;
#!$shell
# ------------------------------------------------------------------------------
# NAME
#   apply-patch
#
# DESCRIPTION
#   This script is generated automatically by the "fcm mkpatch" command. It
#   applies the patch to the current working directory which must be a working
#   copy of a valid project tree that can accept the import of the patches.
#
#   Patch created from $organisation URL: $u
#   Changeset: $rev
# ------------------------------------------------------------------------------

this=`basename \$0`
echo "\$this: Applying patch for changeset $rev."

# Location of the patch, base on the location of this script
cd `dirname \$0` || exit 1
fcm_patch_dir=\$PWD

# Change directory back to the working copy
cd \$OLDPWD || exit 1

# Check working copy does not have local changes
status=`svn status`
if [[ -n \$status ]]; then
  echo "\$this: working copy contains changes, abort." >&2
  exit 1
fi
if [[ -a "#commit_message#" ]]; then
  echo "\$this: existing commit message in "#commit_message#", abort." >&2
  exit 1
fi

# Apply the changes
EOF

    # Script content
    print FILE @before_script if @before_script;
    print FILE "patch -p0 <\${fcm_patch_dir}/patchfile || exit 1\n"
      if $patch_needed;
    print FILE @after_script  if @after_script;

    # Script footer
    print FILE <<EOF;

# Copy in the commit message
cp \${fcm_patch_dir}/log-message "#commit_message#"

echo "\$this: finished normally."
#EOF
EOF

    close FILE or die $out, ': cannot close (', $!, ')';

    # Add executable permission
    chmod 0755, $out;

    # Script to commit the change
    push @script, '${fcm_patches_dir}/' . $rev . '/apply-patch';
    push @script, 'svn commit -F "#commit_message#"';
    push @script, 'rm -f "#commit_message#"';
    push @script, 'svn update';
    push @script, '';

    $CLI_MESSAGE->('PATCH_REV', $rev);
  }

  # Write the main output script if necessary. Otherwise remove output directory
  # ----------------------------------------------------------------------------
  if (@script) {
    # Add line break to each line in @script
    @script = map {($_ ? $_ . ' || exit 1' . "\n" : "\n")} @script;

    # Write script to output
    my $out = File::Spec->catfile ($outdir, 'fcm-import-patch');
    open FILE, '>', $out or die $out, ': cannot open (', $!, ')';

    # Script header
    my $shell = Fcm::Config->instance()->setting(qw/TOOL SHELL/);
    print FILE <<EOF;
#!$shell
# ------------------------------------------------------------------------------
# NAME
#   fcm-import-patch
#
# SYNOPSIS
#   fcm-import-patch TARGET
#
# DESCRIPTION
#   This script is generated automatically by the "fcm mkpatch" command, as are
#   the revision "patches" created in the same directory. The script imports the
#   patches into TARGET, which must either be a URL or a working copy of a valid
#   project tree that can accept the import of the patches.
#
#   Patch created from $organisation URL: $u
# ------------------------------------------------------------------------------

this=`basename \$0`

# Check argument
target=\$1

# First argument must be a URL or working copy
if [[ -z \$target ]]; then
  echo "\$this: the first argument must be a URL or a working copy, abort." >&2
  exit 1
fi

if [[ \$target == svn://*  || \$target == svn+ssh://* || \\
      \$target == http://* || \$target == https://*   || \\
      \$target == file://* ]]; then
  # A URL, checkout a working copy in a temporary location
  fcm_tmp_dir=`mktemp -d \${TMPDIR:=/tmp}/\$this.XXXXXX`
  fcm_working_copy=\$fcm_tmp_dir
  svn checkout -q \$target \$fcm_working_copy || exit 1
else
  fcm_working_copy=\$target
fi

# Location of the patches, base on the location of this script
cd `dirname \$0` || exit 1
fcm_patches_dir=\$PWD

# Change directory to the working copy
cd \$fcm_working_copy || exit 1

# Set the language to avoid encoding problems
export LANG=en_GB

# Commands to apply patches
EOF

    # Script content
    print FILE @script;

    # Script footer
    print FILE <<EOF;
# Remove temporary working copy, if necessary
if [[ -d \$fcm_tmp_dir && -w \$fcm_tmp_dir ]]; then
  rm -rf \$fcm_tmp_dir
fi

echo "\$this: finished normally."
#EOF
EOF

    close FILE or die $out, ': cannot close (', $!, ')';

    # Add executable permission
    chmod 0755, $out;

    # Diagnostic
    $CLI_MESSAGE->('PATCH_DONE', $outdir);

  } else {
    # Remove output directory
    rmtree $outdir or die $outdir, ': cannot remove';

    # Diagnostic
    return _cm_abort(Fcm::Cm::Abort->NULL);
  }

  return 1;
}

# ------------------------------------------------------------------------------
# CLI: fcm add.
sub _cli_command_add {
    my @args = map {($_ eq '--check' || $_ eq '-c' ? () : $_)} @_;
    my %option = (st_check_handler => $CLI_HANDLER_OF{'WC_STATUS_PATH'});
    return (
        @args == @_ ? _svn("add", @args) : cm_check_unknown(\%option, @args)
    );
}

# ------------------------------------------------------------------------------
# CLI: fcm checkout.
sub _cli_command_checkout {
    if (@ARGV) {
        my $target = is_url($ARGV[-1]) ? cwd() : $ARGV[-1];
        if (-d $target && is_wc($target)) {
            return _cm_err(Fcm::Cm::Exception->WC_EXIST, $target);
        }
    }
    return _svn('checkout', @ARGV);
}

# ------------------------------------------------------------------------------
# CLI: fcm delete.
sub _cli_command_delete {
    my @args = map {($_ eq '--check' || $_ eq '-c' ? () : $_)} @_;
    my %option = (st_check_handler => $CLI_HANDLER_OF{'WC_STATUS_PATH'});
    return (
        @args == @_ ? _svn("delete", @args) : cm_check_missing(\%option, @args)
    );
}

# ------------------------------------------------------------------------------
# CLI: fcm switch.
sub _cli_command_switch {
    local(@ARGV) = @_;
    if (grep {$_ eq '--relocate'} @ARGV) {
        return _svn('switch', @ARGV);
    }
    my %option;
    if (!GetOptions(\%option, 'non-interactive', 'revision|r=s', 'quiet|q')) {
        _cli_err();
    }
    if (!$option{'non-interactive'}) {
        $option{st_check_handler} = $CLI_HANDLER_OF{WC_STATUS};
    }
    if (!@ARGV) {
        _cli_err();
    }
    $CLI_MESSAGE->(q{}, join(q{}, cm_switch(\%option, @ARGV)));
}

# ------------------------------------------------------------------------------
# CLI: fcm update.
sub _cli_command_update {
    local(@ARGV) = @_;
    my %option;
    if (!GetOptions(\%option, 'non-interactive', 'revision|r=s', 'quiet|q')) {
        _cli_err();
    }
    if (!$option{'non-interactive'}) {
        $option{st_check_handler} = $CLI_HANDLER_OF{WC_STATUS};
    }
    $CLI_MESSAGE->(q{}, join(q{}, cm_update(\%option, @ARGV)));
}

# ------------------------------------------------------------------------------
# CLI error.
sub _cli_err {
    my ($key, @args) = @_;
    $key ||= 'CLI_USAGE';
    my $message = sprintf($CLI_MESSAGE_FOR_ERROR{$key}, @args);
    die(Fcm::CLI::Exception->new({message => $message}));
}

# ------------------------------------------------------------------------------
# Handles abort exception.
sub _cli_e_handler_of_cm_abort {
    my ($function, $e) = @_;
    if ($e->get_code() eq $e->FAIL) {
        die(sprintf($CLI_MESSAGE_FOR_ABORT{FAIL}, $function));
    }
    else {
        $CLI_MESSAGE->($e->get_code(), $function);
    }
}

# ------------------------------------------------------------------------------
# Handles CM exception.
sub _cli_e_handler_of_cm_exception {
    my ($function, $e) = @_;
    die(sprintf($CLI_MESSAGE_FOR_ERROR{$e->get_code()}, $e->get_targets()));
}

# ------------------------------------------------------------------------------
# Handles CLI exception.
sub _cli_e_handler_of_cli_exception {
    my ($function, $e) = @_;
    $CLI_MESSAGE->('CLI', $e);
    $CLI_MESSAGE->('CLI_HELP', $function);
}

# ------------------------------------------------------------------------------
# The default handler of the "WC_STATUS" event.
sub _cli_handler_of_wc_status {
    my ($name, $target_list_ref, $status_list_ref) = @_;
    if (@{$status_list_ref}) {
        $CLI_MESSAGE->('STATUS', join(q{}, @{$status_list_ref}));
        if ($CLI_PROMPT->($name, 'CONTINUE') ne 'y') {
            return _cm_abort();
        }
    }
    return @{$status_list_ref};
}

# ------------------------------------------------------------------------------
# The default handler of the "WC_STATUS_PATH" event.
sub _cli_handler_of_wc_status_path {
    my ($name, $target_list_ref, $status_list_ref) = @_;
    $CLI_MESSAGE->(q{}, join(q{}, @{$status_list_ref}));
    my @paths = map {chomp(); ($_ =~ $PATTERN_OF{ST_PATH})} @{$status_list_ref};
    my @paths_of_interest;
    while (my $path = shift(@paths)) {
        my %handler_of = (
            a => sub {push(@paths_of_interest, $path, @paths); @paths = ()},
            n => sub {},
            y => sub {push(@paths_of_interest, $path)},
        );
        my $reply = $CLI_PROMPT->(
            {type => 'yna'}, $name, 'RUN_SVN_COMMAND', "$name $path",
        );
        $handler_of{$reply}->();
    }
    return @paths_of_interest;
}

# ------------------------------------------------------------------------------
# Prints help for a given $subcommand.
sub _cli_help {
    my ($key, $exit_val) = @_;
    my $pod
        = File::Spec->catfile(dirname($INC{'Fcm/Cm.pm'}), 'CLI', "fcm-$key.pod");
    my $has_pod = -f $pod;
    if ($has_pod) {
        pod2usage({
            '-exitval' => defined($exit_val) ? $exit_val : 2,
            '-input'   => $pod,
            '-verbose' => 1,
        });
    }
    if (!$has_pod || exists($CLI_MORE_HELP_FOR{$key})) {
        local(@ARGV) = ($key);
        return _svn('help', $key);
    }
}

# ------------------------------------------------------------------------------
# Expands location keywords in a list.
sub _cli_keyword_expand_url {
    my ($arg_list_ref) = @_;
    ARG:
    for my $arg (@{$arg_list_ref}) {
        my ($label, $value) = ($arg =~ $PATTERN_OF{CLI_OPT});
        if (!$label) {
            ($label, $value) = (q{}, $arg);
        }
        if (!$value) {
            next ARG;
        }
        eval {
            $value = Fcm::Util::tidy_url(Fcm::Keyword::expand($value));
        };
        if ($@) {
            if ($value ne 'fcm:revision') {
                die($@);
            }
        }
        $arg = $label . $value;
    }
}

# ------------------------------------------------------------------------------
# Expands revision keywords in -r and --revision options in a list.
sub _cli_keyword_expand_rev {
    my ($arg_list_ref) = @_;
    my @targets;
    for my $arg (@{$arg_list_ref}) {
        if (-e $arg && is_wc($arg) || is_url($arg)) {
            push(@targets, $arg);
        }
    }
    if (!@targets) {
        push(@targets, get_url_of_wc());
    }
    if (!@targets) {
        return;
    }
    my @old_arg_list = @{$arg_list_ref};
    my @new_arg_list = ();
    ARG:
    while (defined(my $arg = shift(@old_arg_list))) {
        my ($key, $value) = $arg =~ $PATTERN_OF{CLI_OPT_REV};
        if (!$key) {
            push(@new_arg_list, $arg);
            next ARG;
        }
        push(@new_arg_list, '--revision');
        if (!$value) {
            $value = shift(@old_arg_list);
        }
        my @revs = grep {defined()} ($value =~ $PATTERN_OF{CLI_OPT_REV_RANGE});
        my ($url, @url_list) = @targets;
        for my $rev (@revs) {
            if ($rev !~ $PATTERN_OF{SVN_REV}) {
                $rev = (Fcm::Keyword::expand($url, $rev))[1];
            }
            if (@url_list) {
                $url = shift(@url_list);
            }
        }
        push(@new_arg_list, join(q{:}, @revs));
    }
    @{$arg_list_ref} = @new_arg_list;
}

# ------------------------------------------------------------------------------
# Prints a message.
sub _cli_message {
    my ($key, @args) = @_;
    for (
        [\*STDOUT, \%CLI_MESSAGE_FOR        , q{}          ],
        [\*STDERR, \%CLI_MESSAGE_FOR_WARNING, q{[WARNING] }],
        [\*STDERR, \%CLI_MESSAGE_FOR_ABORT  , q{[ABORT] }  ],
        [\*STDERR, \%CLI_MESSAGE_FOR_ERROR  , q{[ERROR] }  ],
    ) {
        my ($handle, $hash_ref, $prefix) = @{$_};
        if (exists($hash_ref->{$key})) {
            return printf({$handle} $prefix . $hash_ref->{$key}, @args);
        }
    }
}

# ------------------------------------------------------------------------------
# Wrapper for Fcm::Interactive::get_input.
sub _cli_prompt {
    my %option
        = (type => 'yn', default => 'n', (ref($_[0]) ? %{shift(@_)} : ()));
    my ($name, $key, @args) = @_;
    return Fcm::Interactive::get_input(
        title   => $CLI_PROMPT_PREFIX . $name,
        message => sprintf($CLI_MESSAGE_FOR_PROMPT{$key}, @args),
        %option,
    );
}

# ------------------------------------------------------------------------------
# Check missing status and delete.
sub cm_check_missing {
    my %option = %{shift()};
    my $checker
        = _svn_status_checker('delete', 'MISSING', $option{st_check_handler});
    my @paths = $checker->(\@_);
    if (@paths) {
        run_command([qw{svn delete}, @paths]);
    }
}

# ------------------------------------------------------------------------------
# Check unknown status and add.
sub cm_check_unknown {
    my %option = %{shift()};
    my $checker
        = _svn_status_checker('add', 'UNKNOWN', $option{st_check_handler});
    my @paths = $checker->(\@_);
    if (@paths) {
        run_command([qw{svn add}, @paths]);
    }
}

# ------------------------------------------------------------------------------
# FCM wrapper to SVN switch.
sub cm_switch {
    my %option = %{shift()};
    my ($target, $path) = @_;
    $path ||= cwd();
    if (!-e $path) {
        return _cm_err(Fcm::Cm::Exception->NOT_EXIST, $path);
    }
    if (!is_wc($path)) {
        return _cm_err(Fcm::Cm::Exception->INVALID_WC, $path);
    }

    # Check for merge template in the commit log file in the working copy
    my $path_of_wc = get_wct($path);
    my $ci_mesg = Fcm::CmCommitMessage->new();
    $ci_mesg->dir($path_of_wc);
    $ci_mesg->read_file();
    if (@{$ci_mesg->auto_mesg()}) {
        return _cm_err(
            Fcm::Cm::Exception->SWITCH_UNSAFE,
            $path eq $path_of_wc ? $ci_mesg->base() : $ci_mesg->file(),
        );
    }

    # Check for any local modifications
    if (defined($option{st_check_handler})) {
        my $handler = $CLI_HANDLER_OF{WC_STATUS};
        _svn_status_checker('switch', 'MODIFIED', $handler)->([$path_of_wc]);
    }

    # Invokes "svn switch"
    _svn(
        {METHOD => 'qx', PRINT => !$option{quiet}},
        'switch',
        ($option{'non-interactive'} ? '--non-interactive'       : ()),
        ($option{revision}          ? ('-r', $option{revision}) : ()),
        ($option{quiet}             ? '--quiet'                 : ()),
        _cm_get_source(
            $target,
            Fcm::CmBranch->new(URL => get_url_of_wc($path_of_wc)),
        )->url_peg(),
        ($path_of_wc eq cwd() ? () : $path_of_wc),
    );
}

# ------------------------------------------------------------------------------
# FCM wrapper to SVN update.
sub cm_update {
    my %option = %{shift()};
    my @targets = @_;
    if (!@targets) {
        @targets = (cwd());
    }
    for my $target (@targets) {
        if (!-e $target) {
            return _cm_err(Fcm::Cm::Exception->NOT_EXIST, $target);
        }
        if (!is_wc($target)) {
            return _cm_err(Fcm::Cm::Exception->INVALID_WC, $target);
        }
        $target = get_wct($target);
        if ($target eq cwd()) {
            $target = q{.};
        }
    }
    if (defined($option{st_check_handler})) {
        my ($matcher_keys_ref, $show_updates)
            = defined($option{revision}) ? (['MODIFIED'               ], undef)
            :                              (['MODIFIED', 'OUT_OF_DATE'], 1    )
            ;
        my $matcher = sub {
            for my $key (@{$matcher_keys_ref}) {
                $ST_MATCHER_FOR{$key}->(@_) && return 1;
            }
        };
        _svn_status_checker(
            'update', $matcher, $option{st_check_handler}, $show_updates,
        )->(\@targets);
    }
    if ($option{revision} && $option{revision} !~ $PATTERN_OF{SVN_REV}) {
        $option{revision} = (
            Fcm::Keyword::expand(get_url_of_wc($targets[0]), $option{revision})
        )[1];
    }
    return _svn_update(\@targets, \%option);
}

# ------------------------------------------------------------------------------
# Raises an abort exception.
sub _cm_abort {
    my ($code) = @_;
    $code ||= Fcm::Cm::Abort->USER;
    die(bless({code => $code, message => 'abort'}, 'Fcm::Cm::Abort'));
}

# ------------------------------------------------------------------------------
# Raises a failure.
sub _cm_err {
    my ($code, @targets) = @_;
    die(bless(
        {code => $code, message => "ERROR: $code", targets => \@targets},
        'Fcm::Cm::Exception',
    ));
}

# ------------------------------------------------------------------------------
# Returns the corresponding Fcm::CmBranch instance for $src_url w.r.t. $target.
sub _cm_get_source {
    my ($src_url, $target) = @_;
    my $source = Fcm::CmBranch->new(URL => $src_url);
    if (!$source->is_url()) {
        # Not a full URL, construct full URL based on current URL
        $source->url_peg($target->url_peg());
        my $project = $target->project();
        my ($path) = $src_url =~ qr{\A/*(.*)\z}xms;
        if (index($path, $project) == 0) {
            # Argument contains the full path under the repository root
            $path = substr($path, length($project));
        }
        if ($path =~ $PATTERN_OF{FCM_BRANCH_PATH}) {
            # Argument contains the full branch name
            $path = join(q{/}, $target->project_path(), $path);
        }
        else {
            # Argument contains the shorter branch name
            $path = join(q{/}, $target->project_path(), 'branches', $path);
        }
        $source->path_peg($path);
    }
    # Replace source sub-directory with the target sub-directory
    $source->subdir($target->subdir());
    # Ensure that the branch name exists
    if (!$source->url_exists()) {
        return _cm_err(Fcm::Cm::Exception->INVALID_URL, $src_url);
    }
    # Ensure that the branch name is valid
    if (!$source->branch()) {
        return _cm_err(Fcm::Cm::Exception->INVALID_BRANCH, $src_url);
    }
    # Ensure that the source and target URLs are in the same project
    if ($source->project_url() ne $target->project_url()) {
        return _cm_err(
            Fcm::Cm::Exception->DIFF_PROJECTS,
            $target->url_peg(),
            $source->url_peg(),
        );
    }
    return $source;
}

# ------------------------------------------------------------------------------
# Runs "svn".
sub _svn {
    my @args = @_;
    my %option;
    if (@args && ref($args[0])) {
        %option = %{shift(@args)};
    }
    return run_command(
        ['svn', @args],
        PRINT => ($args[0] ne 'cat' && !grep {$_ eq '--xml'} @args),
        %option,
    );
}

# ------------------------------------------------------------------------------
# Returns the results of "svn status".
sub _svn_status_get {
    my ($target_list_ref, $show_updates) = @_;
    my @targets = (defined($target_list_ref) ? @{$target_list_ref} : ());
    for my $target (@targets) {
        if ($target eq cwd()) {
            $target = q{.};
        }
    }
    my @options = ($show_updates ? qw{--show-updates} : ());
    return _svn({METHOD => 'qx', PRINT => 0}, 'status', @options, @targets);
}

# ------------------------------------------------------------------------------
# Returns a "svn status" checker.
sub _svn_status_checker {
    my ($name, $matcher, $handler, $show_updates) = @_;
    if (!ref($matcher)) {
        $matcher = $ST_MATCHER_FOR{$matcher};
    }
    return sub {
        my ($target_list_ref) = @_;
        my @status = _svn_status_get($target_list_ref, $show_updates);
        if ($show_updates) {
            @status = map {$_ =~ $PATTERN_OF{ST_AGAINST_REV} ? () : $_} @status;
        }
        my @status_of_interest = grep {$matcher->($_)} @status;
        if (defined($handler)) {
            return $handler->($name, $target_list_ref, \@status_of_interest);
        }
        return @status_of_interest;
    }
}

# ------------------------------------------------------------------------------
# Runs "svn update".
sub _svn_update {
    my ($target_list_ref, $option_hash_ref) = @_;
    my %option = (defined($option_hash_ref) ? %{$option_hash_ref} : ());
    _svn(
        {METHOD => 'qx', PRINT => !$option{quiet}},
        'update',
        ($option{'non-interactive'} ? '--non-interactive'       : ()),
        ($option{revision}          ? ('-r', $option{revision}) : ()),
        ($option{quiet}             ? '--quiet'                 : ()),
        (defined($target_list_ref) ? @{$target_list_ref} : ()),
    );
}

# ------------------------------------------------------------------------------
# Abort exception.
package Fcm::Cm::Abort;
use base qw{Fcm::Exception};
use constant {FAIL => 'FAIL', NULL => 'NULL', USER => 'USER'};

sub get_code {
    return $_[0]->{code};
}

# ------------------------------------------------------------------------------
# Resource exception.
package Fcm::Cm::Exception;
our @ISA = qw{Fcm::Cm::Abort};
use constant {
    CHDIR             => 'CHDIR',
    INVALID_BRANCH    => 'INVALID_BRANCH',
    INVALID_PROJECT   => 'INVALID_PROJECT',
    INVALID_TARGET    => 'INVALID_TARGET',
    INVALID_URL       => 'INVALID_URL',
    INVALID_WC        => 'INVALID_WC',
    MERGE_REV_INVALID => 'MERGE_REV_INVALID',
    MERGE_SELF        => 'MERGE_SELF',
    MERGE_UNRELATED   => 'MERGE_UNRELATED',
    MERGE_UNSAFE      => 'MERGE_UNSAFE',
    MKPATH            => 'MKPATH',
    NOT_EXIST         => 'NOT_EXIST',
    PARENT_NOT_EXIST  => 'PARENT_NOT_EXIST',
    RMTREE            => 'RMTREE',
    SWITCH_UNSAFE     => 'SWITCH_UNSAFE',
    WC_EXIST          => 'WC_EXIST',
    WC_INVALID_BRANCH => 'WC_INVALID_BRANCH',
    WC_URL_NOT_EXIST  => 'WC_URL_NOT_EXIST',
};

sub get_targets {
    return @{$_[0]->{targets}};
}

1;
__END__

=pod

=head1 NAME

Fcm::Cm

=head1 SYNOPSIS

    use Fcm::Cm qw{cli};

    # Use as a wrapper to Subversion, and other FCM code management commands
    cli('info', '--revision', 'HEAD', $url);

    use Fcm::Cm qw{cm_check_missing cm_check_unknown cm_switch cm_update};

    # Checks status for "missing" items and "svn delete" them
    $missing_st_handler = sub {
        my ($name, $target_list_ref, $status_list_ref) = @_;
        # ...
        return @paths_of_interest;
    };
    cm_check_missing({st_check_handler => $missing_st_handler}, @targets);

    # Checks status for "unknown" items and "svn add" them
    $unknown_st_handler = sub {
        my ($name, $target_list_ref, $status_list_ref) = @_;
        # ...
        return @paths_of_interest;
    };
    cm_check_unknown({st_check_handler => $unknown_st_handler}, @targets);

    # Sets up a status checker
    $st_check_handler = sub {
        my ($name, $target_list_ref, $status_list_ref) = @_;
        # ...
    };
    # Switches a "working copy" at the "root" level to a new URL target
    cm_switch(
        {
            'non-interactive'  => $non_interactive_flag,
            'quiet'            => $quiet_flag,
            'revision'         => $revision,
            'st_check_handler' => $st_check_handler,
        },
        $target, $path_of_wc,
    );
    # Runs "svn update" on each working copy from their "root" level
    cm_update(
        {
            'non-interactive'  => $non_interactive_flag,
            'quiet'            => $quiet_flag,
            'revision'         => $revision,
            'st_check_handler' => $st_check_handler,
        },
        @targets,
    );

=head1 DESCRIPTION

Wraps the Subversion client and implements other FCM code management
functionalities.

=head1 FUNCTIONS

=over 4

=item cli($function,@args)

Implements the FCM code management CLI. If --help or -h is specified in @args,
it displays help and returns.  Otherwise, it attempts to expand any FCM location
and revision keywords in @args. Calls the relevant FCM code management function
according to $function, or a SVN command if $function is not modified by FCM.

=item cm_check_missing(\%option,@targets)

Use "svn status" to check for missing items in @targets. If @targets is an empty
list, the function adds the current working directory to it. Expects
$option{st_check_handler} to be a CODE reference. Calls
$option{st_check_handler} with ($name, $target_list_ref, $status_list_ref) where
$name is "delete", $target_list_ref is \@targets, and $status_list_ref is an
ARRAY reference to a list of "svn status" output with the "missing" status.
$option{st_check_handler} should return a list of interesting paths, which will
be scheduled for removal using "svn delete".

=item cm_check_unknown(\%option,@targets)

Similar to cm_check_missing(\%option,@targets) but checks for "unknown" items,
which will be scheduled for addition using "svn add".

=item cm_switch(\%option,$target,$path_of_wc)

Invokes "svn switch" at the root of a working copy specified by $path_of_wc (or
the current working directory if $path_of_wc is not specified).
$option{'non-interactive'}, $option{quiet}, $option{revision} determines the
options (of the same name) that are passed to "svn switch". If
$option{st_check_handler} is set, it should be a CODE reference, and will be
called with ('switch', [$path_of_wc], $status_list_ref), where $status_list_ref
is an ARRAY reference to the output returned by "svn status" on $path_of_wc.
This can be used for the application to display the working copy status to the
user before prompting him/her to continue. The return value of
$option{st_check_handler} is ignored.

=item cm_update(\%option,@targets)

Invokes "svn update" at the root of each working copy specified by @targets. If
@targets is an empty list, the function adds the current working directory to
it. $option{'non-interactive'}, $option{quiet}, $option{revision} determines the
options (of the same name) that are passed to "svn update". If
$option{st_check_handler} is set, it should be a CODE reference, and will be
called with ($name, $target_list_ref, $status_list_ref), where $name is
'update', $target_list_ref is \@targets and $status_list_ref is an ARRAY
reference to the output returned by "svn status -u" on the @targets. This can be
used for the application to display the working copy update status to the user
before prompting him/her to continue. The return value of
$option{st_check_handler} is ignored.

=back

=head1 DIAGNOSTICS

The following exceptions can be raised:

=over 4

=item Fcm::Cm::Abort

This exception @ISA L<Fcm::Exception|Fcm::Exception>. It is raised if a command
is aborted for some reason. The $e->get_code() method can be used to retrieve an
error code, which can be one of the following:

=over 4

=item $e->FAIL

The command aborts because of a failure.

=item $e->NULL

The command aborts because it will result in no change.

=item $e->USER

The command aborts because of an action by the user.

=back

=item Fcm::Cm::Exception

This exception @ISA L<Fcm::Abort|Fcm::Abort>. It is raised if a command fails
with a known reason. The $e->get_targets() method can be used to retrieve a list
of targets/resources associated with this exception. The $e->get_code() method
can be used to retrieve an error code, which can be one of the following:

=over 4

=item $e->CHDIR

Fails to change directory to a target.

=item $e->INVALID_BRANCH

A target is not a valid branch URL in the standard FCM project layout.

=item $e->INVALID_PROJECT

A target is not a valid project URL in the standard FCM project layout.

=item $e->INVALID_TARGET

A target is not a valid Subversion URL or working copy.

=item $e->INVALID_URL

A target is not a valid Subversion URL.

=item $e->INVALID_WC

A target is not a valid Subversion working copy.

=item $e->MERGE_REV_INVALID

An invalid revision (target element 0) is specified for a merge.

=item $e->MERGE_SELF

Attempt to merge a URL (target element 0) to its own working copy (target
element 1).

=item $e->MERGE_UNRELATED

The merge target (target element 0) is not directly related to the merge source
(target element 1).

=item $e->MERGE_UNSAFE

A merge source (target element 0) contains changes outside the target
sub-directory.

=item $e->MKPATH

Fail to create a directory (target element 0) recursively.

=item $e->NOT_EXIST

A target does not exist.

=item $e->PARENT_NOT_EXIST

The parent of the target no longer exists.

=item $e->RMTREE

Fail to remove a directory (target element 0) recursively.

=item $e->SWITCH_UNSAFE

A merge template exists in the commit message file (target element 0) in a
working copy target.

=item $e->WC_EXIST

The target working copy already exists.

=item $e->WC_INVALID_BRANCH

The URL of the target working copy is not a valid branch URL in the standard FCM
project layout.

=item $e->WC_URL_NOT_EXIST

The URL of the target working copy no longer exists at the HEAD revision.

=back

=back

=head1 TO DO

Reintegrate with L<Fcm::CmUrl|Fcm::CmUrl> and L<Fcm::CmBranch|Fcm::CmBranch>,
but separate this module into the CLI part and the CM part. Expose the remaining
CM functions when this is done.

Use L<SVN::Client|SVN::Client> to interface with Subversion.

Move C<mkpatch> out of this module.

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
