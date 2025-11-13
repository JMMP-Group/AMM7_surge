# ------------------------------------------------------------------------------
# NAME
#   Fcm::CmUrl
#
# DESCRIPTION
#   This class contains methods for manipulating a Subversion URL in a standard
#   FCM project.
#
# COPYRIGHT
#   (C) Crown copyright Met Office. All rights reserved.
#   For further details please refer to the file COPYRIGHT.txt
#   which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

package Fcm::CmUrl;
@ISA = qw(Fcm::Base);

# Standard pragma
use warnings;
use strict;

# Standard modules
use HTTP::Date;
use XML::DOM;

# FCM component modules
use Fcm::Base;
use Fcm::Keyword;
use Fcm::Util qw/run_command svn_date/;

# Special branches
our %owner_keywords = (Share => 'shared', Config => 'config', Rel => 'release');

# Revision pattern
my $rev_pattern = '\d+|HEAD|BASE|COMMITTED|PREV|\{.+\}';

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $cm_url = Fcm::CmUrl->new ([URL => $url,]);
#
# DESCRIPTION
#   This method constructs a new instance of the Fcm::CmUrl class.
#
# ARGUMENTS
#   URL - URL of a branch
# ------------------------------------------------------------------------------

sub new {
  my $this  = shift;
  my %args  = @_;
  my $class = ref $this || $this;

  my $self = Fcm::Base->new (%args);

  $self->{URL} = (exists $args{URL} ? $args{URL} : '');

  for (qw/ANALYSED BRANCH BRANCH_LIST INFO LIST LOG LOG_RANGE PEGREV RLIST
          PROJECT SUBDIR/) {
    $self->{$_} = undef;
  }

  bless $self, $class;
  return $self;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->url_peg;
#   $cm_url->url_peg ($url);
#
# DESCRIPTION
#   This method returns/sets the current URL@PEG.
# ------------------------------------------------------------------------------

sub url_peg {
  my $self = shift;

  if (@_) {
    if (! $self->{URL} or $_[0] ne $self->{URL}) {
      # Re-set URL
      $self->{URL} = shift;

      # Re-set essential variables
      $self->{$_}  = undef for (qw/ANALYSED RLIST LIST INFO LOG LOG_RANGE/);
    }
  }

  return $self->{URL};
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $flag = $cm_url->is_url ();
#
# DESCRIPTION
#   Returns true if current url is a valid Subversion URL.
# ------------------------------------------------------------------------------

sub is_url {
  my $self = shift;

  # This should handle URL beginning with svn://, http:// and svn+ssh://
  return ($self->url_peg =~ m#^[\+\w]+://#);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $flag = $cm_url->url_exists ([$rev]);
#
# DESCRIPTION
#   Returns true if current url exists (at operative revision $rev) in a
#   Subversion repository.
# ------------------------------------------------------------------------------

sub url_exists {
  my ($self, $rev) = @_;

  my $exists = $self->svnlist (REV => $rev);

  return defined ($exists);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $string = $cm_url->svninfo ([FLAG => $flag], [REV => $rev]);
#
# DESCRIPTION
#   Returns the value of $flag, where $flag is a field returned by "svn info".
#   (If $flag is not set, default to "URL".) Otherwise returns an empty string.
#   If REV is specified, it will be used as the operative revision.
# ------------------------------------------------------------------------------

sub svninfo {
  my $self = shift;
  my %args = @_;

  my $flag = exists $args{FLAG} ? $args{FLAG} : 'URL';
  my $rev  = exists $args{REV}  ? $args{REV}  : undef;

  $rev = ($self->pegrev ? $self->pegrev : 'HEAD') if not $rev;

  return if not $self->is_url;

  # Get "info" for the specified revision if necessary
  if (not exists $self->{INFO}{$rev}) {
    # Invoke "svn info" command
    my @info = &run_command (
      [qw/svn info -r/, $rev, $self->url_peg],
      PRINT   => $self->config->verbose > 2,
      METHOD  => 'qx',
      DEVNULL => 1,
      ERROR   => 'ignore',
    );

    # Store selected information
    for (@info) {
      chomp;

      if (/^(.+?):\s*(.+)$/) {
        $self->{INFO}{$rev}{$1} = $2;
      }
    }
  }

  my $return = exists $self->{INFO}{$rev}{$flag}
               ? $self->{INFO}{$rev}{$flag} : undef;

  return $return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   %logs = $cm_url->svnlog (
#     [REV          => $rev,]
#     [REV          => \@revs,] # reference to a 2-element array
#     [STOP_ON_COPY => 1,]
#   );
#
# DESCRIPTION
#   Returns the logs for the current URL. If REV is a range of revisions or not
#   specified, return a hash where the keys are revision numbers and the values
#   are the entries (which are hash references). If a single REV is specified,
#   return the entry (a hash reference) at the specified REV. Each entry in the
#   returned list is a hash reference, with the following structure:
#
#   $entry = {
#     author => $author,              # the commit author
#     date   => $date,                # the commit date (in seconds since epoch)
#     msg    => $msg,                 # the log message
#     paths  => {                     # list of changed paths
#       $path1  => {                  # a changed path
#         copyfrom-path => $frompath, # copy-from-path
#         copyfrom-rev  => $fromrev,  # copy-from-revision
#         action        => $action,   # action status code
#       },
#       ...     => { ... },           # ... more changed paths ...
#     },
#   }
# ------------------------------------------------------------------------------

sub svnlog {
  my $self = shift;
  my %args = @_;

  my $stop_on_copy  = exists $args{STOP_ON_COPY} ? $args{STOP_ON_COPY} : 0;
  my $rev_arg       = exists $args{REV}          ? $args{REV}          : 0;

  my @revs;

  # Get revision options
  # ----------------------------------------------------------------------------
  if ($rev_arg) {
    if (ref ($rev_arg)) {
      # Revsion option is an array, a range of revisions specified?
      ($revs [0], $revs [1]) = @$rev_arg;

    } else {
      # A single revision specified
      $revs [0] = $rev_arg;
    }

    # Expand 'HEAD' revision
    for my $rev (@revs) {
      next unless uc ($rev) eq 'HEAD';
      $rev = $self->svninfo (FLAG => 'Revision', REV => 'HEAD');
    }

  } else {
    # No revision option specified, get log for all revisions
    $revs [0] = $self->svninfo (FLAG => 'Revision');
    $revs [1] = 1;
  }

  $revs [1] = $revs [0] if not $revs [1];
  @revs     = sort {$b <=> $a} @revs;

  # Check whether a "svn log" run is necessary
  # ----------------------------------------------------------------------------
  my $need_update = ! ($revs [0] == $revs [1] and exists $self->{LOG}{$revs [0]});
  my @ranges      = @revs;
  if ($need_update and $self->{LOG_RANGE}) {
    my %log_range = %{ $self->{LOG_RANGE} };

    if ($stop_on_copy) {
      $ranges [1] = $log_range{UPPER} if $ranges [1] >= $log_range{LOWER_SOC};

    } else {
      $ranges [1] = $log_range{UPPER} if $ranges [1] >= $log_range{LOWER};
    }
  }

  $need_update = 0 if $ranges [0] < $ranges [1];

  if ($need_update) {
    # Invoke "svn log" command for all revisions of the current branch
    # --------------------------------------------------------------------------
    my @command = (
      qw/svn log --xml -v/, ($stop_on_copy ? '--stop-on-copy' : ()),
      '-r' . join (':', @ranges),
      $self->url_peg,
    );

    my $rc;
    my @xml = &run_command (
      \@command,
      PRINT   => $self->config->verbose > 2,
      METHOD  => 'qx',
      DEVNULL => 1,
      ERROR   => 'ignore',
      RC      => \$rc,
    );

    # Parse the XML
    # --------------------------------------------------------------------------
    if (not $rc) {
      my $parser = XML::DOM::Parser->new;
      my $doc    = $parser->parse (join ('', @xml));

      my $entry_list = $doc->getElementsByTagName ('logentry');

      # Record the author, date, message and path change for each revision
      for my $i (0 .. $entry_list->getLength - 1) {
        # Select current entry from node list
        my $entry = $entry_list->item ($i);
        my %this = ();

        # Revision is an attribute of the entry node
        my $rev   = $entry->getAttributeNode ('revision')->getValue;

        # Author, date and log message are children elements of the entry node
        for my $key (qw/author date msg/) {
          # Get data of each node, also convert date to seconds since epoch
          my $node    = $entry->getElementsByTagName ($key)->item (0);
          my $data    = ($node and $node->getFirstChild)
                        ? $node->getFirstChild->getData : '';
          $this{$key} = ($key eq 'date' ? str2time ($data) : $data);
        }

        # Path nodes are grand children elements of the entry node
        my $paths = $entry->getElementsByTagName ('path');

        for my $p (0 .. $paths->getLength - 1) {
          # Select current path node from node list
          my $node = $paths->item ($p);

          # Get data from the path node
          my $path = $node->getFirstChild->getData;
          $this{paths}{$path} = {};

          # Action, copyfrom-path and copyfrom-rev are attributes of path nodes
          for my $key (qw/action copyfrom-path copyfrom-rev/) {
            next unless $node->getAttributeNode ($key); # ensure attribute exists

            $this{paths}{$path}{$key} = $node->getAttributeNode ($key)->getValue;
          }
        }

        $self->{LOG}{$rev} = \%this;
      }
    }

    # Update the range cache
    # --------------------------------------------------------------------------
    # Upper end of the range
    $self->{LOG_RANGE}{UPPER} = $ranges [0]
      if ! $self->{LOG_RANGE}{UPPER} or $ranges [0] > $self->{LOG_RANGE}{UPPER};

    # Lower end of the range, need to take into account the stop-on-copy option
    if ($stop_on_copy) {
      # Lower end of the range with stop-on-copy option
      $self->{LOG_RANGE}{LOWER_SOC} = $ranges [1]
        if ! $self->{LOG_RANGE}{LOWER_SOC} or
           $ranges [1] < $self->{LOG_RANGE}{LOWER_SOC};

      my $low = (sort {$a <=> $b} keys %{ $self->{LOG} }) [0];
      $self->{LOG_RANGE}{LOWER} = $low
        if ! $self->{LOG_RANGE}{LOWER} or $low < $self->{LOG_RANGE}{LOWER};

    } else {
      # Lower end of the range without the stop-on-copy option
      $self->{LOG_RANGE}{LOWER} = $ranges [1]
        if ! $self->{LOG_RANGE}{LOWER} or
           $ranges [1] < $self->{LOG_RANGE}{LOWER};

      $self->{LOG_RANGE}{LOWER_SOC} = $ranges [1]
        if ! $self->{LOG_RANGE}{LOWER_SOC} or
           $ranges [1] < $self->{LOG_RANGE}{LOWER_SOC};
    }
  }

  my %return = ();

  if (! $rev_arg or ref ($rev_arg)) {
    # REV is an array, return log entries if they are within range
    for my $rev (sort {$b <=> $a} keys %{ $self->{LOG} }) {
      next if $rev > $revs [0] or $revs [1] > $rev;

      $return{$rev} = $self->{LOG}{$rev};

      if ($stop_on_copy) {
        last if exists $self->{LOG}{$rev}{paths}{$self->branch_path} and
           $self->{LOG}{$rev}{paths}{$self->branch_path}{action} eq 'A';
      }
    }

  } else {
    # REV is a scalar, return log of the specified revision if it exists
    %return = %{ $self->{LOG}{$revs [0]} } if exists $self->{LOG}{$revs [0]};
  }

  return %return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $string = $cm_branch->display_svnlog ($rev, [$wiki]);
#
# DESCRIPTION
#   This method returns a string for displaying the log of the current branch
#   at a $rev. If $wiki is set, returns a string for displaying in a Trac wiki
#   table.  The value of $wiki should be the Subversion URL of a FCM project
#   associated with the intended Trac system.
# ------------------------------------------------------------------------------

sub display_svnlog {
  my ($self, $rev, $wiki) = @_;
  my $return = '';

  my %log = $self->svnlog (REV => $rev);

  if ($wiki) {
    # Output in Trac wiki format
    # --------------------------------------------------------------------------
    $return .= '|| ' . &svn_date ($log{date}) . ' || ' . $log{author} . ' || ';

    my $trac_url = Fcm::Keyword::get_browser_url($self->url);

    # Get list of tickets from log
    my @tickets;
    while ($log{msg} =~ /(?:(\w+):)?(?:#|ticket:)(\d+)/g) {
      push @tickets, [$1, $2];
    }
    @tickets = sort {
      if ($a->[0] and $b->[0]) {
        $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1];

      } elsif ($a->[0]) {
        1;

      } else {
        $a->[1] <=> $b->[1];
      }
    } @tickets;

    if ($trac_url =~ m#^$wiki(?:/*|$)#) {
      # URL is in the specified $wiki, use Trac link
      $return .= '[' . $rev . '] ||';

      for my $ticket (@tickets) {
        $return .= ' ';
        $return .= $ticket->[0] . ':' if $ticket->[0];
        $return .= '#' . $ticket->[1];
      }

      $return .= ' ||';

    } else {
      # URL is not in the specified $wiki, use full URL
      my $rev_url = $trac_url;
      $rev_url    =~ s{/intertrac/source:.*\z}{/intertrac/changeset:$rev}xms;
      $return    .= '[' . $rev_url . ' ' . $rev . '] ||';

      my $ticket_url = $trac_url;
      $ticket_url    =~ s{/intertrac/source:.*\z}{/intertrac/}xms;

      for my $ticket (@tickets) {
        $return .= ' [' . $ticket_url;
        $return .= $ticket->[0] . ':' if $ticket->[0];
        $return .= 'ticket:' . $ticket->[1] . ' ' . $ticket->[1] . ']';
      }

      $return .= ' ||';
    }

  } else {
    # Output in plain text format
    # --------------------------------------------------------------------------
    my @msg  = split /\n/, $log{msg};
    my $line = (@msg > 1 ? ' lines' : ' line');

    $return .= join (
      ' | ',
      ('r' . $rev, $log{author}, &svn_date ($log{date}), scalar (@msg) . $line),
    );
    $return .= "\n\n";
    $return .= $log{msg};
  }

  return $return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   @list = $cm_url->svnlist ([REV => $rev], [RECURSIVE => 1]);
#
# DESCRIPTION
#   The method returns a list of paths as returned by "svn list". If RECURSIVE
#   is set, "svn list" is invoked with the "-R" option.
# ------------------------------------------------------------------------------

sub svnlist {
  my $self = shift;
  my %args = @_;

  my $recursive = exists $args{RECURSIVE} ? $args{RECURSIVE} : 0;
  my $rev       = exists $args{REV}       ? $args{REV}       : undef;
  my $key       = $recursive ? 'RLIST' : 'LIST';

  # Find out last changed revision of the current URL
  $rev = $self->svninfo (FLAG => 'Last Changed Rev', REV => $rev);
  return () if not $rev;

  # Get directory listing for the current URL at the last changed revision
  if (not exists $self->{$key}{$rev}) {
    my $rc;

    my @list = map {chomp; $_} &run_command (
      [qw/svn list -r/, $rev, ($recursive ? '-R' : ()), $self->url_peg],
      METHOD => 'qx', ERROR => 'ignore', DEVNULL => 1, RC => \$rc,
    );

    $self->{$key}{$rev} = $rc ? undef : \@list;
  }

  return (defined ($self->{$key}{$rev}) ? @{ $self->{$key}{$rev} } : undef);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   @list = $cm_url->branch_list ($rev);
#
# DESCRIPTION
#   The method returns a list of branches in the current project, assuming the
#   FCM naming convention. If $rev if specified, it returns the list of
#   branches at $rev.
# ------------------------------------------------------------------------------

sub branch_list {
  my ($self, $rev) = @_;

  # Current URL must be a valid FCM project
  return if not $self->project;

  # Find out last changed revision of the current URL
  $rev = $self->svninfo (FLAG => 'Revision', REV => $rev);
  return () if not $rev;

  if (not exists $self->{BRANCH_LIST}{$rev}) {
    $self->{BRANCH_LIST}{$rev} = [];

    # Get URL of the project "branches/" sub-directory
    my $url = Fcm::CmUrl->new (URL => $self->project_url . '/branches');

    # List three levels underneath "branches/"
    # First level, i.e. dev, test, pkg, etc
    my @list1 = map {$url->url . '/' . $_} $url->svnlist (REV => $rev);
    @list1    = grep m#/$#, @list1;

    # Second level, i.e. user name, Shared, Rel or Config
    my @list2;
    for (@list1) {
      my $u    = Fcm::CmUrl->new (URL => $_);
      my @list = $u->svnlist (REV => $rev);

      push @list2, map {$u->url . $_} @list;
    }

    # Third level, branch name
    for (@list2) {
      my $u    = Fcm::CmUrl->new (URL => $_);
      my @list = map {s#/*$##; $_} $u->svnlist (REV => $rev);

      push @{ $self->{BRANCH_LIST}{$rev} }, map {$u->url . $_} @list;
    }
  }

  return @{ $self->{BRANCH_LIST}{$rev} };
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $self->_analyse_url ();
#
# DESCRIPTION
#   The method analyses the current URL, breaking it up into the project
#   (substring of URL up to the slash before "trunk", "branches" or "tags"),
#   branch name ("trunk", "branches/<type>/<id>/<name>" or "tags/<name>") and
#   the sub-directory below the top of the project sub-tree. It re-sets the
#   corresponding interal variables.
# ------------------------------------------------------------------------------

sub _analyse_url {
  my $self = shift;
  my ($url, $project, $branch, $subdir, $pegrev);

  # Check that URL is set
  $url    = $self->url_peg;
  return if not $url;
  return if not $self->is_url;

  # Extract from URL the peg revision
  $pegrev = $1 if $url =~ s/@($rev_pattern)$//i;

  if ($url =~ m#^(.*?)/+(trunk|branches|tags)(?:/+(.*))?/*$#) {
    # URL is under the "trunk", a branch or a tag
    $project                 = $1;
    my ($branch_id, $remain) = ($2, $3);

    $remain = '' if not defined $remain;

    if ($branch_id eq 'trunk') {
      # URL under the "trunk"
      $branch = 'trunk';

    } else {
      # URL under a branch or a tag
      $branch = $branch_id;

      # Assume "3 sub-directories", FCM branch naming convention
      for (1 .. 3) {
        if ($remain =~ s#^([^/]+)(?:/+|$)##) {
          $branch .=  '/' . $1;

        } else {
          $branch = undef;
          last;
        }
      }
    }

    $subdir = $remain ? $remain : '' if $branch;

  } else {
    # URL is at some level above the "trunk", a branch or a tag
    # Use "svn ls" to determine whether it is a project URL
    my @list = $self->svnlist (REV => ($pegrev ? $pegrev : 'HEAD'));
    my %lines = map {chomp $_; ($_, 1)} @list;

    # A project URL should have the "trunk", "branches" and "tags" directories
    ($project = $url) =~ s#/*$##
      if $lines{'trunk/'} and $lines{'branches/'} and $lines{'tags/'};
  }

  $self->{PROJECT}  = $project;
  $self->{BRANCH}   = $branch;
  $self->{SUBDIR}   = $subdir;
  $self->{PEGREV}   = $pegrev;
  $self->{ANALYSED} = 1;

  return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->root ();
#
# DESCRIPTION
#   The method returns the repository root of the current URL.
# ------------------------------------------------------------------------------

sub root {
  my $self = shift;

  return $self->svninfo (FLAG => 'Repository Root');
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->project_url_peg ();
#   $cm_url->project_url_peg ($url);
#
# DESCRIPTION
#   The method returns the URL@PEG of the "project" part of the current URL. If
#   an argument is specified, the URL of the "project" part and the peg
#   revision of the current URL are re-set.
# ------------------------------------------------------------------------------

sub project_url_peg {
  my $self = shift;

  if (@_) {
    my $url = shift;

    # Re-construct URL is necessary
    if (! $self->project_url_peg or $url ne $self->project_url_peg) {
      my $pegrev = ($url =~ s/@($rev_pattern)$//i) ? $1 : '';

      $url .= '/' . $self->branch if $self->branch;
      $url .= '/' . $self->subdir if $self->subdir;
      $url .= '@' . $pegrev       if $pegrev;

      $self->url_peg ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->{PROJECT} . ($self->pegrev ? '@' . $self->pegrev : '');
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->project_url ();
#   $cm_url->project_url ($url);
#
# DESCRIPTION
#   The method returns the URL of the "project" part of the current URL. If an
#   argument is specified, the URL of the "project" part of the current URL is
#   re-set.
# ------------------------------------------------------------------------------

sub project_url {
  my $self = shift;

  if (@_) {
    my $url = shift;
    $url =~ s/@($rev_pattern)$//i;

    # Re-construct URL is necessary
    if (! $self->project_url or $url ne $self->project_url) {
      $url .= '/' . $self->branch if $self->branch;
      $url .= '/' . $self->subdir if $self->subdir;

      $self->url ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->{PROJECT};
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $path = $cm_url->project_path ();
#   $cm_url->project_path ($path);
#
# DESCRIPTION
#   The method returns the path of the "project" part of the current URL. If an
#   argument is specified, the path of the "project" part of the current URL is
#   re-set.
# ------------------------------------------------------------------------------

sub project_path {
  my $self = shift;

  # Repository root
  my $root = $self->root;
  $root    = substr (
    $self->project_url,
    0,
    length ($self->project_url) - length ($self->project) - 1
  ) if not $root;

  if (@_) {
    my $path = shift;

    # Re-construct URL is necessary
    if (! $self->project_path or $path ne $self->project_path) {
      $path .= '/' . $self->branch if $self->branch;
      $path .= '/' . $self->subdir if $self->subdir;

      $self->path ($path);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return substr ($self->{PROJECT}, length ($root));
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $name = $cm_url->project ();
#   $cm_url->project ($name);
#
# DESCRIPTION
#   The method returns the basename of the "project" part of the current URL.
#   If an argument is specified, the basename of the "project" part of the
#   current URL is re-set.
# ------------------------------------------------------------------------------

sub project {
  my $self = shift;

  if (@_) {
    my $name = shift;

    # Re-construct URL is necessary
    if (! $self->project or $name ne $self->project) {
      my $url = '';
      if ($self->project) {
        $url =  $self->project;
        $url =~ s#/[^/]+$##;

      } else {
        $url =  $self->root;
      }

      $url .=  '/' . $name;
      $url .=  '/' . $self->branch if $self->branch;
      $url .=  '/' . $self->subdir if $self->subdir;
      $url .=  '@' . $self->pegrev if $self->pegrev;

      $self->url_peg ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  my $name = $self->{PROJECT};
  $name =~ s#^.*/([^/]+)$#$1# if $name;

  return $name;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->branch_url_peg ();
#   $cm_url->branch_url_peg ($url);
#
# DESCRIPTION
#   The method returns the URL@PEG of the "branch" part of the current URL. If
#   an argument is specified, the URL@PEG of the "branch" part of the current
#   URL is re-set.
# ------------------------------------------------------------------------------

sub branch_url_peg {
  my $self = shift;

  if (@_) {
    my $url = shift;

    # Re-construct URL is necessary
    if (! $self->branch_url_peg or $url ne $self->branch_url_peg) {
      my $pegrev = ($url =~ s/@($rev_pattern)$//i) ? $1 : '';

      $url .= '/' . $self->subdir if $self->subdir;
      $url .= '@' . $pegrev       if $pegrev;

      $self->url_peg ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->project_url . '/' . $self->branch .
         ($self->pegrev ? '@' . $self->pegrev : '');
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->branch_url ();
#   $cm_url->branch_url ($url);
#
# DESCRIPTION
#   The method returns the URL of the "branch" part of the current URL. If an
#   argument is specified, the URL of the "branch" part of the current URL is
#   re-set.
# ------------------------------------------------------------------------------

sub branch_url {
  my $self = shift;

  if (@_) {
    my $url = shift;
    $url =~ s/@($rev_pattern)$//i;

    # Re-construct URL is necessary
    if (! $self->branch_url or $url ne $self->branch_url) {
      $url .= '/' . $self->subdir if $self->subdir;

      $self->url ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->project_url . '/' . $self->branch;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $path = $cm_url->branch_path ();
#   $cm_url->branch_path ($path);
#
# DESCRIPTION
#   The method returns the path of the "branch" part of the current URL. If an
#   argument is specified, the path of the "branch" part of the current URL is
#   re-set.
# ------------------------------------------------------------------------------

sub branch_path {
  my $self = shift;

  if (@_) {
    my $path = shift;

    # Re-construct URL is necessary
    if (! $self->branch_path or $path ne $self->branch_path) {
      $path .= '/' . $self->subdir if $self->subdir;

      $self->path ($path);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return ($self->branch ? $self->project_path . '/' . $self->branch : undef);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $branch = $cm_url->branch ();
#   $cm_url->branch ($branch);
#
# DESCRIPTION
#   The method returns the "branch" part of the current URL. If an argument is
#   specified, the "branch" part of the current URL is re-set.
# ------------------------------------------------------------------------------

sub branch {
  my $self = shift;

  if (@_) {
    my $branch = shift;

    # Re-construct URL is necessary
    if (! $self->branch or $branch ne $self->branch) {
      my $url = $self->project_url;
      $url   .= '/' . $branch;
      $url   .= '/' . $self->subdir if $self->subdir;

      $self->url ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->{BRANCH};
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $string = $obj->branch_owner;
#
# DESCRIPTION
#   This method returns the owner of the branch.
# ------------------------------------------------------------------------------

sub branch_owner {
  my $self = shift;
  my $return;

  if ($self->is_branch and $self->branch_url =~ m#/([^/]+)/[^/]+/*$#) {
    my $user = $1;
    $return = $user;
  }

  return $return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $flag = $cm_url->is_trunk ();
#
# DESCRIPTION
#   The method returns true if the the current URL is (a sub-tree of) the trunk.
# ------------------------------------------------------------------------------

sub is_trunk {
  my $self = shift;

  $self->_analyse_url () if not $self->{ANALYSED};

  return ($self->branch and $self->branch eq 'trunk');
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $flag = $cm_url->is_branch ();
#
# DESCRIPTION
#   The method returns true if the the current URL is (a sub-tree of) a branch.
# ------------------------------------------------------------------------------

sub is_branch {
  my $self = shift;

  $self->_analyse_url () if not $self->{ANALYSED};

  return ($self->branch and $self->branch =~ m#^branches/#);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $flag = $cm_url->is_tag ();
#
# DESCRIPTION
#   The method returns true if the the current URL is (a sub-tree of) a tag.
# ------------------------------------------------------------------------------

sub is_tag {
  my $self = shift;

  $self->_analyse_url () if not $self->{ANALYSED};

  return ($self->branch and $self->branch =~ m#^tags/#);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $subdir = $cm_url->subdir ();
#   $cm_url->subdir ($subdir);
#
# DESCRIPTION
#   The method returns the "subdir" part of the current URL. If an argument is
#   specified, the "subdir" part of the current URL is re-set.
# ------------------------------------------------------------------------------

sub subdir {
  my $self = shift;

  if (@_) {
    my $subdir = shift;

    # Re-construct URL is necessary
    if (! $self->subdir or $subdir ne $self->subdir) {
      my $url = $self->project_url;
      $url   .= '/' . $self->branch if $self->branch;
      $url   .= '/' . $subdir if $subdir;

      $self->url ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->{SUBDIR};
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $url = $cm_url->url ();
#   $cm_url->url ($url);
#
# DESCRIPTION
#   The method returns the URL without the "peg revision" part. If an argument
#   is specified, the URL is re-set without modifying the "peg revision" part.
# ------------------------------------------------------------------------------

sub url {
  my $self = shift;

  if (@_) {
    my $url = shift;
    $url    =~ s/@($rev_pattern)$//i;

    # Re-construct URL if necessary
    if (! $self->url or $url ne $self->url) {
      $self->url_peg ($url . ($self->pegrev ? '@' . $self->pegrev : ''));
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  (my $url = $self->url_peg) =~ s/@($rev_pattern)$//i;

  return $url;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $path = $cm_url->path ();
#   $cm_url->path ($path);
#
# DESCRIPTION
#   The method returns the "path" part of the URL (i.e. URL without the
#   "root" part). If an argument is specified, the "path" part of the URL is
#   re-set.
# ------------------------------------------------------------------------------

sub path {
  my $self = shift;

  # Repository root
  my $root = $self->root;
  $root    = substr (
    $self->project_url,
    0,
    length ($self->project_url) - length ($self->project) - 1
  ) if not $root;

  if (@_) {
    my $path = shift;
    $path    =~ s/@($rev_pattern)$//i;

    # Re-construct URL is necessary
    if (! $self->path or $path ne $self->path) {
      my $url = ($root . (substr ($path, 0, 1) eq '/' ? '' : '/') . $path);
      $self->url ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return substr ($self->url, length ($root));
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $path = $cm_url->path_peg ();
#   $cm_url->path_peg ($path);
#
# DESCRIPTION
#   The method returns the PATH@PEG part of the URL (i.e. URL without the
#   "root" part). If an argument is specified, the PATH@PEG part of the URL is
#   re-set.
# ------------------------------------------------------------------------------

sub path_peg {
  my $self = shift;

  # Repository root
  my $root = $self->root;
  $root    = substr (
    $self->project_url,
    0,
    length ($self->project_url) - length ($self->project) - 1
  ) if not $root;

  if (@_) {
    my $path = shift;

    # Re-construct URL is necessary
    if (! $self->path_peg or $path ne $self->path_peg) {
      my $url = ($root . (substr ($path, 0, 1) eq '/' ? '' : '/') . $path);
      $self->url_peg ($url);
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return substr ($self->url_peg, length ($root));
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $rev = $cm_url->pegrev ();
#   $cm_url->pegrev ($rev);
#
# DESCRIPTION
#   The method returns the "peg revision" part of the current URL. If an
#   argument is specified, the "peg revision" part of the current URL is
#   re-set.
# ------------------------------------------------------------------------------

sub pegrev {
  my $self = shift;

  if (@_) {
    my $pegrev = shift;

    # Re-construct URL is necessary
    if (! $self->pegrev or $pegrev ne $self->pegrev) {
      $self->url_peg ($self->url . ($pegrev ? '@' . $pegrev : ''));
    }
  }

  $self->_analyse_url () if not $self->{ANALYSED};

  return $self->{PEGREV};
}

# ------------------------------------------------------------------------------

1;

__END__
