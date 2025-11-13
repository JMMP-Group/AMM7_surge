#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
#
# ============
# Fread_dir.sh
# ============
#
# ---------------------
# Read user directories
# ---------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fread_dir.sh
#
#
# DESCRIPTION
# ===========
#
#
# Read directoires needed from standard input
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fread_dir.sh Directory_NAME YES/NO
#
#
# TODO
# ====
#
# option debug
#
#
# EVOLUTIONS
# ==========
#
# $Id: Fread_dir.sh 3294 2012-01-28 16:44:18Z rblod $
#
#
#
#   * creation
#
#-

if [ "$2" == "YES" ]; then
   echo -n " $1 [Y/n] "
   read answer
   answer=`echo $answer | sed 's/^[yY].*$/y/'`
   if [  -z "$answer" -o "x$answer" = "xy" ]; then
      TAB[$ind]="$1"
      let ind=ind+1
      echo " $1 selected "
      echo "  "
   else
      echo " $1 Not selected "
      echo "  "
   fi
   unset -v answer
else
   echo -n " $1 [y/N] "
   read answer
   answer=`echo $answer | sed 's/^[nN].*$/N/'`
   if [ "x$answer" = "xy" ]; then
      TAB[$ind]="$1"
      let ind=ind+1
      echo " $1 selected "
      echo "  "
   else
      echo " $1 Not selected "
      echo "  "
   fi
   unset -v answer
fi
