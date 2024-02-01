#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
#
# =============
# Fclean_var.sh
# =============
#
# ----------------------------
# Clean environment variables
# ----------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fclean_var.sh
#
#
# DESCRIPTION
# ===========
#
#
# Clean environment variables
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fclean_var.sh
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
# $Id: Fclean_var.sh 4990 2014-12-15 16:42:49Z timgraham $
#
#
#
#   * creation
#
#-
#- Unset variables

unset -v NSTOP
unset -v TAB
unset -v NEW_CONF
unset -v REF_CONF
unset -v CMP_NAM
unset -v NBR_PRC
unset -v NEM_SUBDIR
unset -v MAIN_DIR
unset -v CONFIG_DIR
unset -v TOOLS_DIR
unset -v COMPIL_DIR
unset -v NEMO_DIR
unset -v USEBLD
