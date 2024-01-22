#!/bin/bash
######################################################
# Author : Rachid Benshila for NEMO
# Contact : rblod@locean-ipsl.upmc.fr
#
# Some functions called from makenemo
# Fcheck_config   : config checking 
######################################################
#set -x
set -o posix
#set -u
#set -e
#+
#
# ================
# Fcheck_config.sh
# ================
#
# --------------------------
# Check the configuration
# --------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fcheck_config.sh FILENAME CONFNAME
#
#
# DESCRIPTION
# ===========
#
#
# Check the choice of the configuration:
#
# - Two cases
# - One is explicitely set
# - Nothing set, use the previous in use 
#
# We use TOOLS/cfgs_DIR/cfg.txt to check if the configuration exists.
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fcheck_config.sh
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
# $Id: Fcheck_config.sh 9598 2018-05-15 22:47:16Z nicolasmartin $
#
#
#
#   * creation
#
#-

declare -a ZTAB
if [ ${#2} -eq 0 ]; then
   tail -1  ${CONFIG_DIR}/$1	> ${CONFIG_DIR}/cfg.tmp
   read -a ZTAB < ${CONFIG_DIR}/cfg.tmp
   NEW_CONF=${ZTAB[0]} ; TAB=( ${ZTAB[@]:1} )
   \rm ${CONFIG_DIR}/cfg.tmp
   echo "Warning !!!"
   echo "No configuration specified"
   echo "Use makenemo -n MYCONFIG"
   echo "or  makenemo -h for help"
   echo "Using default configuration : ${NEW_CONF}"
fi
if [ "$1" == cfg.txt ]; then
   cat ${CONFIG_DIR}/$1 | grep "${NEW_CONF} " > ${CONFIG_DIR}/cfg.tmp
   read -a ZTAB < ${CONFIG_DIR}/cfg.tmp
   NEW_CONF=${ZTAB[0]} ; TAB=( ${ZTAB[@]:1} )
   \rm ${CONFIG_DIR}/cfg.tmp
fi

unset -v ZTAB
