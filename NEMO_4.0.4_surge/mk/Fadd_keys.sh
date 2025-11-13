#!/bin/bash
######################################################
# Author : Simona Flavoni for NEMO
# Contact : sflod@locean-ipsl.upmc.fr
#
# Some functions called from makenemo
# Fadd_keys   : add keys in cpp.fcm file  
######################################################
#set -x
set -o posix
#set -u
#set -e
#+
#
# ============
# Fadd_keys.sh
# ============
#
# --------------------
# Add compilation keys
# --------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fadd_keys.sh CONFIG_NAME add_key "LIST_KEYS"
#
#
# DESCRIPTION
# ===========
#
#
# Script to add a set of key when compiling a configuration.
# The list of key to be added has to be enclosed with " ". 
# A 'sed' is performed to modify the CONFIG_NAME/cpp.fcm file to    
# add the new key(s). 
#
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fadd_keys.sh ORCA2_LIM add_key "key_mpp_rep"
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
# $Id: Fadd_keys.sh 2158 2010-10-20 17:30:03Z sflod $
#
#
#
#   * creation
#
#-
 echo "Adding keys in : ${NEW_CONF}" 
 for i in ${list_add_key} ; do
	if [ "$(cat ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm | grep -c "\<$i\>" )" -ne 0 ] ; then
 	   echo "key $i already present in cpp_${NEW_CONF}.fcm" 
	else
	   sed -e "s/$/ ${i}/"  ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm >  ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm.tmp
 	   mv ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm.tmp   ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm
 	   echo "added key $i in ${NEW_CONF}" 
	fi
 done
 
 unset -v list_add_key

