#!/bin/bash
######################################################
# Author : Simona Flavoni for NEMO
# Contact : sflod@locean-ipsl.upmc.fr
#
# Some functions called from makenemo
# Fdel_keys   : del keys in cpp.fcm file  
######################################################
#set -x
set -o posix
#set -u
#set -e
#+
#
# ================
# Fdel_keys.sh
# ================
#
# --------------------------
# Add compilation keys
# --------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fdel_keys.sh CONFIG_NAME del_key "LIST_KEYS"
#
#
# DESCRIPTION
# ===========
#
#
# Add cpp keys when compiling a configuration, key list has to be enclosed with " ".
# We perform a 'sed' on the CONFIG_NAME/CPP.fcm file, containing the list of keys. 
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fdel_keys.sh CONFIG_NAME del_key "key_agrif"
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
# $Id: Fdel_keys.sh 2158 2010-10-20 17:30:03Z sflod $
#
#
#
#   * creation
#
#-

echo "Removing keys in : ${NEW_CONF}"

for i in ${list_del_key} ; do

     if [ "$(echo ${i} | grep -c key_nproc )" -ne 0                                      ]; then
        sed -e "s/key_nproc[ij]=.* //" ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm \
	    >  ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm.tmp
        mv ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm.tmp ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm
        echo " "
     elif [ "$(cat ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm | grep -c "$i" )" -ne 0 ]; then
         sed -e "s/\b${i}\b//" ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm \
	     >  ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm.tmp
         mv ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm.tmp ${CONFIG_DIR}/${NEW_CONF}/cpp_${NEW_CONF}.fcm
         echo "deleted key $i in ${NEW_CONF}"
     fi

done

unset -v list_del_key
