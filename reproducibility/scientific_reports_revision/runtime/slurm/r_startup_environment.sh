#!/usr/bin/env bash
# Project-process scope only. No personal/site R configuration file is changed.
: "${TSR_R_LIB:?Set the explicit isolated project library}"
[[ $TSR_R_LIB == /* ]] || { echo 'TSR_R_LIB must be absolute' >&2;return 2; }
# Prevent inherited/site/user startup files from rewriting package search paths,
# including in ordinary child Rscript calls that do not pass --vanilla.
export R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS="$TSR_R_LIB" R_LIBS_USER="$TSR_R_LIB" R_LIBS_SITE=NULL
