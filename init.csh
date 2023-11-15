#!/bin/csh
source dnastorage.env

git rev-parse --git-dir >& /dev/null

set is_git_status = $?
if (${is_git_status} == 0) then
    git submodule update --init --recursive
else
    git clone https://github.com/kvolkel/schwimmbad
endif

if ( "$argv[1]" == "-no-env" ) then
    exit
endif
