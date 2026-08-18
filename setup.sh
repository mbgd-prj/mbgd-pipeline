#!/bin/bash
# install domclust
MPL_TOPDIR=`pwd`
echo "export \MPL_TOPDIR=$MPL_TOPDIR" > etc/bashrc
cd package/
git clone https://github.com/mbgd-prj/domclust.git
cd domclust
./configure
make
