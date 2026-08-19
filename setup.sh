#!/bin/bash
# install domclust
MPL_TOPDIR=`pwd`
echo "export MPL_TOPDIR=$MPL_TOPDIR" > etc/bashrc
echo "PATH=\$MPL_TOPDIR/bin:\$PATH" >> etc/bashrc
cd package/
if [ ! -f domclust/bin/domclust ]; then
	git clone https://github.com/mbgd-prj/domclust.git
	cd domclust
	./configure --prefix=$MPL_TOPDIR
	make  install
fi
