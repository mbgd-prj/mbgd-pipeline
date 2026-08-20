#!/bin/bash
# install domclust
MPL_TOPDIR=`pwd`
mkdir -p etc
echo "export MPL_TOPDIR=$MPL_TOPDIR" > etc/bashrc
echo "PATH=\$MPL_TOPDIR/bin:\$PATH" >> etc/bashrc
cd package/
if [ ! -f domclust/bin/domclust ]; then
	git clone https://github.com/mbgd-prj/domclust.git
	cd domclust
	./configure --prefix=$MPL_TOPDIR
	make  install
	cd ..
fi
if [ ! -f CGB/lib/CompareMap.jar ]; then
	wget https://mbgd.nibb.ac.jp/CGB/dist/CGB.tgz
	tar xvf CGB.tgz
	rm CGB.tgz
	cp CGB/lib/CompareMap.jar ../lib
fi
