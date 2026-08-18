#!/bin/bash
# install domclust
cd package/
git clone https://github.com/mbgd-prj/domclust.git
cd domclust
./configure
make

