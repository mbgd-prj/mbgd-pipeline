#!/bin/bash
# install database
mkdir -p db
cd db
current_ver=`curl https://mbgd.nibb.ac.jp/dist/current`
curl -O https://mbgd.nibb.ac.jp/dist/${current_ver}/mbgd_prof.tgz 
tar xvfz mbgd_prof.tgz
