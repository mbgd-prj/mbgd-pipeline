#!/bin/bash
# install database
mkdir -p db
cd db
latest_rel=`curl https://mbgd.nibb.ac.jp/dist/latest_rel`
curl -O https://mbgd.nibb.ac.jp/dist/${latest_rel}/${latest_rel}_prof.tgz 
tar xvfz mbgd_prof.tgz
