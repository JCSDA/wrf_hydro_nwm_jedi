#!/bin/bash
cd $1
touch .repo_hash_log
mv .repo_hash_log .repo_hash_log_previous
echo `date` > .repo_hash_log
for repo in fckit atlas oops saber ioda ioda-converters ufo wrf_hydro_nwm_jedi; do cd $repo; echo `pwd`  `git rev-parse HEAD`; cd ..; done >> .repo_hash_log
echo >> .repo_hash_log
head -n1000 .repo_hash_log_previous >> .repo_hash_log
exit 0
