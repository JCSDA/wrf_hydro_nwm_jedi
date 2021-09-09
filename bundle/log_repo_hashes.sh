#!/bin/bash
cd $1
touch .repo_hash_log
mv .repo_hash_log .repo_hash_log_previous
echo `date` > .repo_hash_log

# check if repo directory exists before adding to repo_hash_log
check_repos=(fckit atlas oops saber ioda ioda-converters ufo wrf_hydro_nwm_jedi)
declare -a repos
for repo in "${check_repos[@]}"; do
    [ -d $repo ] && repos=(${repos[@]} "$repo")
done

for repo in "${repos[@]}"; do cd $repo; echo `pwd` `git rev-parse HEAD`; cd ..; done >> .repo_hash_log
echo >> .repo_hash_log
head -n1000 .repo_hash_log_previous >> .repo_hash_log

exit 0
