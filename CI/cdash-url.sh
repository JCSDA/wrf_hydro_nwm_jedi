#!/bin/bash

dir=$1
tag=$(head -1 $dir/TAG)
Done=$(cat $dir/$tag/Done.xml)
buildID=$(echo $Done | grep -o -P '(?<=buildId>).*(?=</build)')
URL=http://cdash.jcsda.org:8080/viewTest.php?buildid=$buildID
echo $URL
