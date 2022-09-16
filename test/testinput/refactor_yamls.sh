#!/bin/bash

for file in *.yaml *.yml
do
  echo "******** $file *********"
  src=$file
  dest=YAML_UPGRADE/$file
  echo "$src --> $dest"
  python3 /Users/afox/Jedi/ioda_internal/ioda-bundle/ioda/tools/refactor-yaml.py -i $src -o $dest
done
