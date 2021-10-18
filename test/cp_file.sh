#!/bin/bash

# Copy first argument into further ones
for ((i = 2; i <= $#; i++ )); do
  cp $1 ${!i}
done

exit $?
