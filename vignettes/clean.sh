#!/usr/bin/env bash

for f in $(ls -I "*.Rmd" -I "*.sh"); do
  echo "... removing '${f}'"
  rm -rf ${f}
done
