#!/usr/bin/env bash

for f in $(ls .); do
  ext="${f##*.}"
  if [ "Rmd" != "${ext}" -a "sh" != "${ext}" ]; then
    cmd="rm -rfv ${f}"
    echo "$cmd"
    eval $cmd
  fi
done
