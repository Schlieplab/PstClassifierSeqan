#!/usr/bin/bash
if [ -d /home/ghal/git-data/tmp ];
  then
  rm -rf /home/ghal/git-data/tmp/*
  cd     /home/ghal/git-data/tmp
else
  mkdir /home/ghal/git-data/tmp
  cd    /home/ghal/git-data/tmp
fi
c=0
csplit -s -z "$1" '/>/' '{*}'
for i in xx* ; do \
  if grep -q chromosome "$i";  then
    n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
    mv "$i" "$c.fa" ; \
    c=$((c+1))
  else
    rm -rf "$i";
  fi
done
