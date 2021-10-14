#!/bin/bash
# Theo Tricou


# data from barlow 2018

wget  https://datadryad.org/stash/downloads/file_stream/13237 \
      https://datadryad.org/stash/downloads/file_stream/13233 \
      https://datadryad.org/stash/downloads/file_stream/13224 \
      https://datadryad.org/stash/downloads/file_stream/13225 \
      https://datadryad.org/stash/downloads/file_stream/13226

gunzip *.fa.gz


# Script from barlow 2018

git clone git@github.com:jacahill/Admixture.git # installation instruction at https://github.com/jacahill/Admixture

/Admixture/D_stat Adm1_rep1_all.fa 235_rep1_all.fa \
  191Y_rep1_all.fa Uamericanus_all.fa 1000000 > ghost_bear.txt

python2 Admixture/D-stat_parser.py ghost_bear.txt 1 > D_ghost_bear.txt

python2 Admixture/weighted_block_jackknife.py ghost_bear.txt 1000000 > Err_ghost_bear.txt

# Simple commande to compute the Z-score from the D-statistic and the error

calc(){ awk "BEGIN { print "$*" }"; }

D=`awk '{print $2}' D_ghost_bear.txt`
E=`awk '{print $2}' Err_ghost_bear.txt`
calc $D/$E > Z_ghost_bear


/Admixture/D_stat Adm1_rep1_all.fa 235_rep1_all.fa \
  NB_rep1_all.fa Uamericanus_all.fa 1000000 > ingroup_bear.txt

python2 ~/GitHub/Admixture/D-stat_parser.py ingroup_bear.txt 1 > D_ingroup_bear.txt

python2 ~/GitHub/Admixture/weighted_block_jackknife.py ingroup_bear.txt 1000000 > Err_ingroup_bear.txt


D=`awk '{print $2}' D_ingroup_bear.txt`
E=`awk '{print $2}' Err_ingroup_bear.txt`
calc $D/$E > Z_ingroup_bear
