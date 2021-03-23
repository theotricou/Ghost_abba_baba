#!/bin/bash
# Theo Tricou


# download at https://datadryad.org/stash/dataset/doi:10.5061/dryad.cr1496b

for i in gz;do
  gzip -d $i
done

for i in *fa; do
  B="$(cut -d'_' -f2 <<<"$i")"
  mkdir $B
  mv $i $B/
done

# cut file in scaffolds
for i in */;do
  cd $i
  csplit *.fa '/scaffold/' '{*}'
  for j in xx*; do
    sed -i "s/scaffold.*/$i/g" $j
  done
  cd ..
done

# concate ali by species
mkdir polar_ali
cd NB
for i in xx*; do
  cat ../Adm1/$i ../235/$i ../NB/$i ../Uamericanus/$i > ../polar_ali/$i # selecet 4 species
done
cd ../polar_ali

rm xx00
file=
for i in  xx*[0-9]; do seaview -convert -sites  $i; done # correct for unaligned extremities
rm *fa
for i in *fst; do
  if [ "$i" == "xx01.fst" ]; then # skip first msa to concate with second msa
    cp $i $i\_n
    file=$i\_n
  else # concate msa with previous reste and cut in windows of 1000000
    seaview -o $i\_n -concatenate $i $file
    msa_split -w 1000000,0 $i\_n --out-root $i 2> log
    rm $file
    # file=$i\_n
    file=`tail -n 2 log | grep -ho "xx.*fa"`
  fi
done
rm *_n *[0-9] log

for i in *fa;do # correcte the space in > NAME and remove non-snp sites
  sed -i "s/> />/" $i
  snp-sites -c $i
done
rm *fa
