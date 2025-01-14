#!/bin/bash

cat DSname.txt | while read line; do
grep ${line} Hs.seq.all.cass.chrom.can.exon.bed >> DSposition.txt
done
