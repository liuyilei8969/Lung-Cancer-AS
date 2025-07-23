#!/bin/bash

cat upevent.txt | while read line; do  # maybe dnevent.txt
grep ${line} Hs.seq.all.cass.chrom.can.exon.bed >> DSposition.txt
done
