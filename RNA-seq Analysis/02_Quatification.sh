#Alignment
olego -o sample_1.sam -j hg19.intron.hmr.bed -t 4 --strand-mode 3 -r olegodir/models/hg.cfg -v hg19.fa /path/to/sample/sample_1.fastq.gz >& /path/to/log/sample_1.log
olego -o sample_2.sam -j hg19.intron.hmr.bed -t 4 --strand-mode 3 -r olegodir/models/hg.cfg -v hg19.fa /path/to/sample/sample_2.fastq.gz >& /path/to/log/sample_2.log

#Generate Bed Files
perl /path/to/quantas/gapless/gapless_huge_file.pl --library-type unstranded -v -sam -uniq --split-size 10000000 -isoform hg19.exon.trio.hmr.nr.bed -E 400 -big --print-singleton -o ${f} ${f}_1.sam ${f}_2.sam

#Quatify Expression
perl /path/to/quantas/countit/summarize_expression_wrapper.pl -weight -big -exon hg19.exon.uniq.core.bed -e2g hg19.exon.uniq.core.id2gene2symbol -v ${f}/pair.gapless.bed ${f}.expr.txt

#Quatify Splicing
perl /path/to/quantas/countit/summarize_splicing_wrapper.pl -weight -big -conf hg19.conf -dbkey hg19 -cass -iret -mutx -taca -alt5 -alt3 -v ${f}/pair.gapless.bed ${f}
