#Generate RPKM and PSI matrix
/path/to/quantas/countit/gen_expression_matrix.pl -base /path/to/expression/data -v expression.conf expression_matrix.txt
/path/to/quantas/quantas/countit/gen_splicing_matrix.pl -base /path/to/AS/data --id2gene2symbol Hs.seq.all.AS.chrom.can.id2gene2symbol -v AS.conf splicing_matrix_cass.txt

#Generate DE DS result
perl /path/to/quantas/countit/test_expr_diff.pl -base  /path/to/expression/data -rpkm -logFC -v expressiondiff.conf expression_diff.txt
perl /path/to/quantas/countit/test_splicing_diff.pl -base /path/to/AS/data --id2gene2symbol Hs.seq.all.AS.chrom.can.id2gene2symbol -v ASdiff.conf splicing_cass_diff.txt

#Analyse DS DE result
awk '$3 <= -1 && $5 <= -0.05 && $9 == 1 {print $2}' splicing_cass_diff.txt > expressiondn.txt
awk '$3 >= 1 && $5 <= -0.05 && $9 == 1 {print $2}' splicing_cass_diff.txt > expressionup.txt
awk '$10 >= 20 && $13 <= -0.1 && $15 <= 0.01 {print $5}' splicing_cass_diff.txt > ASdn.txt
awk '$10 >= 20 && $13 >= 0.1 && $15 <= 0.01 {print $5}' splicing_cass_diff.txt > ASup.txt
