bedtools intersect -a exon.bed -b snp.bed -c -s | awk '{if($7>=2) print;}' | cut -f1-6 | bedtools sort -i stdin > exons_with_2_snps.bed
