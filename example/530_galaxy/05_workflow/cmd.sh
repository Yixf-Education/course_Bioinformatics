bedtools intersect -a chrY_hg19_exon.bed -b chrY_hg19_snp.bed -c -s | awk '{if($7>=10) print;}' | cut -f1-6 | bedtools sort -i stdin > exons_with_10_snps_with_strand.bed

bedtools intersect -a chrY_hg19_exon.bed -b chrY_hg19_snp.bed -c | awk '{if($7>=10) print;}' | cut -f1-6 | bedtools sort -i stdin > exons_with_10_snps_without_strand.bed
