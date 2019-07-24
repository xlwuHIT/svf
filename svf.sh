#python /home/xwu/1sv/svf.py feature -b /data/xwu/bam/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -c /data/xwu/gsd/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz  -v /home/xwu/opsv/NA12878low.vcf.gz -o /home/xwu/1sv/NA12878low_feature.txt -g hg19 -f /home/xwu/reference/hs37d5/hs37d5.fa
python svf.py feature -b /data/xwu/bam/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam -c /data/xwu/gsd/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz  -v /home/xwu/opsv/NA12878high.vcf.gz -o NA12878high_feature.txt -g hg19 -f /home/xwu/reference/hs37d5/hs37d5.fa

