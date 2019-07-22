from SV import *

gsd='/data/xwu/gsd/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz'
lump='/home/xwu/opsv/lumpy/NA12878high_lumpy.vcf.gz'
dell='/home/xwu/opsv/delly/NA12878high_lumpy.vcf.gz'

sv=SV(gsd,'NA12878','hg19','/home/xwu/reference/hs37d5/hs37d5.fa',True)
lsv=SV(lump,'NA12878','hg19','/home/xwu/reference/hs37d5/hs37d5.fa',True)
dsv=SV(dell,'NA12878','hg19','/home/xwu/reference/hs37d5/hs37d5.fa',True)
