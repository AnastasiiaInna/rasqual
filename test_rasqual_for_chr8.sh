#!/bin/sh

#  test_rasqual_for_chr8.sh
#  
#
#  Created by Hryhorzhevska, Anastasiia on 16.12.20.
#
VCF_FILENAME=~/github/rasqual/data/signe_chipseq/VCF_DATA/phased/unimputed_subset/phased_rasual.vcf.gz
COUNT_MTRX_FILENAME=~/github/rasqual/data/signe_chipseq/count_mtrx_dex_for_rasqual.bin
OFSSET_MTRX_FILENAME=~/github/rasqual/data/signe_chipseq/dex_offset.bin
PEAKSET_INFO_FILENAME=~/github/rasqual/data/signe_chipseq/peakset_info_dex_sample_range.txt

CHR=8
CIS_REG_START=11000000
CIS_REG_END=12000000
GENE_START=11245707
GENE_END=11846107

PEAK_START=11445707#11545707
PEAK_END=11546107#11546107

PEAK=$CHR:$PEAK_START:$PEAK_END

L=$(tabix $VCF_FILENAME $CHR:$CIS_REG_START-$CIS_REG_END | wc -l | awk '{print $1}')
M=$(tabix $VCF_FILENAME $CHR:$GENE_START-$GENE_END | wc -l | awk '{print $1}')

tabix $VCF_FIILENAME | bin/rasqual -y $COUNT_MTRX_FILENAME -k $OFSSET_MTRX_FILENAME -n 48 -j 1  -l $L -M $M -s $PEAK_START -e $PEAK_END -t -f $PEAK -z --no-posterior-update
