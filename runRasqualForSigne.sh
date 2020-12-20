#!/bin/bash

if [ ${#@} -ne 5 ]; then
	printf "\nUsage: %s Not enough arguments\n" "${0}" > /dev/stderr;
	exit 1;
fi

VCF_FILENAME="${1}"
COUNT_MTRX_FILENAME="${2}"
OFSSET_MTRX_FILENAME="${3}"
PEAKSET_INFO_FILENAME="${4}"

NUM_TOTAL_SNP=480680
NUM_FEATURE_SNP=480680

RESULT_FILENAME="${5}"

if [ ! -f "${PEAKSET_INFO_FILENAME}" ] ; then
	printf 'Error: the file "%s" does not exist.\n' "${PEAKSET_INFO_FILENAME}" > /dev/stderr;
        exit 1;
fi

if [ -f "${RESULT_FILENAME}" ] ; then
	rm "${RESULT_FILENAME}" 
fi

touch "${RESULT_FILENAME}"

echo Feature_ID rs_ID Chromosome SNP_position Ref Alt Allele_freq HWE_Chi-square Imputation_quality_score Log_10_Benjamini-Hochberg_Q-value Chi_square Effect_size_Pi Sequencing_error_rate Reference_allele_mapping_bias_Phi Overdispersion SNP_ID_within_the_region No_of_feature_SNPs No_of_tested_SNPs No_of_iterations_for_null_hypothesis No_of_iterations_for_alternative_hypothesis Random_location_of_ties Log_likelihood_of_the_null_hypothesis Convergence_status Squared_correlation_between_prior_and_posterior_genotypes_fSNPs Squared_correlation_between_prior_and_posterior_genotypes_rSNP > $RESULT_FILENAME


jj=0
for PEAK in $(cat "${PEAKSET_INFO_FILENAME}") ; do
	jj=$((jj+1))
	CHR=$(awk -F: '{print $1}' <<< "${PEAK}");
	
	PEAK_START=$(awk -F: '{print $2}' <<< "${PEAK}");	
	PEAK_END=$(awk -F: '{print $3}' <<< "${PEAK}");
	
	GENE_START=$PEAK_START;
	GENE_END=$PEAK_END;
	
	CIS_REG_START=$((PEAK_START-5000));
	CIS_REG_END=$((PEAK_END+5000));
	
	L=$(tabix $VCF_FILENAME $CHR:$CIS_REG_START-$CIS_REG_END | wc -l | awk '{print $1}');
	M=$(tabix $VCF_FILENAME $CHR:$GENE_START-$GENE_END | wc -l | awk '{print $1}');
	
	# RASQUAL_RESULT=$(tabix "${VCF_FIILENAME}" | bin/rasqual -y "${COUNT_MTRX_FILENAME}" -k "${OFSSET_MTRX_FILENAME}" -n 47 -j jj  -l "${NUM_TOTAL_SNP}" -m "${NUM_FEATURE_SNP}" -s "${PEAK_START}" -e "${PEAK_END}" -t -f "${PEAK}" -z);
	RASQUAL_RESULT=$(tabix $VCF_FILENAME $CHR:$CIS_REG_START-$CIS_REG_END | bin/rasqual -y $COUNT_MTRX_FILENAME -k $OFSSET_MTRX_FILENAME -n 47 -j jj  -l $L -m $M -s $PEAK_START -e $PEAK_END -t -f $PEAK -z --no-posterior-update);	
	echo "${RASQUAL_RESULT}" >> $RESULT_FILENAME;
	
	printf "Peak start is %s and "${PEAK_START}";
	printf "Peak end is %s \n"${PEAK_END}";
done
