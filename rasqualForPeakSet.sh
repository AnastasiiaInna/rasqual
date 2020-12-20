#!/bin/bash

if [ ${#@} -ne 4 ]; then
	printf "\nUsage: %s Not enough arguments\n" "${0}" > /dev/stderr;
	exit 1;
fi

VCF_FIILENAME="${1}"
COUNT_MTRX_FILENAME="${2}"
OFSSET_MTRX_FILENAME="${3}"
PEAKSET_INFO_FILENAME="${4}"

NUM_TOTAL_SNP=480680
NUM_FEATURE_SNP=480680

RESULT_FILENAME=rasqual_result.txt

if [ ! -f "${PEAKSET_INFO_FILENAME}" ] ; then
	printf 'Error: the file "%s" does not exist.\n' "${PEAKSET_INFO_FILENAME}" > /dev/stderr;
        exit 1;
fi

touch "${RESULT_FILENAME}"

jj=0
for PEAK in $(cat "${PEAKSET_INFO_FILENAME}") ; do
	jj=$((jj+1))
	PEAK_START=$(awk -F: '{print $2}' <<< "${PEAK}");	
	PEAK_END=$(awk -F: '{print $3}' <<< "${PEAK}");

	RASQUAL_RESULT=$(tabix "${VCF_FIILENAME}" | bin/rasqual -y "${COUNT_MTRX_FILENAME}" -k "${OFSSET_MTRX_FILENAME}" -n 47 -j jj  -l "${NUM_TOTAL_SNP}" -m "${NUM_FEATURE_SNP}" -s "${PEAK_START}" -e "${PEAK_END}" -t -f "${PEAK}" -z);
	echo "${RASQUAL_RESULT}" >> $RESULT_FILENAME;
	
	printf "Peak start is %s and "${PEAK_START}";
	printf "Peak end os %s \n"${PEAK_END}";
done
