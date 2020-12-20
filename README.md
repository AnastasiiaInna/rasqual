# **Allele-specific ChiP-seq analysis**

_This README is created by Anastasiia Hryhorzhevska_


The main problems that could appear while running RASQUAL on mac:

1. RASQUAL requires **CLAPACK** and **GSL** libraries to be installed in the environment. Even though **macOS** already ships with **BLAS** and **LAPACK** implementations in its **vecLib** framework, don't try to adjust RASQUAL code to be run 0using in-biuld libraries. It will take more time and effort. Just download it from http://www.netlib.org/clapack/, then follow the steps from RASQUAL github repo (Section "Installation tips for CLAPACK and GSL").

2. Always tabix the **vcf** files before running any RASQUAL function that requires it as an input parameter. To do it, run

```sh
tabix -pf vcf my_awesome_vcf.vcf.gz
```

3. VCF must be fased.

4. RASQUAL requires custom **vcf** files containing the allele specific counts of the target cellular trait at all SNPs. The files have to contain an additional subfield, "AS", located within the genotype field consisting of two integers, the reference and alternative allele counts, separated by a comma. To obtain the new **vcf**, one should run the function **createASVCF.sh** that is provided. However, for macOS, it should be slightly changed : 

Replace line nr 69 in the file **createASVCF.sh**

```sh
ASVCF_CURRENT_TMP_DIR=$(mktemp -t "${ASVCF_TMP_DIR:=${TMPDIR}}" -d "ASVCF_${VCF_FILENAME##*/}_XXXXXX");
```

with the line :

```sh
ASVCF_CURRENT_TMP_DIR=$(mktemp -d "ASVCF_${VCF_FILENAME##*/}_XXXXXX");
```

5. **Very important** : the order of the column names (sample names) must be in the same in each file, i.e. order of samples of count matrix is the same as that in vcf file.

## **Data preparation**

### **Get VCF files**

The **vcf** data are stored on cluster in : 

`/binder/anastasiia/chip-seq_signe/imputed`

1. Convert from **.haps** to **.vcf.gz** : 

Input data are **haps** files that were prepared by Jade and are stored in :

`binder/jade/Signe_genotypes_rephasing/04_shape`

To make **vcf** files we need to get the absolute paths of **haps** files and write them out to the text file. To do this, run :

```sh
pwd=/binder/jade/Signe_genotypes_rephasing/04_shape/

ls -d "$pwd"550k610k_samples_info06_maf_hwe_bggt_*.shapeit.phased.haps > haps_imputed_list.txt
```

Once we have the list of abolute paths, we can run the following function to create **vcf** files :

```sh
chmod +x convertHaps2vcf.sh

./convertHaps2vcf.sh  haps_imputed_list.txt '/home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed'
```

#### **`convertHaps2vcf.sh`**
>```sh
>#!/bin/bash
>
># The function to obtain VCFs from HAPs
>
># param 1 : text file with abolute paths to haps
># param 2 : location where to write out the vcf
>
>if [ ${#@} -ne 2 ]; then
>       printf "\nUsage: %s Not enough arguments\n" "${0}" > /dev/stderr;
>       exit 1;
>fi
>
>HAPS_FILENAME="${1}"
>OUTPUT_DIR="${2}"
>CHR_ID=0
>
>for HAP in $(cat "${HAPS_FILENAME}") ; do
>   CHR_ID=$((CHR_ID+1))
>   sbatch --wrap="shapeit -convert \
>                 --input-hap ${HAP} \
>                 --output-vcf '${OUTPUT_DIR}/chr${CHR_ID}.phased.vcf'"
>done
>```


2. Extract subset of individuals from vcfs

Input data : 
- list of samples to be extracted : one line - one sample 

`chipseq_dex_samples.txt` and `chipseq_veh_samples.txt`

- list of absolute paths to **vcf** files that were obtained in previous step : `vcf_list.txt`

To get get the absolute paths of **vcf** files and write them out to the text file, run :

```sh
pwd=/home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/

ls -d "$pwd"*.vcf> vcf_imputed_list.txt
```

Note :
**.txt** files should be in unix format. To be sure that they are, run the command : 

```sh 
dos2unix file.txt
```

Once we have the list of abolute paths to **vcf**, we can run the following function to extract the subset :

```sh
chmod +x extractSubsetVCF.sh
```

>**For dex**
```sh
./extractSubsetVCF.sh vcf_imputed_list.txt chipseq_dex_samples.txt '/home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/imputed_subset/dex'
```

> **For veh**
```sh
./extractSubsetVCF.sh vcf_imputed_list.txt chipseq_veh_samples.txt '/home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/imputed_subset/veh'
```

#### **`extractSubsetVCF.sh`**
>```sh
>#!/bin/bash
>
># Funtion to extract subset of samples from the given list of VCF files
>
># Param 1 : text file with abolute paths to vcf files
># Param 2 : text file with the list of samples to be extracted
># Param 3 : location where to write out the result
>
>if [ ${#@} -ne 3 ]; then
>       printf "\nUsage: %s Not enough arguments\n" "${0}" > /dev/stderr;
>       exit 1;
>fi
>
>VCF_FILENAME="${1}"
>SUBSET_LIST_FILENAME="${2}"
>OUTPUT_DIR="${3}"
>CHR_ID=0
>
>for FILE in $(cat "${VCF_FILENAME}") ; do
>    CHR_ID=$((CHR_ID+1));
>    bgzip ${FILE};
>    tabix -p vcf ${FILE}.gz;
>    bcftools view -Oz -S ${SUBSET_LIST_FILENAME} ${FILE}.gz > "${OUTPUT_DIR}"/"chr${CHR_ID}.phased.vcf.gz";
>    tabix -p vcf "${OUTPUT_DIR}"/"chr${CHR_ID}.phased.vcf.gz";
>done
>```

3. Merge multiple **vcf** files into single file

>This could be done before extracting the subset of samples.

The input files must be sorted by chr and position.

>**For dex**
```sh
cd /home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/imputed_subset/dex

bcftools concat -o phased_imputed_dex.vcf chr1.phased.vcf.gz chr2.phased.vcf.gz chr3.phased.vcf.gz chr4.phased.vcf.gz chr5.phased.vcf.gz chr6.phased.vcf.gz chr7.phased.vcf.gz chr8.phased.vcf.gz chr9.phased.vcf.gz chr10.phased.vcf.gz chr11.phased.vcf.gz chr12.phased.vcf.gz chr13.phased.vcf.gz chr14.phased.vcf.gz chr15.phased.vcf.gz chr16.phased.vcf.gz chr17.phased.vcf.gz chr18.phased.vcf.gz chr19.phased.vcf.gz chr20.phased.vcf.gz chr21.phased.vcf.gz chr22.phased.vcf.gz

bgzip phased_imputed_dex.vcf.gz
```

>**For veh**
```sh
cd /home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/imputed_subset/veh

bcftools concat -o phased_imputed_veh.vcf chr1.phased.vcf.gz chr2.phased.vcf.gz chr3.phased.vcf.gz chr4.phased.vcf.gz chr5.phased.vcf.gz chr6.phased.vcf.gz chr7.phased.vcf.gz chr8.phased.vcf.gz chr9.phased.vcf.gz chr10.phased.vcf.gz chr11.phased.vcf.gz chr12.phased.vcf.gz chr13.phased.vcf.gz chr14.phased.vcf.gz chr15.phased.vcf.gz chr16.phased.vcf.gz chr17.phased.vcf.gz chr18.phased.vcf.gz chr19.phased.vcf.gz chr20.phased.vcf.gz chr21.phased.vcf.gz chr22.phased.vcf.gz

bgzip phased_imputed_veh.vcf.gz
```
So far all actions were made on cluster. If you want to run RASQUAL locally on you PC, copy the obtained **vcf** files by running the following script :

>**For dex**
```sh
scp -Cp -r -o ProxyCommand="ssh cluser@biocomp.psych.mpg.de nc slurmgate 22"\
                            ahryhorzhevska@slurmgate:/home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/imputed_subset/dex/phased_imputed_dex.vcf.gz\
                            /Users/anastasiia_hry/github/rasqual/data/signe_chipseq/VCF_DATA/phased/imputed/dex/
```
>**For veh**
```sh
scp -Cp -r -o ProxyCommand="ssh cluser@biocomp.psych.mpg.de nc slurmgate 22"\
                            ahryhorzhevska@slurmgate:/home/ahryhorzhevska/mpip/datasets/signe_genotype_data/phased/imputed/imputed_subset/veh/phased_imputed_veh.vcf.gz\
                            /Users/anastasiia_hry/github/rasqual/data/signe_chipseq/VCF_DATA/phased/imputed/veh/
```
### **Calculate count matrices**

The count matrices for dex and veh cases are calculated by running R package *DiffBind*. The steps are descrideb in the attached Markdown code book. The result of the code is two text files that contains the count data for individual samples (no header). Each row represents the peakset and each column - sample. Please, rememebr that the order of the columns must be the same as in **vcf** file. The R code sorts the columns accordingly to the given order. You should have the list of ordered samples in a text file. 

The output files are : 

`count_mtrx_dex_for_rasqual.txt` and `count_mtrx_veh_for_rasqual.txt`.

The code also writes out the report of **DiffBind** analysis into **csv** format with statistics and annotations (if required).

## **Run RASQUAL**

Input data : 

- `phased_imputed_dex.vcf.gz`
- `phased_imputed_veh.vcf.gz`
- `count_mtrx_dex_for_rasqual.txt` 
- `count_mtrx_veh_for_rasqual.txt`
- `bam_list_dex.txt`
- `bam_list_veh.txt`
-`peakset_info_dex.txt`
-`peakset_info_veh.txt`

To run RASQUAL, please clone the git repo [RASQUAL](https://github.com/natsuhiko/rasqual), following the steps in teh begining of this README and place the count matrices in the folder `rasqual/data/signe_chipseq/` and **vcf** files along with **bam** list files into the folder `data/signe_chipseq/VCF_DATA/phased/imputed/dex/` for dex or `rasqual/data/signe_chipseq/VCF_DATA/phased/imputed/veh/` for veh.

To create a specific **vcf** file qith **AS** field, we need text file that contain absolute path to **bam** files. The order of the samples must be the same as that in the **vcf**.

1. Create new **vcf** file

```sh
RASQUALDIR=/Users/anastasiia_hry/github/rasqual/
cd $RASQUALDIR/src/ASVCF  
```

>**For dex**
```sh
tabix -pf vcf $RASQUALDIR/data/signe_chipseq/VCF_DATA/phased/imputed/dex/phased_dex_samples.vcf.gz

bash ./createASVCF.sh paired_end \
                      ../../data/signe_chipseq/VCF_DATA/bam_list_dex.txt\
                      ../../data/signe_chipseq/VCF_DATA/phased/unimputed_subset/phased_dex_samples.vcf.gz\
                      ../../data/signe_chipseq/VCF_DATA/phased/unimputed_subset/phased_dex_samples_as.vcf.gz \
                      atac
```

>**For veh**
```sh
tabix -p vcf $RASQUALDIR/data/signe_chipseq/VCF_DATA/phased/imputed/veh/phased_imputed_veh.vcf.gz 

bash ./createASVCF.sh paired_end \
                      ../../data/signe_chipseq/VCF_DATA/bam_list_veh.txt\
                      ../../data/signe_chipseq/VCF_DATA/phased/imputed/veh/phased_imputed_veh.vcf.gz\
                      ../../data/signe_chipseq/VCF_DATA/phased/imputed/veh/phased_imputed_veh_ac.vcf.gz\ 
                      atac
```
2. Create offset and binary files

>**For dex**

> Create offset file

```sh
cd $RASQUALDIR 

R --vanilla --quiet  --args data/signe_chipseq/count_mtrx_dex_for_rasqual.txt < R/makeOffset.R > log 

mv data/yourK.txt data/signe_chipseq/offset_dex.txt 
```
> Create binary files
```sh
R --vanilla --quiet --args data/signe_chipseq/count_mtrx_dex_for_rasqual.txt data/signe_chipseq/offset_dex.txt < R/txt2bin.R > log
```

>**For veh**

> Create offset file

```sh
cd $RASQUALDIR 

R --vanilla --quiet  --args data/signe_chipseq/count_mtrx_veh_for_rasqual.txt < R/makeOffset.R > log 

mv data/yourK.txt data/signe_chipseq/offset_veh.txt 
```
> Create binary files
```sh
R --vanilla --quiet --args data/signe_chipseq/count_mtrx_veh_for_rasqual.txt data/signe_chipseq/offset_veh.txt < R/txt2bin.R > log
```

3. Run RASQUAL

Once all data are prepared, run :

>**For dex**
```sh
DEXVCFDIR = ...
tabix -pf vcf $DEXVCFDIR/phased_imputed_dex_ac.vcf.gz

./runRasqualForSigne.sh \
      data/signe_chipseq/VCF_DATA/phased/imputed/dex/phased_imputed_dex_ac.vcf.gz\
      data/signe_chipseq/count_mtrx_dex_for_rasqual.bin\
      data/signe_chipseq/offset_dex.bin\
      data/signe_chipseq/peakset_info_dex.txt\
      rasqual_result_imputed_dex.txt
```
>**For veh**
```sh
VEHVCFDIR=...
tabix -pf vcf $VEHVCFDIR/phased_imputed_veh_ac.vcf.gz

./runRasqualForSigne.sh \
      data/signe_chipseq/VCF_DATA/phased/imputed/veh/phased_imputed_veh_ac.vcf.gz\
      data/signe_chipseq/count_mtrx_veh_for_rasqual.bin\
      data/signe_chipseq/offset_veh.bin\
      data/signe_chipseq/peakset_info_veh.txt\
      rasqual_result_imputed_veh.txt
```
#### **`runRasqualForSigne.sh`**
>```sh
>#!/bin/bash
>
># Funtion to run RASQUAL for the data
>
># Param 1 : newly create VCF file with AS field
># Param 2 : count matrix in binary format
># Param 3 : offset matrix in binary format
># Param 4 : text file that contains the name of all peaskets in the following format : CHR_NR:PEAK_START:PEAK_END
># Parma 5 : the name of text file where the result are recorded. 
>
>if [ ${#@} -ne 4 ]; then
>   printf "\nUsage: %s Not enough arguments\n" "${0}" > /dev/stderr;
>    exit 1;
>fi
>
>VCF_FILENAME="${1}"
>COUNT_MTRX_FILENAME="${2}"
>OFSSET_MTRX_FILENAME="${3}"
>PEAKSET_INFO_FILENAME="${4}"
>RESULT_FILENAME="${5}"
>
>if [ ! -f "${PEAKSET_INFO_FILENAME}" ] ; then
>        printf 'Error: the file "%s" does not exist.\n' "${PEAKSET_INFO_FILENAME}" > /dev/stderr;
>        exit 1;
>fi
>
>if [ -f "${RESULT_FILENAME}" ] ; then
>        rm "${RESULT_FILENAME}"
>fi
>
>touch "${RESULT_FILENAME}"
>
>echo Feature_ID rs_ID Chromosome SNP_position Ref Alt Allele_freq HWE_Chi-square Imputation_quality_score Log_10_Benjamini-Hochberg_Q-value Chi_square Effect_size_Pi Sequencing_error_rate Reference_all$
>
>jj=0
>for PEAK in $(cat "${PEAKSET_INFO_FILENAME}") ; do
>    jj=$((jj+1))
>    CHR=$(awk -F: '{print $1}' <<< "${PEAK}");
>    
>    PEAK_START=$(awk -F: '{print $2}' <<< "${PEAK}");    
>    PEAK_END=$(awk -F: '{print $3}' <<< "${PEAK}");
>    
>    GENE_START=$PEAK_START;
>    GENE_END=$PEAK_END;
>    
>    CIS_REG_START=$((PEAK_START-5000));
>    CIS_REG_END=$((PEAK_END+5000));
>    
>    L=$(tabix $VCF_FILENAME $CHR:$CIS_REG_START-$CIS_REG_END | wc -l | awk '{print $1}');
>    M=$(tabix $VCF_FILENAME $CHR:$GENE_START-$GENE_END | wc -l | awk '{print $1}');
>   
>    # RASQUAL_RESULT=$(tabix "${VCF_FIILENAME}" | bin/rasqual -y "${COUNT_MTRX_FILENAME}" -k "${OFSSET_MTRX_FILENAME}" -n 47 -j jj  -l "${NUM_TOTAL_SNP}" -m "${NUM_FEATURE_SNP}" -s "${PEAK_START}" -e "${PEAK_END}" -t -f "${PEAK}" -z);
>    RASQUAL_RESULT=$(tabix $VCF_FILENAME $CHR:$CIS_REG_START-$CIS_REG_END | bin/rasqual -y $COUNT_MTRX_FILENAME -k $OFSSET_MTRX_FILENAME -n 47 -j jj  -l $L -m $M -s $PEAK_START -e $PEAK_END -t -f $PEAK -z --no-posterior-update);    
>    echo "${RASQUAL_RESULT}" >> $RESULT_FILENAME;
>    
>    printf "Peak start is %s and "${PEAK_START}";
>    printf "Peak end is %s \n"${PEAK_END}";
> done
```
