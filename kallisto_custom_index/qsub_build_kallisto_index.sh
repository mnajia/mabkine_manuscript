#! /bin/bash
#Build a custom kallisto index
#Mohamad Najia


output_dir=/broad/blainey_lab/Mo/exps/20251010_MN_ACS_NZ_SS_novaseq/kallisto_custom_index
index_name=kallisto_HER2_CAR_custom_index.idx
cdna=$output_dir/hg19.annot.cdna.HER2.CAR.custom
k=31

cd $output_dir
cat HER2.CAR.GFP.cdna.fa /broad/blainey_lab/Mo/genomes/hg19.annot.cdna > $cdna


qsub -N kallisto_index -o $output_dir/out.log -e $output_dir/err.log -m ea -M mnajia@broadinstitute.org /broad/blainey_lab/Mo/scripts/build_kallisto_index.sh $index_name $cdna $k
