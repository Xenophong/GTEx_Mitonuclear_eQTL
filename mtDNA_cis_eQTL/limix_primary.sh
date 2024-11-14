#!/usr/bin/env Rscript

#' -----------------------------------------------------------------------------
#' Script Title: Setup and Run mtDNA cis eQTL Analysis
#' Description: This script sets up the necessary directory structure for the 
#'              analysis of mitochondrial DNA (mtDNA) cis-expression quantitative 
#'              trait loci (eQTL). It prepares directories for various tissues, 
#'              runs the Limix software for QTL analysis, and processes the output 
#'              for downstream analysis.
#' Author: Xenofon Giannoulis
#' Email: xenofon.giannoulis@helmholtz-muenchen.de
#' -----------------------------------------------------------------------------

# Create output directory for cis eQTL
cases=$(echo "3_logTPM_scale/")
tissues=$(echo "Brain_Substantia_nigra" "Uterus" "Cells_EBV-transformed_lymphocytes" "Brain_Spinal_cord_cervical_c-1" "Minor_Salivary_Gland" "Brain_Amygdala" "Vagina" "Brain_Anterior_cingulate_cortex_BA24" "Ovary" "Small_Intestine_Terminal_Ileum" "Brain_Hippocampus" "Brain_Putamen_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Hypothalamus" "Brain_Frontal_Cortex_BA9" "Brain_Caudate_basal_ganglia" "Artery_Coronary" "Liver" "Spleen" "Brain_Nucleus_accumbens_basal_ganglia" "Prostate" "Brain_Cortex" "Brain_Cerebellum" "Adrenal_Gland" "Pituitary" "Pancreas" "Stomach" "Colon_Sigmoid" "Testis" "Esophagus_Gastroesophageal_Junction" "Colon_Transverse" "Heart_Atrial_Appendage" "Artery_Aorta" "Heart_Left_Ventricle" "Breast_Mammary_Tissue" "Esophagus_Muscularis" "Adipose_Visceral_Omentum" "Cells_Cultured_fibroblasts" "Esophagus_Mucosa" "Skin_Not_Sun_Exposed_Suprapubic" "Lung" "Nerve_Tibial" "Artery_Tibial" "Adipose_Subcutaneous" "Thyroid" "Skin_Sun_Exposed_Lower_leg" "Whole_Blood" "Muscle_Skeletal")
dir="/path/to/your/output/directory"
for case in $cases
do 
    echo $case
    mkdir $dir/$case
    for tissue in $tissues 
    do 
        echo $tissue 
        mkdir $dir/$case/$tissue
        mkdir $dir/$case/SFO
        mkdir $dir/$case/TFB
        mkdir $dir/$case/src
    done
done

#############################    ##########################################################################################################
############################# primary ##########################################################################################################
#############################    ##########################################################################################################

# 1.1 Setup and Run Limix
tissues=$(echo "Brain_Substantia_nigra" "Uterus" "Cells_EBV-transformed_lymphocytes" "Brain_Spinal_cord_cervical_c-1" "Minor_Salivary_Gland" "Brain_Amygdala" "Vagina" "Brain_Anterior_cingulate_cortex_BA24" "Ovary" "Small_Intestine_Terminal_Ileum" "Brain_Hippocampus" "Brain_Putamen_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Hypothalamus" "Brain_Frontal_Cortex_BA9" "Brain_Caudate_basal_ganglia" "Artery_Coronary" "Liver" "Spleen" "Brain_Nucleus_accumbens_basal_ganglia" "Prostate" "Brain_Cortex" "Brain_Cerebellum" "Adrenal_Gland" "Pituitary" "Pancreas" "Stomach" "Colon_Sigmoid" "Testis" "Esophagus_Gastroesophageal_Junction" "Colon_Transverse" "Heart_Atrial_Appendage" "Artery_Aorta" "Heart_Left_Ventricle" "Breast_Mammary_Tissue" "Esophagus_Muscularis" "Adipose_Visceral_Omentum" "Cells_Cultured_fibroblasts" "Esophagus_Mucosa" "Skin_Not_Sun_Exposed_Suprapubic" "Lung" "Nerve_Tibial" "Artery_Tibial" "Adipose_Subcutaneous" "Thyroid" "Skin_Sun_Exposed_Lower_leg" "Whole_Blood" "Muscle_Skeletal")
wdir="/path/to/your/working/directory"
dir="/path/to/your/output/directory"
kinship="$wdir/kin/kinships.all.mitonuc.grm.txt"
cov=$(echo "allcovs")

for case in $cases
do 
    echo $case
    for tissue in $tissues 
    do 
        echo $tissue 

        bfile="$wdir/mtgeno_allv8/per_tissue/$tissue.v8.eur.maf01"
        outdir="$dir/$case/$tissue"
        annofile="$wdir/anno/GTEX_MT_37genes_anno.txt"
        covfile="$wdir/covariates/limix_in/${tissue}_${cov}.txt"

        phenofile="/path/to/your/phenotype/file/${case}/${tissue}_mtgenes.txt"
        out_ind="$dir/$case/src"

        echo -e "/path/to/python/executable/python \
        /path/to/limix/run_QTL_analysis.py \
        --plink $bfile -af $annofile -pf $phenofile -od $outdir/ -cf $covfile -kf $kinship \
        -gm standardize -maf 0.01 -hwe 0 -np 1000 -cr 0.9 --cis --minimum_test_samples 60 \
        --no_chromosome_filter --window 20000" >> $out_ind/$case.primary.limix.sh

    done 
done

# 1.2 Setup and run POST PROCESS I
cases=$(echo "3_logTPM_scale/")
tissues=$(echo "Brain_Substantia_nigra" "Uterus" "Cells_EBV-transformed_lymphocytes" "Brain_Spinal_cord_cervical_c-1" "Minor_Salivary_Gland" "Brain_Amygdala" "Vagina" "Brain_Anterior_cingulate_cortex_BA24" "Ovary" "Small_Intestine_Terminal_Ileum" "Brain_Hippocampus" "Brain_Putamen_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Hypothalamus" "Brain_Frontal_Cortex_BA9" "Brain_Caudate_basal_ganglia" "Artery_Coronary" "Liver" "Spleen" "Brain_Nucleus_accumbens_basal_ganglia" "Prostate" "Brain_Cortex" "Brain_Cerebellum" "Adrenal_Gland" "Pituitary" "Pancreas" "Stomach" "Colon_Sigmoid" "Testis" "Esophagus_Gastroesophageal_Junction" "Colon_Transverse" "Heart_Atrial_Appendage" "Artery_Aorta" "Heart_Left_Ventricle" "Breast_Mammary_Tissue" "Esophagus_Muscularis" "Adipose_Visceral_Omentum" "Cells_Cultured_fibroblasts" "Esophagus_Mucosa" "Skin_Not_Sun_Exposed_Suprapubic" "Lung" "Nerve_Tibial" "Artery_Tibial" "Adipose_Subcutaneous" "Thyroid" "Skin_Sun_Exposed_Lower_leg" "Whole_Blood" "Muscle_Skeletal")
wdir='/path/to/your/output/directory/3_logTPM_scale'
python="/path/to/python/executable/python"
pp="/path/to/limix/minimal_postprocess.py"

for tissue in $tissues 
do 
    echo $tissue 
    indir="$wdir/$tissue/"
    outdir1="$wdir/SFO/$tissue."
    outdir2="$wdir/TFB/$tissue."
    echo -e "$python $pp -id $indir -od $outdir1 -sfo" >> $wdir/SFO/sfo.sh
    echo -e "$python $pp -id $indir -od $outdir2 -tfb" >> $wdir/TFB/tfb.sh
done
