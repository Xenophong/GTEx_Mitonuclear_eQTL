
#!/usr/bin/env Rscript

#' -----------------------------------------------------------------------------
#' Script Title: QTL Analysis for mtDNA and Covariates
#' Description: This script performs QTL analysis on mitochondrial DNA (mtDNA)
#' and associated covariates across various tissue types. It generates commands 
#' for the Limix QTL analysis and performs post-processing on the results.
#' Author: Xenofon Giannoulis
#' -----------------------------------------------------------------------------

# List of tissues to be analyzed
tissues_tested=$(echo "Artery_Aorta" "Brain_Spinal_cord_cervical_c-1" "Minor_Salivary_Gland" "Vagina" "Brain_Putamen_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Hypothalamus" "Brain_Frontal_Cortex_BA9" "Brain_Caudate_basal_ganglia" "Artery_Coronary" "Spleen" "Brain_Nucleus_accumbens_basal_ganglia" "Prostate" "Brain_Cortex" "Brain_Cerebellum" "Pituitary" "Pancreas" "Stomach" "Colon_Sigmoid" "Testis" "Esophagus_Gastroesophageal_Junction" "Colon_Transverse" "Heart_Atrial_Appendage" "Heart_Left_Ventricle" "Breast_Mammary_Tissue" "Esophagus_Muscularis" "Adipose_Visceral_Omentum" "Cells_Cultured_fibroblasts" "Esophagus_Mucosa" "Skin_Not_Sun_Exposed_Suprapubic" "Lung" "Nerve_Tibial" "Artery_Tibial" "Adipose_Subcutaneous" "Thyroid" "Skin_Sun_Exposed_Lower_leg" "Whole_Blood" "Muscle_Skeletal")

# Working and output directories
wdir="/no-backup/xenofon/Emperor/22_GTEX/master"
odir='/no-backup/xenofon/Ascend/pipeline/mtDNA_cis_eQTL/secondary/3_logTPM_scale/'
kinship="$wdir/kin/kinships.all.mitonuc.grm.txt"
cov="allcovs"

# Loop through each tissue and generate QTL analysis commands
for tissue in $tissues_tested 
do 
    echo $tissue 

    outdir="$odir/$tissue"
    bfile="$wdir/mtgeno_allv8/per_tissue/$tissue.v8.eur.maf01"
    out="/no-backup/xenofon/Ascend/pipeline/mtDNA_cis_eQTL/secondary/3_logTPM_scale/src"
    phenofile="/no-backup/xenofon/Ascend/master/phenotypes/mtgenes/Load/3_logTPM_scale/${tissue}_mtgenes.txt"
    annofile="$wdir/anno/GTEX_MT_37genes_anno.txt"
    covfile="$wdir/covariates/limix_in/${tissue}_${cov}.txt"
    fvc="/no-backup/xenofon/Ascend/pipeline/mtDNA_cis_eQTL/secondary/3_logTPM_scale/fvc/${tissue}_fvc_input.txt"
    
    # Write command to run QTL analysis
    echo -e "/home/hpc/xenofon.giannoulis/anaconda3/envs/limix_qtl/bin/python \
    /home/hpc/xenofon.giannoulis/soft/limix/limix_qtl-master/Limix_QTL/run_QTL_analysis.py \
    --plink $bfile -af $annofile -pf $phenofile -od $outdir/ -cf $covfile -kf $kinship \
    -gm standardize -maf 0.05 -hwe 0 -np 1000 -cr 0.9 --cis --minimum_test_samples 60 \
    --no_chromosome_filter --window 20000 \
    --feature_variant_covariate $fvc" >> $out/secondary.limix.sh
done 

# POST PROCESS II
tissues_tested=$(echo "Artery_Aorta" "Brain_Spinal_cord_cervical_c-1" "Minor_Salivary_Gland" "Vagina" "Brain_Putamen_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Hypothalamus" "Brain_Frontal_Cortex_BA9" "Brain_Caudate_basal_ganglia" "Artery_Coronary" "Spleen" "Brain_Nucleus_accumbens_basal_ganglia" "Prostate" "Brain_Cortex" "Brain_Cerebellum" "Pituitary" "Pancreas" "Stomach" "Colon_Sigmoid" "Testis" "Esophagus_Gastroesophageal_Junction" "Colon_Transverse" "Heart_Atrial_Appendage" "Heart_Left_Ventricle" "Breast_Mammary_Tissue" "Esophagus_Muscularis" "Adipose_Visceral_Omentum" "Cells_Cultured_fibroblasts" "Esophagus_Mucosa" "Skin_Not_Sun_Exposed_Suprapubic" "Lung" "Nerve_Tibial" "Artery_Tibial" "Adipose_Subcutaneous" "Thyroid" "Skin_Sun_Exposed_Lower_leg" "Whole_Blood" "Muscle_Skeletal")
wdir='/no-backup/xenofon/Ascend/pipeline/mtDNA_cis_eQTL/secondary/3_logTPM_scale/'
python="/home/hpc/xenofon.giannoulis/anaconda3/envs/limix_qtl/bin/python"
pp="/home/hpc/xenofon.giannoulis/soft/limix/limix_qtl-master/Limix_QTL/post_processing/minimal_postprocess.py"

# Loop through each tissue and generate post-processing commands
for tissue in $tissues_tested 
do 
    echo $tissue 
    indir="$wdir/$tissue/"
    outdir1="$wdir/SFO/$tissue."
    outdir2="$wdir/TFB/$tissue."
    echo -e "$python $pp -id $indir -od $outdir1 -sfo" >> $wdir/SFO/sfo.sh
    echo -e "$python $pp -id $indir -od $outdir2 -tfb" >> $wdir/TFB/tfb.sh
done
