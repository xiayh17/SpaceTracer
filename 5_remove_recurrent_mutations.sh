#!/bin/bash

export src_path=$1
export pred_sites=$2
export feature_path=$3
export outpath=$4
export germ_path=$5

export genome_fa=$6
export RNA_editing_path=$7
export gtexGene=$8
export dbSNP=$9
export imprinted_genes=${10}
export PON_file=${11}

prefix="pred"
mkdir -p ${outpath}
python ${src_path}/others/check_polymers.py \
        --mutations  ${pred_sites} \
        --fasta ${genome_fa}\
        --out_file ${outpath}/${prefix}_filter_poly.identifier.txt

bedtools intersect -a <(cat ${outpath}/${prefix}_filter_poly.identifier.txt|sed 's/_/\t/g'|awk '{print $1,$2,$2,$3,$4}' OFS="\t") \
                   -b ${imprinted_genes} -v |sort -u|awk '{print $1, $2,$4,$5}' OFS="_"  > ${outpath}/${prefix}.filter_poly.filter_imprinted.identifier.txt

perl -e 'open(IN1, shift); open(IN2, shift); while(<IN1>){chomp; @A=split/\s+/; $hash{$A[0]}=1} while(<IN2>){chomp; @B=split/\s+/; print $_,"\n" if exists $hash{$B[0]}}' <(sed '1iidentifier' ${outpath}/${prefix}.filter_poly.filter_imprinted.identifier.txt) ${feature_path}>  ${outpath}/${prefix}.filter_poly.filter_imprinted.feature.txt


# getASE
python ${src_path}/others/get_ASE.py \
        --germline ${germ_path}  \
        --outdir ${outpath}/ \
        --gtexGene ${gtexGene} \
        --dbSNP ${dbSNP}

ase_file="${outpath}/candidate_ASE_sites.txt"

python ${src_path}/others/filter_ASE.py \
    --ase_file ${ase_file} \
    --features ${outpath}/${prefix}.filter_poly.filter_imprinted.feature.txt \
    --outdir ${outpath} \
    --RNA_editing /storage/douyanmeiLab/yangzhirui/Reference/COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed \
    --prefix "${prefix}.filter_poly.filter_imprinted"

awk '{print $1}' ${PON_file} | sort -u >  ${outpath}/temp_exclude.txt
grep -vFf  ${outpath}/temp_exclude.txt ${outpath}/${prefix}.filter_poly.filter_imprinted.filter_ase.filter_editing.identifier.txt > ${outpath}/${prefix}.filter_poly.filter_imprinted.filter_ase.filter_editing.filter_PON.identifier.txt
cp ${outpath}/${prefix}.filter_poly.filter_imprinted.filter_ase.filter_editing.filter_PON.identifier.txt ${outpath}/${prefix}.FINAL.txt