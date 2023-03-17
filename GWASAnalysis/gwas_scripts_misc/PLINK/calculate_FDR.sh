#!/bin/bash
#SBATCH --job-name=fdr
#SBATCH --output=logs/fdr.out
#SBATCH --error=logs/fdr.err
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=85gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp=90gb

d=$1
svtype=$2

meta_comb_dir=${d}/results/${svtype}/meta_combined_fdr/
nperm=10

mkdir ${meta_comb_dir}
echo -e "PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR" \
> ${TMPDIR}/eQTLs.txt

# real results
while read line
do
    sp=$line
    echo "Merging and sorting SVs of species $sp"

    cmd="sort -m -k1,1g -T $TMPDIR -S 70G --buffer-size=1000 "
    for input in ${d}/results/dSV/meta/${sp}:*/*.meta_res.eQTLs.txt.gz; do
        cmd="$cmd <(gunzip -c '$input')"
    done
    echo $cmd
    eval "$cmd" | gzip -c > ${TMPDIR}/${sp}.eQTLs.txt.gz   
done < ${d}/data/${svtype}_species.txt

echo "Merging and sorting all per species files"
cmd="sort -m -k1,1g -T $TMPDIR -S 70G --buffer-size=1000 "
for input in ${TMPDIR}/*.eQTLs.txt.gz; do
    cmd="$cmd <(gunzip -c '$input')"
done
echo $cmd
eval "$cmd" | head -1000001  >> ${TMPDIR}/eQTLs.txt 
gzip -c ${TMPDIR}/eQTLs.txt  > ${meta_comb_dir}/eQTLs.txt.gz 

rm ${TMPDIR}/*


echo "permutations"
for p in `seq 1 $nperm`
do
    echo "Perm ${p}"
    echo -e "PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR" \
    > ${TMPDIR}/eQTLs.perm${p}.txt

    while read line
    do
        sp=$line
        echo "Merging and sorting permutation $p results of species $sp"
        cmd="sort -m -k1,1g -T $TMPDIR -S 70G --buffer-size=1000 "
        for input in ${d}/results/dSV/meta/${sp}:*/*.meta_res.eQTLs.perm${p}.txt.gz; do
            cmd="$cmd <(gunzip -c '$input')"
        done
        echo $cmd
        eval "$cmd" | gzip -c > ${TMPDIR}/${sp}.eQTLs.perm${p}.txt.gz    
  
    done < ${d}/data/${svtype}_species.txt
    
    echo "Merging and sorting all permutation $p for all species"
    cmd="sort -m -k1,1g -T $TMPDIR -S 70G --buffer-size=1000 "
    for input in ${TMPDIR}/*.eQTLs.perm${p}.txt.gz ; do
        cmd="$cmd <(gunzip -c '$input')"
    done
    echo $cmd
    eval "$cmd" | head -1000001  >> ${TMPDIR}/eQTLs.perm${p}.txt
    gzip -c ${TMPDIR}/eQTLs.perm${p}.txt  > ${meta_comb_dir}/PermutedEQTLsPermutationRound${p}.txt.gz
    rm ${TMPDIR}/*
done


java -Xmx80g -Xms80g -jar ~/tools/eqtl-mapping-pipeline-1.3.9-SNAPSHOT/eqtl-mapping-pipeline.jar \
--mode util \
--fdr \
--in  ${meta_comb_dir} \
--threshold 0.05 \
--perm 10 --nreqtls 1000000
