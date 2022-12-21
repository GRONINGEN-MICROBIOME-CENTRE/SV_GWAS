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

d=$1
svtype=$2

meta_comb_dir=${d}/results_all_summary_stats/${svtype}/meta_combined_fdr/
nperm=10

mkdir ${meta_comb_dir}
echo -e "PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR" \
> ${TMPDIR}/eQTLs.txt

sort -m -k1,1g -T $TMPDIR -S 70G --buffer-size=1000 ${d}/results_all_summary_stats/${svtype}/meta/*/*.meta_res..eQTLs.txt \
>> ${TMPDIR}/eQTLs.txt

head -1000001 ${TMPDIR}/eQTLs.txt | gzip -c > ${meta_comb_dir}/eQTLs.txt.gz

# https://github.com/molgenis/systemsgenetics/wiki/QTL-mapping-pipeline
java -Xmx80g -Xms80g -jar ${d}/eqtl-mapping-pipeline-1.3.9-SNAPSHOT/eqtl-mapping-pipeline.jar \
--mode util \
--fdr \
--in  ${meta_comb_dir} \
--threshold 0.05 \
--perm 10 --nreqtls 1000000