#!/bin/bash
#SBATCH --mail-user=yunkai.gao@petermac.org
#SBATCH --kill-on-invalid-dep=no
#SBATCH --time=00:30:00			#TO BE CHANGED
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL

#SBATCH --job-name="allc2"
#SBATCH --output="out/GSEA.out"
#SBATCH --error="out/GSEA.err"

echo $(date)
module load java/1.8.0_171-jdk java/1.8.0_171-jre

java -Xmx5000m -cp /home/ygao/Honours/GSEA/gsea3.0.jar xtools.gsea.Gsea \
-res /home/ygao/Honours/GSEA/Data/expressionData_for_GSEA.txt \
-cls /home/ygao/Honours/GSEA/Data/phenotype.label.cls#MOC_versus_BEN \
-gmx /home/ygao/Honours/GSEA/run.separate/all_c2_functional/c2.all.v6.2.symbols.gmt \
-collapse false \
-mode Max_probe \
-norm meandiv \
-nperm 1000 \
-permute phenotype \
-rnd_type no_balance \
-scoring_scheme weighted \
-rpt_label my_analysis \
-metric Signal2Noise \
-sort real \
-order descending \
-create_gcts false \
-create_svgs true \
-include_only_symbols true \
-make_sets true \
-median false \
-num 100 \
-plot_top_x 1000 \
-rnd_seed timestamp \
-save_rnd_lists false \
-set_max 500 \
-set_min 15 \
-zip_report false \
-out /home/ygao/Honours/GSEA/run.separate/all_c2_functional/out/$(date +%b%d) \
-gui false