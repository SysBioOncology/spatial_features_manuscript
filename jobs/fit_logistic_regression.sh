#!/usr/bin/env bash
#SBATCH -J launch_fitting_logistic_regression
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

job_min=1
cpus_per_task=8

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/spatial_features_manuscript

output_dir="${work_dir}/output/models/grid"
model_name="spatial_features"
input_file="${work_dir}/data/CPTAC_prepped_spatial_features_data.rds"
param_grid=${work_dir}/data/metric_weight_grid.csv
# Hyper parameter tunign
# min_lambda=-6
min_alpha=0.1
n_alpha=5
n_lambda=100

# Repeated cross-validation
k_folds=5
n_repeats=10

# number of unique combinations to try
searchtype="grid"
# Determine job array limits
# A. Determine number of files
# job_max=$(ls -d -- $input_dir/* | wc -l) 2>/dev/null
# B. Number of lines in a file
job_max=$(wc -l < "${param_grid}")
# job_max=4
# job_max=2
echo $job_max

# Launch job array
sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J fitting_logistic_regression
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus_per_task}
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "spatial_features_manuscript"

# Read line in param grid file
line=\$(awk "NR==\${SLURM_ARRAY_TASK_ID}" ${param_grid})

# Select optimization metric and class weight accordingly
optimization_metric=\$(cut -d "," -f1 <<< \$line)
class_weight=\$(cut -d "," -f2 <<< \$line)

echo "Class weight: \$class_weight"
echo "Optimizatio metric: \$optimization_metric"

Rscript ${work_dir}/scripts/logistic_regression.R \
    --input_file $input_file \
    --output_dir $output_dir \
    --optimization_metric \$optimization_metric \
    --min_alpha $min_alpha \
    --n_alpha $n_alpha \
    --n_lambda $n_lambda \
    --k_folds $k_folds \
    --n_repeats $n_repeats \
    --model_name $model_name \
    --class_weight \$class_weight \
    --n_cores \${SLURM_CPUS_PER_TASK} \
    --searchtype $searchtype \

EOF

# done
