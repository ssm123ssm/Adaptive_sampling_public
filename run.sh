set -euo pipefail

Rscript simple_framework_simulation.R \
  --replicates=100 \
  --particles=1000 \
  --design-mc=1000 \
  --validation-mc=5000 \
  --selection-confidence=0.95 \
  --workers=4 \
  --seed=20260729

Rscript plot_framework_behaviour.R

OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
Rscript sir_variant_simulation.R \
  --replicates=100 \
  --particles=1200 \
  --design-mc=600 \
  --validation-mc=3000 \
  --selection-confidence=0.95 \
  --workers=6 \
  --seed=20260801
