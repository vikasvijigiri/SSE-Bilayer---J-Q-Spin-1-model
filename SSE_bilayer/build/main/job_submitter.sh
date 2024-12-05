

#!/bin/bash

# Loop 10 times
for ((i=1; i<=10; i++)); do
  echo "Running sbatch job.sh - Iteration $i"
  sbatch jobscript_dummy.sh
done

