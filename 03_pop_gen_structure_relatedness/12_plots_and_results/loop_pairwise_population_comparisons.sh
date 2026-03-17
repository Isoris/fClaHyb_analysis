###
# loop_pairwise_population_comparisons.sh does this:
# Reads population names from a text file and prints all unique pairwise
# population combinations for downstream pairwise analyses.
###
#!/bin/bash
mapfile -t pops < ./02_info/pop.txt
num_pops=${#pops[@]}

for ((i=0; i<num_pops; i++)); do
  pop1="${pops[i]}"
  for ((j=i+1; j<num_pops; j++)); do
    pop2="${pops[j]}"
    echo "$pop1 $pop2"
    # Do your stuff
  done
done
