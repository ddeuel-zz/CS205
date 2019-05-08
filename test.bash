grain_sizes=(100 150 200 250 300)
num_steps_exponents=(3 4 5)

OUTPUT_FILE="acc_timings.csv"

# Write header
echo "Elapsed time (s), Grain size, Number of steps" > ${OUTPUT_FILE}
  
for nm_stps_expnnt in ${num_steps_exponents[*]}; do
	for grn_sz in ${grain_sizes[*]}; do
    num_steps=$((10**nm_stps_expnnt))
  	output_line=$(./grain_acc ${grn_sz} ${num_steps} 1)
    echo "${output_line}" >> ${OUTPUT_FILE}
  done
done