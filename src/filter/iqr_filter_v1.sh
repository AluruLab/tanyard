declare -A iqrth=(
  [flower]=0.45
  [leaf]=0.5
  [root]=0.475
  [rosette]=0.425
  [seed]=0.525
  [seedling1wk]=0.425
  [seedling2wk]=0.4
  [shoot]=0.4
  [wholeplant]=0.425
)

for tissue in "${!iqrth[@]}"; do
  echo "Running for $tissue with threshold ${iqrth[$tissue]}"
  echo Rscript iqr_filter.R $1/${tissue}/${tissue}-final-mas5-qn.csv $1/${tissue}/${tissue}-final-mas5-qn-iqr-v1.csv ${iqrth[$tissue]}
  Rscript iqr_filter.R $1/${tissue}/${tissue}-final-mas5-qn.csv $1/${tissue}/${tissue}-final-mas5-qn-iqr-v1.csv ${iqrth[$tissue]}
done

