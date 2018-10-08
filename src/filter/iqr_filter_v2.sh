declare -A iqrth=(
  [flower]=0.475
  [leaf]=0.525
  [root]=0.5
  [rosette]=0.45
  [seed]=0.55
  [seedling1wk]=0.45
  [seedling2wk]=0.425
  [shoot]=0.425
  [wholeplant]=0.45
)

for tissue in "${!iqrth[@]}"; do
  echo "Running for $tissue with threshold ${iqrth[$tissue]}"
  echo Rscript iqr_filter.R $1/${tissue}/${tissue}-final-mas5-qn.csv $1/${tissue}/${tissue}-final-mas5-qn-iqr-v2.csv ${iqrth[$tissue]}
  Rscript iqr_filter.R $1/${tissue}/${tissue}-final-mas5-qn.csv $1/${tissue}/${tissue}-final-mas5-qn-iqr-v2.csv ${iqrth[$tissue]}
done

