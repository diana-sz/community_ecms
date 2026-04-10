#!/bin/bash
# IMPORTANT: Run this script from the ecmtool directory where main.py is located

models=(
  "Rhodoferax_mEmC_13_ex"
  "Rhodoferax_mEmC_24_ex"
  "Rhodoferax_mEmC_5_ex"
  "Rhodoferax_mEmC_6_ex"
  "Rhodoferax_mE_1_ex"
  "Geobacter_mEmC_12_ex"
  "Geobacter_mEmC_34_ex"
  "Geobacter_mEmC_56_ex"
  "Geobacter_mE_1_ex"
)

for model in "${models[@]}"; do
  python main.py \
    --model_path "../core_models/uranium/models/${model}.xml" \
    --compress true \
    --remove_infeasible false \
    --out_path "../core_models/uranium/results/${model%_ex}_full_conversions_ex_compress.csv"
done