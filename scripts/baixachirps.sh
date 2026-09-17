#!/bin/bash

BASE_URL="https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p05"

for year in $(seq 2003 2026); do
  FILE="chirps-v2.0.${year}.days_p05.nc"
  echo "Downloading $FILE..."
  wget -c "${BASE_URL}/${FILE}"
done
