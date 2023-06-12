#!/bin/bash
# Clara Gimeno Jésus, June 2022

# Bash script to download and run scripts to couple JULES model runs and MOSART routing model

cd jules-mosart/bin

# Adding some geotiff files into your folder 
cp /home/clara/JULES_output/jules-mosart-main/jules-mosart/bin/aux ~/$HOME/jules-mosart/jules-mosart-main/bin/aux

# Before running, you can change the rahu-config.yaml file and change the output directory

# Regrid jules output
./01_regrid_jules_output.sh rahu-config.yaml

# Aggregate jules output
./02_aggregate_jules_output.sh rahu-config.yaml

# Resample jules output
./03_resample_jules_output.sh rahu-config.yaml

# Compute jules annual mean runoffs
./04_compute_jules_annual_mean_runoff.sh rahu-config.yaml

# Combine jules runoff
./05_combine_jules_runoff.sh rahu-config.yaml

# Upscale flow direction
./06_upscale_flow_direction.sh rahu-config.yaml

# Create mosart input
./07_create_mosart_input.sh rahu-config.yaml

## Running mosart
# Copy 
# cp /mnt/scratch/scratch/data/ForClara/config_rahu.yaml /$HOME/jules-mosart

# Modify yaml files with appropriate paths

# Run python script, make sure to make changes to config_rahu.yaml file first
./mosart-rahu.py

# Pull the simulated discharge time series at one of the gauging stations from the mosart output
./extract-qsim.py

# Make plots
./process-qobs.R
