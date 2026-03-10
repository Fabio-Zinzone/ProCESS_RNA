#!/bin/bash

SIF_PATH="/home/hp/app_sif/process_rr_abc.sif"
PROJECT_DIR=$(pwd)

export APPTAINERENV_R_MAKEVARS_USER="/mnt/Makevars"

apptainer exec --cleanenv \
	--bind "$PROJECT_DIR":"/mnt" \
	--pwd "/mnt" \
	"$SIF_PATH" \
	Rscript --vanilla /mnt/pro_rr.R
