# inputs: DIC and DAPI channels
# outputs: cell mask and nuclear mask

#conda activate omnipose

# cell masks
omnipose --dir tests/input/DIC --use_gpu --pretrained_model cyto2 --save_tif --in_folders
# this is an issue:
#mv tests/input/DIC/masks tests/output/

# nuclear masks
#omnipose --dir output/Zprojection --use_gpu --pretrained_model nuclei --save_tif --in_folders --omni --cluster
# cannot get nuclei model to work for now (pytorch not installed properly)
# using bact_phase_omni as fallback
omnipose --dir tests/input/DAPI/Zprojection --use_gpu --pretrained_model nuclei --save_tif
# best model so far:
