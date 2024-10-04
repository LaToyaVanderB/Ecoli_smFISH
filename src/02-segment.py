import sys
import logging
import json
import numpy as np
from pathlib import Path
import time
import re

# from cellpose_omni.core import use_gpu
from omnipose.gpu import use_gpu
from cellpose_omni import io, transforms
import cellpose_omni
from cellpose_omni import models
from cellpose_omni.models import MODEL_NAMES

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

configfile = sys.argv[1]
logging.info(f'reading config file {configfile}.')

with open(configfile, 'r') as f:
    config = json.load(f)

# This checks to see if you have set up your GPU properly.
# CPU performance is a lot slower, but not a problem if you
# are only processing a few images.
use_GPU = use_gpu()

# define parameters
chans = [0,0] #this means segment based on first channel, no second channel
params = {'channels': chans, # always define this with the model
          'rescale': None, # upscale or downscale your images, None = no rescaling
          'mask_threshold': -2, # erode or dilate masks with higher or lower values between -5 and 5
          'flow_threshold': 0, # default is .4, but only needed if there are spurious masks to clean up; slows down output
          'transparency': True, # transparency in flow output
          'omni': True, # we can turn off Omnipose mask reconstruction, not advised
          'cluster': True, # use DBSCAN clustering
          'resample': True, # whether or not to run dynamics on rescaled grid or original grid
          'verbose': False, # turn on if you want to see more output
          'tile': False, # average the outputs from flipped (augmented) images; slower, usually not needed
          'niter': None, # default None lets Omnipose calculate # of Euler iterations (usually <20) but you can tune it for over/under segmentation
          'augment': False, # Can optionally rotate the image and average network outputs, usually not needed
          'affinity_seg': False, # new feature, stay tuned...
         }

# DIC images
model = models.CellposeModel(gpu=use_GPU, model_type="cyto2")
params['min_size'] = 200

ticall = time.time()
# for f in files:
n = 0
for exp in config['experiments']:
    for img in exp['images']:
        n = n + 1
        f = img['dicfile']
        logging.info(f"segmenting {f} [{n}/{config['nr_images']}]")
        tic = time.time()
        mask, flow, style = model.eval(io.imread(f), **params)
        maskfile = f.replace('DIC.', f'DIC_masks.')
        maskfile_latest = f.replace('DIC.', f'DIC_masks_minsize={params["min_size"]}_maskth={params["mask_threshold"]}.')
        io.imwrite(maskfile_latest, mask)
        Path(maskfile).unlink(missing_ok=True)
        Path(maskfile).symlink_to(maskfile_latest)
        img['time']['02-segment-DIC'] = time.time() - tic
        logging.info(f"writing mask to {maskfile_latest}")
net_time = time.time() - ticall
logging.info(f"total DIC segmentation time: {net_time:.2f}s")

# DAPI images
model = models.CellposeModel(gpu=use_GPU, model_type="nuclei")
params['min_size'] = 10

ticall = time.time()
# for f in files:
n = 0
for exp in config['experiments']:
    for img in exp['images']:
        n = n + 1
        f = img['DAPI']['maxprojfile']
        logging.info(f"segmenting {f} [{n}/{config['nr_images']}]")
        tic = time.time()
        mask, flow, style = model.eval(io.imread(f), **params)
        maskfile = f.replace('DAPI_max_proj.', f'DAPI_masks.')
        maskfile_latest = f.replace('DAPI_max_proj.', f'DAPI_masks_minsize={params["min_size"]}_maskth={params["mask_threshold"]}.')
        io.imwrite(maskfile_latest, mask)
        Path(maskfile).unlink(missing_ok=True)
        Path(maskfile).symlink_to(maskfile_latest)
        img['time']['02-segment-DAPI'] = time.time() - tic
        logging.info(f"writing mask to {maskfile_latest}")
net_time = time.time() - ticall
logging.info(f"total DAPI segmentation time: {net_time:.2f}s")

logging.info(f"writing to config file: {configfile}")
with open(configfile, "w") as f:
    json.dump(config, f)