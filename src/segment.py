import sys
import logging
import json
from pathlib import Path
import time

from omnipose.gpu import use_gpu
from cellpose_omni import io, transforms
from cellpose_omni import models
# from cellpose_omni.models import MODEL_NAMES

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

inputdir = sys.argv[1]
logging.info(f'input dir: {inputdir}')

# This checks to see if you have set up your GPU properly.
# CPU performance is a lot slower, but not a problem if you
# are only processing a few images.
use_GPU = use_gpu()

# define parameters
chans = [0,0] #this means segment based on first channel, no second channel
params = {'channels': chans, # always define this with the model
          'rescale': None, # upscale or downscale your images, None = no rescaling
          'mask_threshold': 0.0, # erode or dilate masks with higher or lower values between -5 and 5
          'flow_threshold': 0.0,
          'diameter': 0.0,
          'invert': False,
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
params['min_size'] = 0

ticall = time.time()
# for f in files:
n = 0
files = list(Path(inputdir).glob('*.tif'))
for f in files:
        n = n + 1
        logging.info(f"segmenting {f} [{n}/{len(files)}]")
        tic = time.time()
        mask, flow, style = model.eval(io.imread(f), **params)
        maskfile = f.parent / (f.stem + '_masks' + f.suffix)
        maskfile_latest = f.parent / ( f.stem + f'_masks_minsize={params["min_size"]}_maskth={params["mask_threshold"]}' + f.suffix )
        io.imwrite(maskfile_latest, mask)
        maskfile.unlink(missing_ok=True)
        maskfile.symlink_to(maskfile_latest)
        logging.info(f"time: {time.time() - tic}")
        logging.info(f"writing mask to {maskfile_latest}")
net_time = time.time() - ticall
logging.info(f"total DIC segmentation time: {net_time:.2f}s")

