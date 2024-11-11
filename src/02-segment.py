import sys
import logging
import json
import numpy as np
from pathlib import Path
import time

from omnipose.gpu import use_gpu
from cellpose_omni import io, transforms
from cellpose_omni import models

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

configfile = sys.argv[1]
logging.info(f'reading config file {configfile}.')

with open(configfile, 'r') as f:
    config = json.load(f)

# This checks to see if you have set up your GPU properly.
# CPU performance is a lot slower, but not a problem if you
# are only processing a few images.
use_GPU = use_gpu()


process_DIC = True
if process_DIC is True:

    # models for DIC images
    model_type1 = "cyto2"
    model1 = models.CellposeModel(gpu=use_GPU, model_type=model_type1)
    # possible channels options for model1: [2,0] (no nuclei)  or [2,3] (nuclei)

    model_type2 = "cyto2_omni"
    model2 = models.CellposeModel(gpu=use_GPU, model_type=model_type2)
    # possible channels options for model2: [1,2] (nuclei)

    model_menu = {
        "cyto2": model1,
        "cyto2_omni": model2
    }

    # defaults if nothing specified in config
    model = model1
    chans = [2,0]
    diameter = 0
    mask = 0
    flow = 0
    min_size = 0

    # pass parameters to model
    params = {'channels': chans, # always define this with the model
              'rescale': None, # upscale or downscale your images, None = no rescaling
              'mask_threshold': mask, # erode or dilate masks with higher or lower values between -5 and 5
              'flow_threshold': flow,
              'min_size': min_size,
              'diameter': diameter,
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

    ticall = time.time()
    # for f in files:
    n = 0
    for exp in config['experiments']:
        for img in exp['images']:
            if ('segment' not in img) or (img['segment'] is True):
                n = n + 1
                f = img['grgbfile']
                logging.info(f"segmenting {f} [{n}/{config['nr_images']}]")
                tic = time.time()

                if 'segmentation' in img:
                    params['diameter'] = img['segmentation']['diameter']
                    params['mask_threshold'] = img['segmentation']['mask']
                    params['flow_threshold'] = img['segmentation']['flow']
                    model = model_menu[img['segmentation']['model_type']]
                    params['channels'] = img['segmentation']['channels']
                    print(f"model_type from config-selected: {img['segmentation']['model_type']} {params['channels']}")

                mask, flow, style = model.eval(np.load(f), **params)
                maskfile_latest = Path(f).parent / f'DIC_masks_model={img["segmentation"]["model_type"]}_chan={str(params["channels"]).replace(" ", "")}_diameter={params["diameter"]}_minsize={params["min_size"]}_mask={params["mask_threshold"]}_flow={params["flow_threshold"]}.tif'
                io.imwrite(maskfile_latest, mask)
                Path(img['cellmaskfile']).unlink(missing_ok=True)
                Path(img['cellmaskfile']).symlink_to(maskfile_latest.parts[-1])
                img['time']['02-segment-DIC'] = time.time() - tic
                logging.info(f"writing mask to {maskfile_latest}")
    net_time = time.time() - ticall
    logging.info(f"total DIC segmentation time: {net_time:.2f}s")


# DAPI images
#
process_DAPI = True
if process_DAPI is True:
    model_type = "nuclei"
    model = models.CellposeModel(gpu=use_GPU, model_type=model_type)
    chans = [0, 0]
    min_size = 10
    # min_size = 0

    # pass parameters to model
    params = {'channels': chans, # always define this with the model
              'rescale': None, # upscale or downscale your images, None = no rescaling
              'mask_threshold': 0.0, # erode or dilate masks with higher or lower values between -5 and 5
              'flow_threshold': 0.0,
              'min_size': min_size,
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

    ticall = time.time()
    # for f in files:
    n = 0
    for exp in config['experiments']:
        for img in exp['images']:
            if ('segment' not in img) or (img['segment'] is True):
                n = n + 1
                f = img['DAPI']['maxprojfile']
                logging.info(f"segmenting {f} [{n}/{config['nr_images']}]")
                tic = time.time()
                mask, flow, style = model.eval(io.imread(f), **params)
                maskfile_latest = Path(f).parent / f'DAPI_masks_model={model_type}_chan={str(params["channels"]).replace(" ", "")}_diameter={params["diameter"]}_minsize={params["min_size"]}_mask={params["mask_threshold"]}_flow={params["flow_threshold"]}.tif'
                io.imwrite(maskfile_latest, mask)
                Path(img['nuclearmaskfile']).unlink(missing_ok=True)
                Path(img['nuclearmaskfile']).symlink_to(maskfile_latest.parts[-1])
                img['time']['02-segment-DAPI'] = time.time() - tic
                logging.info(f"writing mask to {maskfile_latest}")
    net_time = time.time() - ticall
    logging.info(f"total DAPI segmentation time: {net_time:.2f}s")

logging.info(f"writing to config file: {configfile}")
with open(configfile, "w") as f:
    json.dump(config, f)