import sys
import os
import json
import logging
from skimage import io
import numpy as np
from bigfish.detection import decompose_dense
from bigfish.stack import remove_background_gaussian


logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

configfile = sys.argv[1]
logging.info(f'reading config file {configfile}')

with open(configfile, 'r') as config_file:
    config = json.load(config_file)

scale = (200, 65, 65)
spot_radius = (1250, 170, 170)
sigma = (0.75, 2.3, 2.3)

# how many images do we have
nr_images = sum([len(exp['images']) for exp in config['experiments']])

n = 0
for exp in list(config["experiments"]):
    for img in list(exp["images"]):
        n = n + 1
        logging.info(f'processing image: {img["basename"]}.{img["format"]} [{n}/{nr_images}]')
        for ch in config["channels"]:
            mrna = ch['mrna']
            if mrna != "DAPI":

                logging.info(f'..mrna: {mrna}')

                dic_data = io.imread(img['dicfile'])
                mrna_data = io.imread(img[mrna]['rnafile'])
                dapi_data = io.imread(img['DAPI']['rnafile'])
                cell_mask_data = io.imread(img['cellmaskfile'])
                nuclear_mask_data = io.imread(img['nuclearmaskfile'])
                spot_data = np.load(img[mrna]['spotsfile'])

                mrna_filtered = remove_background_gaussian(mrna_data, sigma=sigma)

                spots, dense_regions, reference_spot = decompose_dense(
                    mrna_filtered,
                    spot_data,
                    voxel_size=scale,
                    spot_radius=spot_radius,
                    alpha=0.5, beta=2, gamma=1
                )

                img[mrna]['decompspotsfile'] = os.path.join(config['outputdir'], img['stem'], f'{mrna}_decomposed_spots.npy')
                img[mrna]['ddregionsfile'] = os.path.join(config['outputdir'], img['stem'], f'{mrna}_ddregions.npy')
                img[mrna]['refspotfile'] = os.path.join(config['outputdir'], img['stem'], f'{mrna}_rf_spots.tif')
                np.save(img[mrna]['decompspotsfile'], spots)
                np.save(img[mrna]['ddregionsfile'], dense_regions)
                io.imsave(img[mrna]['refspotfile'], reference_spot)
                logging.info(f"....decomposed spot data: {img[mrna]['decompspotsfile']}")
                logging.info(f"....dense regions data: {img[mrna]['ddregionsfile']}")
                logging.info(f"....reference spot: {img[mrna]['refspotfile']}")


logging.info(f'output config file: {configfile}')
with open(configfile, "w") as f:
    json.dump(config, f)

logging.info("done.")

