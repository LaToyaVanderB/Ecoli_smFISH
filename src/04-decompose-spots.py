import os
import json
import time
import logging
from skimage import io
import numpy as np
from bigfish.detection import decompose_dense
from bigfish.stack import remove_background_gaussian
import argparse

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file (process multiple images).")
    parser.add_argument("-only", "--only", help="VSI image file (process only one image. Requires a config file anyway).")
    args = parser.parse_args()
    logging.info(f'{args}')
    configfile = args.config
    only = args.only

    tic0 = time.time()
    logging.info(f'reading config file {configfile}')

    with open(configfile, 'r') as config_file:
        config = json.load(config_file)

    scale = (200, 65, 65)
    spot_radius = (1250, 170, 170)
    sigma = (0.75, 2.3, 2.3)

    n = 0
    for exp in config["experiments"]:
        for img in exp["images"]:
            if (only is None) or (img['sourcefile'] == only):
                n = n + 1
                tic = time.time()
                logging.info(f'processing image: {img["basename"]}.{img["format"]} [{n}/{config['nr_images']}]')
                for ch in config["channels"]:
                    mrna = ch['mrna']
                    if mrna != "DAPI":

                        logging.info(f'..mrna: {mrna}')

                        dic_data = io.imread(img['dicfile'])
                        mrna_data = io.imread(img[mrna]['rnafile'])
                        # dapi_data = io.imread(img['DAPI']['rnafile'])
                        # cell_mask_data = io.imread(img['cellmaskfile'])
                        # nuclear_mask_data = io.imread(img['nuclearmaskfile'])
                        spot_data = np.load(img[mrna]['spotsfile'])

                        # was already computed in 03-detect-spots, reuse
                        mrna_filtered = remove_background_gaussian(mrna_data, sigma=sigma)
                        # mrna_filtered = np.load(img[mrna]['filteredmrnafile'])

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

                img['time']['04-decompose-spots'] = time.time() - tic

    logging.info(f'writing to config file: {configfile}')
    with open(configfile, "w") as f:
        json.dump(config, f)

    logging.info(f"processed {n} image{'s' if n > 1 else ''} in {time.time() - tic0:0.2f}s")
    logging.info(f"done.")

