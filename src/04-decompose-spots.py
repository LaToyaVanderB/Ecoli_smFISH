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

    scale = tuple(config['scale'])
    spot_radius = tuple(config['spot_radius'])
    sigma = tuple(config['sigma'])

    n = 0
    found = False
    for exp in config["experiments"]:
        if found is False:
            for img in exp["images"]:
                if ('detect_spots' not in img) or (img['detect_spots'] is True):
                    if (only is None) or (img['stem'] == only):
                        n = n + 1
                        tic = time.time()
                        logging.info(f'processing image: {img["basename"]}.{img["format"]} [{n}/{config['nr_images']}]')

                        crop = False
                        if 'crop' in img:
                            crop = True
                            ymin, xmin, ymax, xmax = img['crop']
                            logging.info(f'cropping: top-left corner ({ymin}, {xmin}), bottom-right corner ({ymax}, {xmax}))')

                        for ch in config["channels"]:
                            mrna = ch['mrna']
                            if mrna != "DAPI":

                                logging.info(f'..mrna: {mrna}')

                                mrna_data = io.imread(img[mrna]['rnafile']) # CROP HERE
                                if crop is True:
                                    mrna_data = mrna_data[:, ymin:ymax, xmin:xmax]

                                # Commented out because not needed here:
                                # dic_data = io.imread(img['dicfile'])
                                # dapi_data = io.imread(img['DAPI']['rnafile'])
                                # cell_mask_data = io.imread(img['cellmaskfile'])
                                # nuclear_mask_data = io.imread(img['nuclearmaskfile'])

                                spot_data = np.load(img[mrna]['spotsfile'])[:, 0:3]

                                # was already computed in 03-detect-spots, reuse
                                # mrna_filtered = remove_background_gaussian(mrna_data, sigma=sigma)
                                mrna_filtered = np.load(img[mrna]['filteredmrnafile'])

                                # we should probably only do spot decomposition on spots that are in cells
                                spots, dense_regions, reference_spot = decompose_dense(
                                    mrna_filtered,
                                    spot_data,
                                    voxel_size=scale,
                                    spot_radius=spot_radius,
                                    alpha=0.5, # alpha impacts the number of spots per candidate region
                                    beta=2,    # beta impacts the number of candidate regions to decompose
                                    gamma=1    # gamma the filtering step to denoise the image
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

                        if only and (img['basename'] == only):
                            found = True
                            break


    logging.info(f'writing to config file: {configfile}')
    with open(configfile, "w") as f:
        json.dump(config, f)

    logging.info(f"processed {n} image{'s' if n > 1 else ''} in {time.time() - tic0:0.2f}s")
    logging.info(f"done.")

