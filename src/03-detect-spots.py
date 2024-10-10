import json
import sys
import os.path
import time
from pathlib import Path
import logging
from skimage import io
import numpy as np
from bigfish.detection import detect_spots
from bigfish.stack import remove_background_gaussian
from bigfish.stack import get_in_focus_indices
from bigfish.stack import compute_focus
from scipy.signal import savgol_filter
from typing import Tuple

def find_high_density_patch(mask: np.ndarray, patch_size: Tuple = (200, 200), attempts: int = 20):
    """

    randomly samples patches on the mask image and returns the coordinates of the top-left corner
    of the densest patch

    :param mask: segmentation image expected to have dimension (h * w)
    :param patch_size: height and width of the patch
    :param attempts: how many patches to try

    :return: coordinates of top left corner of densest patch found
    :rtype: Tuple[int, int]

    """
    h, w = mask.shape
    h_patch, w_patch = patch_size

    cell_pixels = 0
    selected_patch = (None, None)  # top left corner
    for attempt in range(attempts):

        row_sample = np.random.randint(0, h - h_patch)
        col_sample = np.random.randint(0, w - w_patch)

        sample_patch = mask[row_sample:row_sample + h_patch, col_sample:col_sample + w_patch]
        if np.sum(sample_patch > 0) > cell_pixels:
            cell_pixels = np.sum(sample_patch > 0)
            selected_patch = (row_sample, col_sample)

    return selected_patch


def find_in_focus_indices(focus: np.ndarray, adjustment_bottom: int = 5, adjustment_top: int = 10):
    """

    find the in-focus indices of calculated focus scores

    :param focus: series of values representing max intensity along z-axis
    :param adjustment_bottom: controls by how much the resulting range should be padded (bottom)
    :param adjustment_top: controls by how much the resulting range should be padded (top)

    :return: low and high z-level between which the spots are in focus
    :rtype: Tuple[int,int]
    """

    # find the inflection points of the smoothed curve
    ifx_1 = min([np.diff(focus).argmax(), np.diff(focus).argmin()])
    ifx_2 = max([np.diff(focus).argmax(), np.diff(focus).argmin()])

    # add a little cushion to one side.
    ifx_1 -= adjustment_bottom
    ifx_2 += adjustment_top

    return ifx_1, ifx_2


logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

configfile = sys.argv[1]
logging.info(f'reading config file {configfile}')
with open(configfile, 'r') as f:
    config = json.load(f)

# parameters - adjust if necessary
scale = (200, 65, 65)
spot_radius = (1250, 170, 170)
sigma = (0.75, 2.3, 2.3)
patch_size = (200, 200)
detection_threshold = None  # set to None for automatic determination by bigFISH
# detection threshold is usually optimised manually and might need to be set differently per channel
# debug
# detection_threshold = 37

n = 0
for exp in config['experiments']:
    for img in exp['images']:
        n = n + 1
        tic = time.time()
        logging.info(f'processing image: {img['basename']}.{img['format']} [{n}/{config['nr_images']}]')

        # find high density region
        cell_mask_data = io.imread(img['cellmaskfile'])
        selected_patch = find_high_density_patch(cell_mask_data, patch_size=patch_size)
        logging.info(f'..selected patch: {selected_patch}')

        for ch in list(config['channels']):
            mrna = ch['mrna']
            if mrna != "DAPI":
                logging.info(f'..mrna: {mrna}')

                mrna_data = io.imread(img[mrna]['rnafile'])
                img_patch = mrna_data[:,
                            selected_patch[0]:selected_patch[0] + patch_size[0],
                            selected_patch[1]:selected_patch[1] + patch_size[1]
                            ]

                focus = compute_focus(img_patch)
                projected_focus = np.max(focus, axis=(1, 2))
                projected_focus_smoothed = savgol_filter(projected_focus, 16, 2, 0)
                ifx_1, ifx_2 = find_in_focus_indices(projected_focus_smoothed, adjustment_bottom=5, adjustment_top=10)
                ifx_1 = max(ifx_1, 0)
                ifx_2 = min(ifx_2, mrna_data.shape[0])
                img["z_max_focus"] = int(np.argmax(projected_focus_smoothed))
                img["ifx_1"] = int(ifx_1)
                img["ifx_2"] = int(ifx_2)

                logging.info(f'....in focus indices: [{ifx_1}, {ifx_2}] (max focus at slice {img["z_max_focus"]})')

                # need to save this if we compute it
                mrna_filtered = remove_background_gaussian(mrna_data, sigma=sigma)
                img[mrna]['filteredmrnafile'] = os.path.join(config['outputdir'], img['stem'], f'{mrna}_filtered.npy')
                np.save(img[mrna]['filteredmrnafile'], mrna_filtered)
                logging.info(f'....saving filtered mRNA image to file {img[mrna]['filteredmrnafile']}')


                spots, threshold = detect_spots(
                    mrna_filtered[ifx_1:ifx_2, ...],
                    threshold=detection_threshold,
                    voxel_size=scale,
                    spot_radius=spot_radius,
                    return_threshold=True
                )

                # restore z-level
                spots[:, 0] = spots[:, 0] + ifx_1

                # adjustable out-of focus filtering
                spots = spots[spots[:, 0] > ifx_1 + 2]

                img[mrna]['spotsfile_latest'] = os.path.join(config['outputdir'], img['stem'], f'{mrna}_spots_thr{threshold}.npy')
                np.save(img[mrna]['spotsfile_latest'], spots)
                img[mrna]['spotsfile'] = os.path.join(config['outputdir'], img['stem'], f'{mrna}_spots.npy')
                Path(img[mrna]['spotsfile']).unlink(missing_ok=True)
                Path(img[mrna]['spotsfile']).symlink_to(img[mrna]['spotsfile_latest'])
                logging.info(f'....saving {spots.shape[0]} spots to spots file {img[mrna]["spotsfile"]}')

        img['time']['03-detect-spots'] = time.time() - tic

        logging.info(f"....writing image parameters to {Path(config['outputdir']) / img['stem'] / 'img.json'}")
        with open(Path(config['outputdir']) / img['stem'] / 'img.json', 'w') as f:
            json.dump(img, f)

logging.info(f'writing to config file: {configfile}')
with open(configfile, "w") as f:
    json.dump(config, f)
logging.info("done.")
