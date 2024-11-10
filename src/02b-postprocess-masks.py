import json
import time
import logging
from skimage import io
import numpy as np
import argparse
from pathlib import Path


def remove_cell_masks_with_no_dapi(cell_masks, dapi_masks):
    cells_to_discard = []
    for label in np.unique(cell_masks):
        dapi_found = False
        cell_coordinates = np.where(cell_masks == label)
        for y, x in zip(cell_coordinates[0], cell_coordinates[1]):
            if dapi_masks[y, x] != 0:
                dapi_found = True
                break
        if dapi_found == False:
            cells_to_discard.append(label)
    print(f'removing {len(cells_to_discard)} cells with no dapi signal')
    cell_masks_dapi = copy(cell_masks)
    cell_masks_dapi[np.isin(cell_masks, cells_to_discard)] = 0
    return cell_masks_dapi, cells_to_discard


def filter_cells_by_area(masks, props, min_cell_area=200):
    discarded = props.query('area < @min_cell_area')['label']
    selected = set(props['label']) - set(discarded)
    masks_selected = copy(masks)
    masks_selected[np.isin(masks, list(discarded))] = 0
    masks_discarded = copy(masks)
    masks_discarded[np.isin(masks, list(selected))] = 0

    return masks_selected, masks_discarded, len(np.unique(masks_selected)), len(np.unique(masks_discarded)), len(np.unique(masks))


def filter_cells_by_shape(masks, props, min_clump_area=1000, max_clump_eccentricity=0.8):
    discarded = props.query('area > @min_clump_area').query('eccentricity < @max_clump_eccentricity')['label']
    selected = set(props['label']) - set(discarded)
    masks_selected = copy(masks)
    masks_selected[np.isin(masks, list(discarded))] = 0
    masks_discarded = copy(masks)
    masks_discarded[np.isin(masks, list(selected))] = 0

    return masks_selected, masks_discarded, len(np.unique(masks_selected)), len(np.unique(masks_discarded)), len(np.unique(masks))


def postprocess_mask(img, min_cell_area=200, min_clump_area=1000, min_clump_eccentricity=0.8):
    cell_masks = io.imread(img['cellmaskfile'])
    dapi_masks = io.imread(img['nuclearmaskfile'])
    # cell_masks_with_dapi = remove_cell_masks_with_no_dapi(cell_masks, dapi_masks)
    # io.imsave(Path(img['rootdir']) / 'DIC_masks_with_dapi.tif', cell_masks_with_dapi)
    imgdir = Path(config['outputdir']) / img['stem']
    cell_masks_with_dapi = io.imread(imgdir / 'DIC_masks_with_dapi.tif')

    masks_selected, masks_discarded, nr_selected, nr_discarded, nr_start = filter_cells_by_area(cell_masks_with_dapi, min_cell_area=200)
    io.imsave(imgdir / 'DIC_masks_discarded_by_area.tif', discarded)
    io.imsave(imgdir / 'DIC_masks_selected_by_area.tif', selected)
    logging.info(f'{imgdir}: discarded {nr_discarded} out of {nr_start} cells by area')

    masks_selected, masks_discarded, nr_selected, nr_discarded, nr_start = filter_cells_by_shape(masks_selected, min_clump_area=1000, max_clump_eccentricity=0.8)
    io.imsave(imgdir / 'DIC_masks_discarded_by_shape.tif', discarded)
    io.imsave(imgdir / 'DIC_masks_selected_by_shape.tif', selected)
    logging.info(f'{imgdir}: discarded {nr_discarded} out of {nr_start} cells by shape')


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
                    if (only is None) or (img['basename'] == only):
                        n = n + 1
                        tic = time.time()
                        logging.info(f'processing image: {img["basename"]}.{img["format"]} [{n}/{config['nr_images']}]')
                        postprocess_mask(img)

                        img['time']['02b-postprocess-masks'] = time.time() - tic

                        if only and (img['basename'] == only):
                            found = True
                            break


    logging.info(f'writing to config file: {configfile}')
    with open(configfile, "w") as f:
        json.dump(config, f)

    logging.info(f"processed {n} image{'s' if n > 1 else ''} in {time.time() - tic0:0.2f}s")
    logging.info(f"done.")

