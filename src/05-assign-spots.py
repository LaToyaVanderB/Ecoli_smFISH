import sys
import os
import json
import time
import logging
from skimage import io
from skimage.segmentation import clear_border
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops
import argparse
from pathlib import Path

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)


def preprocess_spot_data(spot_data, dense_data):
    # spot_data has the form:
    # z, y, x

    # dense data has the form
    # z, y, x, mRNA counts, -- other information --

    # let's introduce mRNA counts of 1 for the spots:
    spot_data_padded = np.pad(spot_data, ((0,0),(0,1)), mode='constant', constant_values=1)

    # discard other information and merge
    spot_data_combined = np.concatenate([spot_data_padded, dense_data[:,:4]], axis=0)
    return spot_data_combined


def count_spots(mask, nuclear_mask, spot_data, cells):
    for z, y, x, number in spot_data:
        cell_id = mask[y, x]
        nucleus = nuclear_mask[y, x]

        if number == 1:
            cells[cell_id]['spots'] += number
        else:
            cells[cell_id]['dense_regions'] += 1
            cells[cell_id]['decomposed_RNAs'] += number

            # if the spot sits in the nucleus,
            # also increase nascent RNAs and transcription sites
            if nucleus > 0:
                cells[cell_id]['tx'] += 1
                cells[cell_id]['nascent_RNAs'] += number
    return cells


def count_nuclei(mask, nuclear_mask, cells):
    # count nuclei per cell - hyphae may have multiple ones!
    for nucleus in regionprops(nuclear_mask):
        y, x = nucleus.centroid
        cell_id = mask[int(y), int(x)]
        cells[cell_id]['nuclei'] += 1
    return cells


def spot_assignment(mask, expanded_mask, nuclear_mask, spot_data, dense_data):
    cells = {}

    for cell_id in np.unique(expanded_mask):
        cells[cell_id] = {
            'nuclei': 0,
            'spots': 0,
            'dense_regions': 0,
            'decomposed_RNAs': 0,
            'tx': 0,
            'nascent_RNAs': 0,
        }

    spot_data_combined = preprocess_spot_data(spot_data, dense_data)

    cells = count_spots(expanded_mask, nuclear_mask, spot_data_combined, cells)
    cells = count_nuclei(expanded_mask, nuclear_mask, cells)

    # remove spots on background
    del cells[0]

    # convert to dataframe, collect object information and merge
    df = pd.DataFrame(cells).T.reset_index().rename(columns={'index': 'label'})
    df['total_RNAs'] = df['spots'] + df['decomposed_RNAs'] - df['dense_regions']

    props = pd.DataFrame(regionprops_table(mask, properties=['label', 'bbox', 'area', 'eccentricity', 'axis_minor_length', 'axis_major_length', 'orientation', 'perimeter', 'solidity']))
    props_expanded = pd.DataFrame(regionprops_table(expanded_mask, properties=['label', 'bbox', 'area', 'eccentricity', 'axis_minor_length', 'axis_major_length', 'orientation', 'perimeter', 'solidity']))
    props = props.merge(props_expanded, on='label', how='right', suffixes=('', '_expanded'))
    df = props.merge(df, on='label')

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file (process multiple images).")
    parser.add_argument("-only", "--only",
                        help="VSI image file (process only one image. Requires a config file anyway).")
    args = parser.parse_args()
    logging.info(f'{args}')
    configfile = args.config
    only = args.only

    tic0 = time.time()
    logging.info(f'reading config file {configfile}')

    with open(configfile, 'r') as config_file:
        config = json.load(config_file)

    n = 0
    for exp in config['experiments']:
        for img in exp['images']:
            if ('detect_spots' not in img) or (img['detect_spots'] is True):
                if (only is None) or (img['stem'] == only):

                    n = n + 1
                    tic = time.time()

                    logging.info(f'processing image: {img['basename']}.{img['format']} [{n}/{config['nr_images']}]')

                    crop = False
                    if 'crop' in img:
                        crop = True
                        ymin, xmin, ymax, xmax = img['crop']
                        logging.info(f'cropping: top-left corner ({ymin}, {xmin}), bottom-right corner ({ymax}, {xmax}))')

                    img['resultfile'] = os.path.join(config['outputdir'], img['stem'], 'results.csv')
                    img['cellfile'] = os.path.join(config['outputdir'], img['stem'], 'cells.csv')
                    alldfs = pd.DataFrame()

                    for ch in config['channels']:
                        mrna = ch['mrna']
                        if mrna != "DAPI":
                            logging.info(f'..channel: {mrna}')

                            dic_data = io.imread(img['dicfile'])
                            mrna_data = io.imread(img[mrna]['rnafile'])
                            dapi_data = io.imread(img['DAPI']['rnafile'])
                            # cell_mask_data = io.imread(Path(img['cellmaskfile']).resolve())
                            cell_mask_data = io.imread(os.path.join(config['outputdir'], img['stem'], 'DIC_masks_selected_by_shape.tif'))
                            expanded_cell_mask_data = io.imread(os.path.join(config['outputdir'], img['stem'], 'DIC_masks_expanded.tif'))

                            nuclear_mask_data = io.imread(Path(img['nuclearmaskfile']).resolve())

                            # all the image files need to be cropped:
                            if crop is True:
                                dic_data = dic_data[ymin:ymax, xmin:xmax]
                                mrna_data = mrna_data[ymin:ymax, xmin:xmax]
                                dapi_data = dapi_data[ymin:ymax, xmin:xmax]
                                cell_mask_data = cell_mask_data[ymin:ymax, xmin:xmax]
                                expanded_cell_mask_data = expanded_cell_mask_data[ymin:ymax, xmin:xmax]
                                nuclear_mask_data = nuclear_mask_data[ymin:ymax, xmin:xmax]

                            # the spot files only contain cropped data:
                            spot_data = np.load(img[mrna]['spotsfile'])[:, 0:3]
                            dense_data = np.load(img[mrna]['ddregionsfile'])

                            # we also need to ignore the cell masks that touch the border
                            cell_mask_data = clear_border(cell_mask_data)
                            expanded_cell_mask_data = clear_border(expanded_cell_mask_data)


                            df = spot_assignment(cell_mask_data, expanded_cell_mask_data, nuclear_mask_data, spot_data, dense_data)
                            df.rename(columns={'label': 'image_cell_id'}, inplace=True)
                            cell_columns = ['image_cell_id',
                                            'bbox-0', 'bbox-1', 'bbox-2', 'bbox-3', 'area', 'eccentricity', 'axis_minor_length', 'axis_major_length', 'orientation', 'perimeter', 'solidity',
                                            'bbox-0_expanded', 'bbox-1_expanded', 'bbox-2_expanded', 'bbox-3_expanded', 'area_expanded', 'eccentricity_expanded', 'axis_minor_length_expanded', 'axis_major_length_expanded', 'orientation_expanded', 'perimeter_expanded', 'solidity_expanded',
                                            'nuclei']
                            df_cells = df[cell_columns]
                            rna_columns = ['image_cell_id', 'spots', 'dense_regions', 'decomposed_RNAs', 'tx', 'nascent_RNAs', 'total_RNAs']
                            df_rnas = df.loc[:, rna_columns]
                            df_rnas['mrna'] = mrna
                            df_rnas['strain'] = exp['strain']
                            df_rnas['condition'] = exp['condition']
                            df_rnas['seqnr'] = img['seqnr']
                            df_rnas.to_csv(f'{img["resultfile"]}.{mrna}', index=False)
                            df_rnas.drop(columns=['mrna', 'condition', 'strain', 'seqnr'], inplace=True)
                            df_rnas.rename(columns={
                                'spots': 'spots' + f'_{mrna}',
                                'dense_regions': 'dense_regions' + f'_{mrna}',
                                'decomposed_RNAs': 'decomposed_RNAs' + f'_{mrna}',
                                'tx': 'tx' + f'_{mrna}',
                                'nascent_RNAs': 'nascent_RNAs' + f'_{mrna}',
                                'total_RNAs': 'total_RNAs' + f'_{mrna}',
                            }, inplace=True)
                            if alldfs.empty:
                                df_cells.to_csv(img["cellfile"], index=False)
                                alldfs = df_cells
                            alldfs = alldfs.merge(df_rnas, on=['image_cell_id'], how='left')

                    alldfs['strain'], alldfs['condition'], alldfs['seqnr'] = exp['strain'].replace('Ecoli', 'MG1655'), exp['condition'], img['seqnr']

                    logging.info(f'saving data to {img['resultfile']}')
                    alldfs.to_csv(img['resultfile'], index=False)

                    img['time']['05-assign-spots'] = time.time() - tic

    logging.info(f'writing to config file: {configfile}')
    with open(configfile, "w") as f:
        json.dump(config, f)

    logging.info(f"processed {n} image{'s' if n > 1 else ''} in {time.time() - tic0:0.2f}s")
    logging.info(f"done.")


