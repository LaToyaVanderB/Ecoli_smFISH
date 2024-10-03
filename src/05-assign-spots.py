import sys
import os
import json
import logging
from skimage import io
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

configfile = sys.argv[1]
logging.info(f'reading config file {configfile}')

with open(configfile, 'r') as config_file:
    config = json.load(config_file)

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
            cells[cell_id]['spots_per_cell'] += number
        else:
            cells[cell_id]['dense_regions_per_cell'] += 1
            cells[cell_id]['decomposed_RNAs'] += number

            # if the spot sits in the nucleus,
            # also increase nascent RNAs and transcription sites
            if nucleus > 0:
                cells[cell_id]['tx_per_cell'] += 1
                cells[cell_id]['nascent_RNAs'] += number
    return cells


def count_nuclei(mask, nuclear_mask, cells):
    # count nuclei per cell - hyphae may have multiple ones!
    for nucleus in regionprops(nuclear_mask):
        y, x = nucleus.centroid
        cell_id = mask[int(y), int(x)]
        cells[cell_id]['nuclei'] += 1
    return cells


def spot_assignment(mask, nuclear_mask, spot_data, dense_data):
    cells = {}

    for cell_id in np.unique(mask):
        cells[cell_id] = {
            'spots_per_cell': 0,
            'dense_regions_per_cell': 0,
            'decomposed_RNAs': 0,
            'tx_per_cell': 0,
            'nascent_RNAs': 0,
            'nuclei': 0
        }

    spot_data_combined = preprocess_spot_data(spot_data, dense_data)

    cells = count_spots(mask, nuclear_mask, spot_data_combined, cells)
    cells = count_nuclei(mask, nuclear_mask, cells)

    # remove spots on background
    del cells[0]

    # convert to dataframe, collect object information and merge
    df = pd.DataFrame(cells).T.reset_index().rename(columns={'index': 'label'})
    df['total_RNAs_per_cell'] = df['spots_per_cell'] + df['decomposed_RNAs'] - df['dense_regions_per_cell']

    props = pd.DataFrame(regionprops_table(mask, properties=['label', 'bbox', 'area', 'eccentricity']))
    df = props.merge(df, on='label')

    return df


# how many images do we have
nr_images = sum([len(exp['images']) for exp in config['experiments']])

n = 0
for exp in list(config['experiments']):
    for img in list(exp['images']):
        n = n + 1
        logging.info(f'processing image: {img['basename']}.{img['format']} [{n}/{nr_images}]')
        img['resultfile'] = os.path.join(config['outputdir'], img['stem'], 'results.csv')
        alldfs = []

        for ch in config['channels']:
            mrna = ch['mrna']
            if mrna != "DAPI":
                logging.info(f'..channel: {mrna}')

                dic_data = io.imread(img['dicfile'])
                mrna_data = io.imread(img[mrna]['rnafile'])
                dapi_data = io.imread(img['DAPI']['rnafile'])
                cell_mask_data = io.imread(img['cellmaskfile'])
                nuclear_mask_data = io.imread(img['nuclearmaskfile'])
                spot_data = np.load(img[mrna]['decompspotsfile'])
                dense_data = np.load(img[mrna]['ddregionsfile'])

                df = spot_assignment(cell_mask_data, nuclear_mask_data, spot_data, dense_data)
                columns = df.columns
                df['strain'], df['condition'], df['seqnr'], df['mRNA'] = exp['strain'].replace('Ecoli', 'MG1655'), exp['condition'], img['seqnr'], mrna
                df = df[list(df.columns[-4:]) + list(df.columns[:-4])]
                alldfs.append(df)

        logging.info(f'saving data to {img['resultfile']}')
        pd.concat(alldfs).to_csv(img['resultfile'])

logging.info(f'output config file: {configfile}')
with open(configfile, "w") as f:
    json.dump(config, f)

logging.info("done.")


