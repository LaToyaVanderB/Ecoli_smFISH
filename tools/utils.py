import numpy as np
from copy import copy
from skimage.segmentation import expand_labels
from skimage.measure import regionprops_table
import pandas as pd
from typing import Tuple

from skimage.measure import regionprops


# with eternal thanks to the interwebs:
def translate_image(img, tx, ty):
    N, M = img.shape
    img_translated = np.zeros_like(img)

    img_translated[max(tx, 0):M + min(tx, 0), max(ty, 0):N + min(ty, 0)] = img[-min(tx, 0):M - max(tx, 0),
                                                                             -min(ty, 0):N - max(ty, 0)]

    img_translated_cropped = img[-min(tx, 0):M - max(tx, 0), -min(ty, 0):N - max(ty, 0)]
    return img_translated, img_translated_cropped


# very dumb implementation, needs fixing
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


def filter_cells_by_area(masks, props, min_cell_area=2000):
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


def expand_masks(masks, nr_of_pixels):
    return expand_labels(masks, distance=nr_of_pixels)


def get_regionprops(masks):
    properties = ['label', 'bbox', 'area', 'eccentricity', 'centroid', 'axis_major_length', 'axis_minor_length',
                  'orientation']
    return pd.DataFrame(regionprops_table(masks, properties=properties))


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
