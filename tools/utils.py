import numpy as np
from copy import copy

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