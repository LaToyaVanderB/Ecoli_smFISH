import sys
import logging
from glob import glob
from pathlib import Path
from skimage import io
import numpy as np
import json
import napari

# need to do this in napari console otherwise the import does not work
# because life is hard and then you die
# import sys; sys.path.insert(0, "")

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

modes = {
    'all': 'layers_ordered',
    'light': 'layers_light',
    'seg': 'layers_segmentation',
    'cells': 'layers_cells',
}

layers_ordered = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'name': 'DIC_masks', 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'name': 'DIC', 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'name': 'DAPI_masks', 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'name': 'DAPI', 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'name': 'DAPI_max_proj', 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rpoD', 'type': 'image', 'properties': { 'name': 'rpoD', 'colormap': 'magenta', 'visible': True, 'blending': 'additive' } },
    { 'name': 'rpoD_filtered', 'type': 'rna_filtered', 'properties': { 'name': 'rpoD_filtered', 'blending': 'additive', 'visible': False } },
    { 'name': 'rpoD_spots', 'type': 'spots', 'properties': { 'name': 'rpoD_spots', 'blending': 'additive', 'visible': False } },
    { 'name': 'rnlAB', 'type': 'image', 'properties': { 'name': 'rnlAB', 'colormap': 'cyan', 'visible': True, 'blending': 'additive' } },
    { 'name': 'rnlAB_filtered', 'type': 'rna_filtered', 'properties': { 'name': 'rnlAB_filtered', 'blending': 'additive', 'visible': False } },
    { 'name': 'rnlAB_spots', 'type': 'spots', 'properties': { 'name': 'rnlAB_spots', 'blending': 'additive', 'visible': False } },
    { 'name': 'hipBA', 'type': 'image', 'properties': { 'name': 'hipBA', 'colormap': 'yellow', 'visible': True, 'blending': 'additive' } },
    { 'name': 'hipBA_filtered', 'type': 'rna_filtered', 'properties': { 'name': 'hipBA_filtered', 'blending': 'additive', 'visible': False} },
    { 'name': 'hipBA_spots', 'type': 'spots', 'properties': { 'name': 'hipBA_spots', 'blending': 'additive', 'visible': False} }
]

layers_light = [
    { 'name': 'DIC', 'type': 'image', 'properties': { 'name': 'DIC', 'colormap': 'grey', 'visible': True, 'blending': 'translucent_no_depth', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'name': 'DAPI', 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
]

layers_segmentation = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'name': 'DIC_masks', 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'name': 'DIC', 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'name': 'DAPI_masks', 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'name': 'DAPI', 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'name': 'DAPI_max_proj', 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
]

layers_cells = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'name': 'DIC_masks', 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'name': 'DIC', 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'name': 'DAPI_masks', 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'name': 'DAPI', 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'name': 'DAPI_max_proj', 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
]

def import_layers(imagedir, mode=all, viewer=None):
    files = Path(imagedir).glob('*')
    stems = { Path(f).stem: f for f in files }
    logging.debug(stems)

    if mode == 'light':
        layers = layers_light
    elif mode == 'seg':
        layers = layers_segmentation
    elif mode == 'cells':
        layers = layers_cells
    else:
        layers = layers_ordered

    for layer in layers:
        logging.debug(f'layer: {layer}')
        if layer['name'] in stems.keys():
            logging.debug(f"stem: {stems[layer['name']]}")

            if viewer is not None:
                if layer['type'] == 'image':
                    viewer.add_image(io.imread(stems[layer['name']]), **layer['properties'])
                elif layer['type'] == 'labels':
                    viewer.add_labels(io.imread(stems[layer['name']]), **layer['properties'])
                elif layer['type'] == 'rna_filtered':
                    viewer.add_image(np.load(stems[layer['name']]), **layer['properties'])
                elif layer['type'] == 'spots':
                    viewer.add_points(np.load(Path(stems[layer['name']]).resolve()), **layer['properties'])

    if viewer is not None:
        viewer.title = Path(imagedir).parts[-1]
        for l in viewer.layers:
            l.scale = np.ones(len(l.scale))
        with open(Path(imagedir) / "img.json") as f:
            viewer.dims.set_point(0, json.load(f)['z_max_focus'])

if __name__ == '__main__':
    import_layers(imagedir=sys.argv[1], mode=sys.argv[2], viewer=napari.Viewer())