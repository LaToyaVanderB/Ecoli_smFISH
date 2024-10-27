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
    { 'name': 'DIC_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image',
      'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image',
      'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image',
      'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },

    { 'name': 'rpoD', 'type': 'image',
      'properties': { 'colormap': 'magenta', 'visible': True, 'blending': 'additive' } },
    { 'name': 'rpoD_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'magenta', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rpoD_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'red', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'rnlAB', 'type': 'image',
      'properties': { 'colormap': 'cyan', 'visible': True, 'blending': 'additive' } },
    { 'name': 'rnlAB_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'cyan', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rnlAB_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'red', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'hipBA', 'type': 'image',
      'properties': { 'colormap': 'yellow', 'visible': True, 'blending': 'additive' } },
    { 'name': 'hipBA_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'yellow', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'hipBA_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'red', 'face_color': 'transparent', 'opacity': 0.5 } },
]

layers_spots = [
    { 'name': 'DIC_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image',
      'properties': { 'colormap': 'grey', 'visible': False, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image',
      'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image',
      'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },

    { 'name': 'rpoD', 'type': 'image',
      'properties': { 'colormap': 'magenta', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rpoD_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'magenta', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rpoD_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'red', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'rnlAB', 'type': 'image',
      'properties': { 'colormap': 'cyan', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rnlAB_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'cyan', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rnlAB_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'red', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'hipBA', 'type': 'image',
      'properties': { 'colormap': 'yellow', 'visible': False, 'blending': 'additive' } },
    { 'name': 'hipBA_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'yellow', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'hipBA_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'red', 'face_color': 'transparent', 'opacity': 0.5 } },
]


layers_input = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rpoD', 'type': 'image', 'properties': { 'colormap': 'magenta', 'visible': True, 'blending': 'additive' } },
    { 'name': 'rnlAB', 'type': 'image', 'properties': { 'colormap': 'cyan', 'visible': True, 'blending': 'additive' } },
    { 'name': 'hipBA', 'type': 'image', 'properties': { 'colormap': 'yellow', 'visible': True, 'blending': 'additive' } },
]

layers_light = [
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'translucent_no_depth', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
]

layers_segmentation = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
]

layers_cells = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'translate':  [-10, 10], 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
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
    elif mode == 'input':
        layers = layers_input
    elif mode == 'spots':
        layers = layers_spots
    else:
        layers = layers_ordered

    for layer in layers:
        logging.debug(f'layer: {layer}')
        if layer['name'] in stems.keys():
            logging.debug(f"stem: {stems[layer['name']]}")

            if viewer is not None:
                if layer['type'] == 'image':
                    l = viewer.add_image(io.imread(stems[layer['name']]), **layer['properties'])
                    l.name = layer['name']

                elif layer['type'] == 'labels':
                    l = viewer.add_labels(io.imread(stems[layer['name']]), **layer['properties'])
                    l.name = layer['name']

                elif layer['type'] == 'rna_filtered':
                    rna_filtered = np.load(stems[layer['name']])
                    l = viewer.add_image(rna_filtered, **layer['properties'])
                    l.name = layer['name']
                    l_mp = viewer.add_image(np.max(rna_filtered, axis=0), **layer['properties'])
                    l_mp.name = l.name + '_max_proj'

                elif layer['type'] == 'spots':
                    spots = np.load(Path(stems[layer['name']]).resolve())[:, 0:3]
                    l = viewer.add_points(spots, **layer['properties'])
                    l.name = layer['name']
                    l_mp = viewer.add_points(np.delete(spots, 0, 1), **layer['properties'])
                    l_mp.name = l.name + "_max_proj"
                    l_mp.border_color = 'blue'

    if viewer is not None:
        viewer.title = Path(imagedir).parts[-1]
        for l in viewer.layers:
            l.scale = np.ones(len(l.scale))
        with open(Path(imagedir) / "img.json") as f:
            viewer.dims.set_point(0, json.load(f)['z_max_focus'])


if __name__ == '__main__':
    import_layers(imagedir=sys.argv[1], mode=sys.argv[2], viewer=napari.Viewer())

