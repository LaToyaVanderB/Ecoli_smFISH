import sys
import logging
from glob import glob
from pathlib import Path
from skimage import io
import numpy as np
import json
import napari
import re
import pandas as pd

# need to do this in napari console otherwise the import does not work
# import sys; sys.path.insert(0, "")

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)



layers_ordered = [
    { 'name': 'DIC_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image',
      'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive' } },
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
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'rnlAB', 'type': 'image',
      'properties': { 'colormap': 'cyan', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rnlAB_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'cyan', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rnlAB_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'hipBA', 'type': 'image',
      'properties': { 'colormap': 'yellow', 'visible': False, 'blending': 'additive' } },
    { 'name': 'hipBA_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'yellow', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'hipBA_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },
]

layers_spots = [
    { 'name': 'DIC_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image',
      'properties': { 'colormap': 'grey', 'visible': False, 'blending': 'additive' } },
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
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'rnlAB', 'type': 'image',
      'properties': { 'colormap': 'cyan', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rnlAB_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'cyan', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rnlAB_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'hipBA', 'type': 'image',
      'properties': { 'colormap': 'yellow', 'visible': False, 'blending': 'additive' } },
    { 'name': 'hipBA_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'yellow', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'hipBA_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },
]

layers_input = [
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DIC_untranslated', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': False, 'blending': 'additive', 'translate':  [-10, 10] } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rpoD', 'type': 'image', 'properties': { 'colormap': 'magenta', 'visible': True, 'blending': 'additive' } },
    { 'name': 'rnlAB', 'type': 'image', 'properties': { 'colormap': 'cyan', 'visible': True, 'blending': 'additive' } },
    { 'name': 'hipBA', 'type': 'image', 'properties': { 'colormap': 'yellow', 'visible': True, 'blending': 'additive' } },
]

layers_light = [
    { 'name': 'DIC_untranslated', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': False, 'blending': 'translucent_no_depth', 'translate':  [-10, 10] } },
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
]

layers_spotlight = [
    { 'name': 'DIC_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DIC', 'type': 'image',
      'properties': { 'colormap': 'grey', 'visible': False, 'blending': 'additive' } },
    { 'name': 'DAPI_masks', 'type': 'labels',
      'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI_max_proj', 'type': 'image',
      'properties': { 'colormap': 'blue', 'visible': False, 'blending': 'additive' } },

    { 'name': 'rpoD_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'magenta', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rpoD_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'rnlAB_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'cyan', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rnlAB_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },

    { 'name': 'hipBA_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'yellow', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'hipBA_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },
]

layers_cells = [
    { 'name': 'DIC', 'type': 'image', 'properties': { 'colormap': 'grey', 'visible': True, 'blending': 'additive'} },
    { 'name': 'DIC_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
]

spot_filename_pattern = re.compile(rf'_spots_thr=(?P<threshold>[0-9.]+)_ifx1=(?P<ifx1>[0-9]+)_ifx2=(?P<ifx2>[0-9]+).npy')

mode_layers = {
    'all': layers_ordered,
    'light': layers_light,
    'cells': layers_cells,
    'input': layers_input,
    'spots': layers_spots
}

def import_layers(imagedir, mode=all, viewer=None):
    files = Path(imagedir).glob('*')
    stems = { Path(f).stem: f for f in files }
    logging.debug(stems)

    layers = mode_layers[mode] if mode in mode_layers else layers_ordered

    crop = False
    with open(Path(imagedir) / "img.json") as f:
        img_json = json.load(f)
        crop = [0, 0, 0, 0]
        if 'crop' in img_json:
            crop = img_json['crop']

    if viewer is not None:
        for l in viewer.layers:
            l.visible = False

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

                    if mode != 'spotlight':
                        l = viewer.add_image(rna_filtered, **layer['properties'])
                        l.name = layer['name']
                        l.translate = [0, crop[0], crop[1]]

                    l_mp = viewer.add_image(np.max(rna_filtered, axis=0), **layer['properties'])
                    l_mp.name = layer['name'] + '_max_proj'
                    l_mp.translate = [crop[0], crop[1]]

                elif layer['type'] == 'spots':
                    spots_all = np.load(stems[layer['name']].resolve())
                    print(stems[layer['name']].resolve())
                    spots = spots_all[:, 0:3]

                    spot_features = pd.DataFrame()
                    if spots_all.shape[1] == 5:
                        spot_features = pd.DataFrame(spots_all, columns=['z', 'y', 'x', 'intensity', 'filtered_intensity'])
                    elif spots_all.shape[1] == 6:
                        spot_features = pd.DataFrame(spots_all, columns=['z', 'y', 'x', 'intensity', 'filtered_intensity', 'label'])
                        spot_features['in_cell'] = spot_features.apply(lambda s: False if s['label'] == 0 else True, axis=1)
                        layer['properties']['border_color'] = 'in_cell'
                        layer['properties']['border_color_cycle'] = ['cyan', 'red'] if spot_features.iloc[0]['in_cell'] == True else ['red', 'cyan']
                    layer['properties']['features'] = spot_features

                    print(spots_all.shape, spot_features.shape)

                    match = spot_filename_pattern.search(str(stems[layer['name']].resolve()))
                    d = match.groupdict()
                    l_name = f'{layer["name"]} thr={float(d["threshold"]):.2f} [{d["ifx1"]}, {d["ifx2"]}]'

                    if mode != 'spotlight':
                        l = viewer.add_points(spots, **layer['properties'])
                        print(l.data.shape)
                        l.name = l_name
                        l.translate = [0, crop[0], crop[1]]

                    l_mp = viewer.add_points(np.delete(spots, 0, 1), **layer['properties'])
                    l_mp.name = l_name + "_max_proj"
                    l_mp.border_width = 0.05
                    l_mp.translate = [crop[0], crop[1]]



    if viewer is not None:
        viewer.title = Path(imagedir).parts[-1]
        for l in viewer.layers:
            l.scale = np.ones(len(l.scale))
        with open(Path(imagedir) / "img.json") as f:
            img_json = json.load(f)
            focus = 20
            if 'rpoD' in img_json and 'z_max_focus' in img_json['rpoD']:
                focus = img_json['rpoD']['z_max_focus']
            viewer.dims.set_point(0, focus)


if __name__ == '__main__':
    import_layers(imagedir=sys.argv[1], mode=sys.argv[2], viewer=napari.Viewer())

