import sys
import logging
from copy import copy
from pathlib import Path
from skimage import io
import numpy as np
import json
import napari
import re
import pandas as pd
from skimage.segmentation import expand_labels
from skimage.measure import regionprops_table

# need to do this in napari console otherwise the import does not work
# import sys; sys.path.insert(0, "")

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)



layers_ordered = [
    { 'name': 'DIC_masks_selected_by_shape', 'type': 'labels',
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
    { 'name': 'rpoD_decomposed_spots', 'type': 'spots',
     'properties': {'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc',
                    'size': 5, 'border_width': 0.1, 'border_color': 'magenta', 'face_color': 'transparent',
                    'opacity': 1}},
    { 'name': 'rpoD_ddregions', 'type': 'spots',
     'properties': {'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc',
                    'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'cyan',
                    'opacity': 0.5}},

    { 'name': 'rnlAB', 'type': 'image',
      'properties': { 'colormap': 'cyan', 'visible': False, 'blending': 'additive' } },
    { 'name': 'rnlAB_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'cyan', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'rnlAB_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'magenta', 'face_color': 'transparent', 'opacity': 0.5 } },
    { 'name': 'rnlAB_decomposed_spots', 'type': 'spots',
      'properties': {'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc',
                    'size': 5, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent',
                    'opacity': 1}},
    { 'name': 'rnlAB_ddregions', 'type': 'spots',
      'properties': {'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc',
                    'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'cyan',
                    'opacity': 0.5}},

    { 'name': 'hipBA', 'type': 'image',
      'properties': { 'colormap': 'yellow', 'visible': False, 'blending': 'additive' } },
    { 'name': 'hipBA_filtered_padded', 'type': 'rna_filtered',
      'properties': { 'colormap': 'yellow', 'blending': 'additive', 'visible': False, 'contrast_limits': [0, 6000] } },
    { 'name': 'hipBA_spots', 'type': 'spots',
      'properties': { 'blending': 'additive', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 } },
    {'name': 'hipBA_decomposed_spots', 'type': 'spots',
     'properties': {'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc',
                    'size': 5, 'border_width': 0.1, 'border_color': 'yellow', 'face_color': 'transparent',
                    'opacity': 1}},
    {'name': 'hipBA_ddregions', 'type': 'spots',
     'properties': {'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc',
                    'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'cyan',
                    'opacity': 0.5}},
]

layers_input = [
    { 'name': 'DIC_masks_selected_by_shape', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
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

layers_spots = [
    { 'name': 'DIC_masks_selected_by_shape', 'type': 'labels',
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
    { 'name': 'DIC_masks_selected_by_shape', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
    { 'name': 'DAPI_max_proj', 'type': 'image', 'properties': { 'colormap': 'blue', 'visible': True, 'blending': 'additive' } },
    { 'name': 'DAPI_masks', 'type': 'labels', 'properties': { 'visible': False, 'blending': 'additive', 'opacity': 0.2 } },
]

spot_filename_pattern = re.compile(rf'_spots_thr=(?P<threshold>[0-9.]+)_ifx1=(?P<ifx1>[0-9]+)_ifx2=(?P<ifx2>[0-9]+).npy')

mode_layers = {
    'all': layers_ordered,
    'light': layers_light,
    'cells': layers_cells,
    'input': layers_input,
    'spots': layers_spots,
}

# this is horrendous, do it properly later
def import_layers_world(imagedir, mode='all', viewer=None):
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
            crop[0] = 0.065 * crop[0]
            crop[1] = 0.065 * crop[1]

    if viewer is not None:
        for l in viewer.layers:
            l.visible = False

    for layer in layers:
        logging.debug(f'layer: {layer}')
        if layer['name'] in stems.keys():
            logging.debug(f"stem: {stems[layer['name']]}")

            if viewer is not None:

                if layer['type'] == 'image':
                    img_data = io.imread(stems[layer['name']])
                    if img_data.ndim == 3:
                        layer['properties']['scale'] = [0.2, 0.065, 0.065]
                    elif img_data.ndim == 2:
                        layer['properties']['scale'] = [0.065, 0.065]
                    l = viewer.add_image(img_data, **layer['properties'])
                    l.name = layer['name']

                elif layer['type'] == 'labels':
                    layer['properties']['scale'] = [0.065, 0.065]
                    masks = io.imread(stems[layer['name']])
                    if 'DIC' in layer['name']:
                        props = pd.DataFrame(regionprops_table(masks, properties=['label', 'area', 'eccentricity', 'perimeter', 'solidity', 'axis_minor_length', 'orientation', 'axis_major_length']))
                        props = props.set_index('label', drop=False).set_index('label', drop=False)
                        features = pd.DataFrame(index=np.arange(np.max(masks)))
                        features = pd.concat([features, props], axis=1).convert_dtypes().round(3)
                        layer['properties']['features'] = features

                    l = viewer.add_labels(masks, **layer['properties'])
                    l.name = layer['name']

                elif layer['type'] == 'rna_filtered':
                    layer['properties']['scale'] = [0.2, 0.065, 0.065]

                    rna_filtered = np.load(stems[layer['name']])

                    if mode != 'spots':
                        l = viewer.add_image(rna_filtered, **layer['properties'])
                        l.name = layer['name']
                        l.translate = [0, crop[0], crop[1]]

                    layer['properties']['scale'] = [0.065, 0.065]
                    l_mp = viewer.add_image(np.max(rna_filtered, axis=0), **layer['properties'])
                    l_mp.name = layer['name'] + '_max_proj'
                    l_mp.translate = [crop[0], crop[1]]

                elif layer['type'] == 'spots':
                    spots_all = np.load(stems[layer['name']].resolve())
                    # print(stems[layer['name']].resolve())
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
                    # get cell features as well

                    l_name = layer['name']
                    match = spot_filename_pattern.search(str(stems[layer['name']].resolve()))
                    if match is not None:
                        d = match.groupdict()
                        l_name = f'{layer["name"]} thr={float(d["threshold"]):.2f} [{d["ifx1"]}, {d["ifx2"]}]'

                    if mode != 'spots':
                        layer['properties']['scale'] = [0.2, 0.065, 0.065]
                        l = viewer.add_points(spots, **layer['properties'])
                        # print(l.data.shape)
                        l.name = l_name
                        l.translate = [0, crop[0], crop[1]]

                    layer['properties']['scale'] = [0.065, 0.065]
                    l_mp = viewer.add_points(np.delete(spots, 0, 1), **layer['properties'])
                    l_mp.name = l_name + "_max_proj"
                    l_mp.border_width = 0.05
                    l_mp.translate = [crop[0], crop[1]]



    if viewer is not None:
        viewer.title = Path(imagedir).parts[-1]
        # for l in viewer.layers:
        #     l.scale = np.ones(len(l.scale))
        #     print(l.name, l.scale, l.translate)
        for l in viewer.layers:
            print(l.scale, l.translate)
        #     l.scale = 0.065 * l.scale
        #     l.translate = 0.065 * l.translate

        with open(Path(imagedir) / "img.json") as f:
            img_json = json.load(f)
            focus = 20
            if 'rpoD' in img_json and 'z_max_focus' in img_json['rpoD']:
                focus = img_json['rpoD']['z_max_focus']
            viewer.dims.set_point(0, focus)


def import_layers(imagedir, mode='all', viewer=None):
    files = Path(imagedir).glob('*')
    stems = { Path(f).stem: f for f in files }
    logging.debug(stems)

    layers = mode_layers[mode] if mode in mode_layers else layers_ordered

    # crop = False
    # with open(Path(imagedir) / "img.json") as f:
    #     img_json = json.load(f)
    #     crop = [0, 0, 0, 0]
    #     if 'crop' in img_json:
    #         crop = img_json['crop']
    #         crop[0] = 0.065 * crop[0]
    #         crop[1] = 0.065 * crop[1]

    if viewer is not None:
        for l in viewer.layers:
            l.visible = False

    for layer in layers:
        logging.debug(f'layer: {layer}')
        if layer['name'] in stems.keys():
            logging.debug(f"stem: {stems[layer['name']]}")

            if viewer is not None:

                if layer['type'] == 'image':
                    img_data = io.imread(stems[layer['name']])
                    # if img_data.ndim == 3:
                    #     layer['properties']['scale'] = [0.2, 0.065, 0.065]
                    # elif img_data.ndim == 2:
                    #     layer['properties']['scale'] = [0.065, 0.065]
                    l = viewer.add_image(img_data, **layer['properties'])
                    l.name = layer['name']

                elif layer['type'] == 'labels':
                    # layer['properties']['scale'] = [0.065, 0.065]
                    masks = io.imread(stems[layer['name']])
                    if 'DIC' in layer['name']:
                        props = pd.DataFrame(regionprops_table(masks, properties=['label', 'area', 'eccentricity', 'perimeter', 'solidity', 'axis_minor_length', 'orientation', 'axis_major_length']))
                        props = props.set_index('label', drop=False).set_index('label', drop=False)
                        features = pd.DataFrame(index=np.arange(np.max(masks)))
                        features = pd.concat([features, props], axis=1).convert_dtypes().round(3)
                        layer['properties']['features'] = features

                    l = viewer.add_labels(masks, **layer['properties'])
                    l.name = layer['name']

                elif layer['type'] == 'rna_filtered':
                    # layer['properties']['scale'] = [0.2, 0.065, 0.065]

                    rna_filtered = np.load(stems[layer['name']])

                    if mode != 'spots':
                        l = viewer.add_image(rna_filtered, **layer['properties'])
                        l.name = layer['name']
                        # l.translate = [0, crop[0], crop[1]]

                    # layer['properties']['scale'] = [0.065, 0.065]
                    l_mp = viewer.add_image(np.max(rna_filtered, axis=0), **layer['properties'])
                    l_mp.name = layer['name'] + '_max_proj'
                    # l_mp.translate = [crop[0], crop[1]]

                elif layer['type'] == 'spots':
                    spots_all = np.load(stems[layer['name']].resolve())
                    # print(stems[layer['name']].resolve())
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
                    # get cell features as well

                    l_name = layer['name']
                    match = spot_filename_pattern.search(str(stems[layer['name']].resolve()))
                    if match is not None:
                        d = match.groupdict()
                        l_name = f'{layer["name"]} thr={float(d["threshold"]):.2f} [{d["ifx1"]}, {d["ifx2"]}]'

                    if mode != 'spots':
                        # layer['properties']['scale'] = [0.2, 0.065, 0.065]
                        l = viewer.add_points(spots, **layer['properties'])
                        # print(l.data.shape)
                        l.name = l_name
                        # l.translate = [0, crop[0], crop[1]]

                    # layer['properties']['scale'] = [0.065, 0.065]
                    l_mp = viewer.add_points(np.delete(spots, 0, 1), **layer['properties'])
                    l_mp.name = l_name + "_max_proj"
                    l_mp.border_width = 0.05
                    # l_mp.translate = [crop[0], crop[1]]



    if viewer is not None:
        viewer.title = Path(imagedir).parts[-1]
        # for l in viewer.layers:
        #     l.scale = np.ones(len(l.scale))
        #     print(l.name, l.scale, l.translate)
        for l in viewer.layers:
            print(l.scale, l.translate)
        #     l.scale = 0.065 * l.scale
        #     l.translate = 0.065 * l.translate

        with open(Path(imagedir) / "img.json") as f:
            img_json = json.load(f)
            focus = 20
            if 'rpoD' in img_json and 'z_max_focus' in img_json['rpoD']:
                focus = img_json['rpoD']['z_max_focus']
            viewer.dims.set_point(0, focus)


def expand_masks(distance, viewer):
    mask_layer_data = copy(viewer.layers['DIC_masks'].data)
    mask_layer_name = viewer.layers['DIC_masks'].name
    expanded = expand_labels(mask_layer_data, distance)
    expanded_mask_layer = viewer.add_labels(expanded, name=f'{mask_layer_name}_expanded_{distance}px')
    expanded_mask_layer.contour = 1
    # need to crop the masks here and translate them


    for mrna in ['rpoD', 'rnlAB', 'hipBA']:
        spots = [l for l in viewer.layers if re.match(f'{mrna}_spots.*max_proj', l.name)][0]
        print(spots, spots.translate[0], spots.translate[1])
        expanded_spots_data = copy(spots.data)
        expanded_spots_properties = { 'blending': 'translucent_no_depth', 'visible': False, 'out_of_slice_display': True, 'symbol': 'disc', 'size': 10, 'border_width': 0.1, 'border_color': 'cyan', 'face_color': 'transparent', 'opacity': 0.5 }
        expanded_spots_properties['features'] = copy(spots.features)
        expanded_spots_properties['features']['label'] = [expanded_mask_layer.data[y+int(spots.translate[0]), x+int(spots.translate[1])] for (y, x) in expanded_spots_data]
        expanded_spots_properties['features']['in_cell'] = expanded_spots_properties['features'].apply(lambda s: False if s['label'] == 0 else True, axis=1)
        expanded_spots_properties['border_color'] = 'in_cell'
        expanded_spots_properties['border_color_cycle'] = ['cyan', 'red'] if expanded_spots_properties['features'].iloc[0]['in_cell'] == True else ['red', 'cyan']
        expanded_spots_properties['name'] = f'{mrna} spots after {distance}px expansion'
        # print(f'data={expanded_spots_data}')
        # print(f'expanded_spots_properties={expanded_spots_properties}')
        # print(type(expanded_spots_properties['features']))

        p = viewer.add_points(expanded_spots_data, **expanded_spots_properties)
        p.translate = spots.translate

        viewer.layers['DIC_masks'].visible = False

if __name__ == '__main__':
    import_layers(imagedir=sys.argv[1], mode=sys.argv[2], viewer=napari.Viewer())

