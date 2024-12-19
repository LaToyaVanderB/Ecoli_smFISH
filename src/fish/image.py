import logging
from pathlib import Path
from bioio import BioImage
import bioio_bioformats
from bioio.writers import OmeTiffWriter
from skimage import io
import numpy as np
import pandas as pd

from tools.utils import translate_image, remove_cell_masks_with_no_dapi, filter_cells_by_area, filter_cells_by_shape
import re
import json, jsonpickle

from omnipose.gpu import use_gpu
from cellpose_omni import io, transforms
from cellpose_omni import models

from skimage.measure import regionprops_table
from skimage.segmentation import expand_labels

# This checks to see if you have set up your GPU properly.
# CPU performance is a lot slower, but not a problem if you
# are only processing a few images.
use_GPU = use_gpu()

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p',
                    level=logging.INFO)


class Image:
    """
    Class for a FISH image.
    """

    def __init__(self, cfg, vsi_file, cell_file):
        """
        Class for a FISH image.

        Parameters
        ----------
        vsi_file : the VSI file with the fluorescent channels.
        cell_file : the file with the cell outline data, e.g. a 2D DIC picture.
        savedir : location to save the image.
        config : the configuration object.

        Returns
        -------
        a FISH image object.
        """

        self.cfg = cfg

        try:
            with open(Path(cfg.inputdir) / vsi_file, 'r') as file:
                self.vsi_file = vsi_file
                logging.info(f"found vsi file {vsi_file}")
        except FileNotFoundError:
            logging.warning(f"could not find vsi file {vsi_file}")

        try:
            with open(Path(cfg.inputdir) / cell_file, 'r') as file:
                self.cell_file = cell_file
                logging.info(f"found cell file {cell_file}")
        except FileNotFoundError:
            logging.warning(f"could not find cell file {cell_file}")

        try:
            stem = re.sub(r'_CY[A-Z0-9\s,._]+DAPI', '', Path(vsi_file).stem)
            self.stem = stem
            savepath = str(Path(cfg.outputdir) / stem)
            Path(savepath).mkdir(parents=True, exist_ok=True)
            self.savepath = savepath
            logging.info(f"created output dir {savepath}")
        except:
            logging.warning(f"failed to create output dir {savepath}")

        self.time = {}


    def read_image(self):
        img = {}
        logging.info(f'reading fluorescence image: {self.vsi_file}')
        image = BioImage(Path(self.cfg.inputdir) / self.vsi_file, reader=bioio_bioformats.Reader)
        logging.info(f'..found scenes: {image.scenes}')
        for s in image.scenes:
            if s == 'macro image':
                logging.info(f"....ignoring scene '{s}'")
            else:
                image.set_scene(s)
                logging.info(f"....reading scene '{s}' (shape {image.shape})")
                logging.info(f"....found channels {image.channel_names}")

                for ch in enumerate(s.replace(r'001 ', '').split(", ")):
                    if ch[1] in self.cfg.filter2mRNA.keys():
                        filter = ch[1]
                        mrna = self.cfg.filter2mRNA[ch[1]]
                        img[mrna] = {}
                        logging.info(f"......reading in channel: C={ch[0]} filter={filter} mrna={mrna}")
                        img[mrna]['data'] = image.get_image_data("ZYX", C=ch[0])
                        self.mrna = img

                    else:
                        logging.warning(f"......ignoring unknown channel {filter}")


    def read_cells(self):
        logging.info(f'reading cell image: {self.cell_file}')
        cells = io.imread(Path(self.cfg.inputdir) / self.cell_file)
        self.cells = {}
        self.cells['data'] = cells


    def crop(self, y1, x1, y2, x2):
        self.cropped = ''


    def align(self):
        logging.info(f"aligning DIC image by: {self.cfg.translate_dic_by}")
        dy = self.cfg.translate_dic_by[0]
        dx = self.cfg.translate_dic_by[1]

        # translate and crop DIC
        img_translated = translate_image(self.cells['data'], dy, dx)
        self.cells['aligned'] = img_translated[1]

        # crop fluorescent channels
        for ch in self.mrna.keys():
            data = self.mrna[ch]['data']
            Y, X = data.shape[1], data.shape[2]
            self.mrna[ch]['aligned'] = data[:, max(dy, 0):Y + min(dy, 0), max(dx, 0):X + min(dx, 0)]


    def create_grgb(self):
        logging.info(f'creating GRGB composite image')
        cells = self.cells['aligned']
        dapimaxproj = np.max(self.mrna['DAPI']['aligned'], axis=0)
        self.grgb = np.stack([np.zeros(cells.shape), cells, dapimaxproj, np.zeros(cells.shape)], axis=0)


    def save_layers(self):
        logging.info(f'saving layers to: {self.savepath}')
        # fluorescent channels
        for ch in self.mrna.keys():
            data = self.mrna[ch]['aligned']
            savepath = Path(self.savepath) / f'{ch}.tif'
            OmeTiffWriter.save(data, savepath, dim_order="ZYX", channel_names=[ch])
            logging.info(f'..saving {ch} channel to {savepath}')

        # cell channel
        data = self.cells['aligned']
        savepath = Path(self.savepath) / 'DIC.tif'
        io.imwrite(savepath, data)
        logging.info(f'..saving DIC channel to {savepath}')

        # grgb file for segmentation
        np.save(Path(self.savepath) / 'grgb.npy', self.grgb)


    def save_metadata(self):
        logging.info(f'saving metadata to: {self.savepath}/img.json')
        keys = [
            'vsi_file',
            'cell_file',
            'stem',
            'savepath',
            'segmentation',
            'time'
        ]
        metadata = { i: getattr(self, i) for i in keys if hasattr(self, i)}

        with open(Path(self.savepath) / "img.json", "w") as f:
            json.dump(metadata, f, default=vars, indent=4)


    def save(self, suffix):
        logging.info(f'saving image to: {self.savepath}/img.{suffix}.pkl')
        with open(Path(self.savepath) / f'img.{suffix}.pkl', "w") as f:
            f.write(jsonpickle.encode(self))


    def unpickle(pickle_file):
        with open(pickle_file, "r+") as f:
            return jsonpickle.decode(f.read())



    def segment(self, **params):
        self.cell_masks = ''
        self.dapi_masks = ''
        self.params = params
        self.regionprops = ''

        self.segment_cells()
        self.segment_dapi()

    # redo for list of images ot avoid loading the models for each image
    def segment_cells(self):
        logging.info(f'segmenting DIC image')

        # pass parameters to model
        params = {
            'channels': [1, 2],  # always define this with the model
            'rescale': None,  # upscale or downscale your images, None = no rescaling
            'mask_threshold': 0,  # erode or dilate masks with higher or lower values between -5 and 5
            'flow_threshold': 0,
            'min_size': 200,
            'diameter': 0,
            'invert': False,
            'transparency': True,  # transparency in flow output
            'omni': True,  # we can turn off Omnipose mask reconstruction, not advised
            'cluster': True,  # use DBSCAN clustering
            'resample': True,  # whether or not to run dynamics on rescaled grid or original grid
            'verbose': False,  # turn on if you want to see more output
            'tile': False,  # average the outputs from flipped (augmented) images; slower, usually not needed
            'niter': None,
            # default None lets Omnipose calculate # of Euler iterations (usually <20) but you can tune it for over/under segmentation
            'augment': False,  # Can optionally rotate the image and average network outputs, usually not needed
            'affinity_seg': False,  # new feature, stay tuned...
        }

        # best models for DIC images
        model1 = models.CellposeModel(gpu=use_GPU, model_type="cyto2")
        # possible channels options for model1: [2,0] (no nuclei)  or [2,3] (nuclei)

        model2 = models.CellposeModel(gpu=use_GPU, model_type="cyto2_omni")
        # possible channels options for model2: [1,2] (nuclei)

        model_menu = {
            "cyto2": model1,
            "cyto2_omni": model2
        }

        # default
        model_type = "cyto2_omni"
        model = model_menu[model_type]

        # overwrite default segmentation parameters with img.json values
        if hasattr(self, 'segmentation'):
            params['diameter'] = self.segmentation['cells']['diameter']
            params['mask_threshold'] = self.segmentation['cells']['mask']
            params['flow_threshold'] = self.segmentation['cells']['flow']
            model = model_menu[self.segmentation['cells']['model_type']]
            params['channels'] = self.segmentation['cells']['channels']
            logging.info(f"model_type from config-selected: {self.segmentation['cells']['model_type']} {params['channels']}")
        else:
            self.segmentation = { 'cells': {}, 'dapi': {} }

        mask, flow, style = model.eval(self.grgb, **params)
        self.cell_masks = mask
        self.segmentation['cells'] = params
        maskfile_latest = Path(self.savepath) / f'DIC_masks.model={model_type}_chan={str(params["channels"]).replace(" ", "")}_diameter={params["diameter"]}_minsize={params["min_size"]}_mask={params["mask_threshold"]}_flow={params["flow_threshold"]}.tif'
        io.imwrite(maskfile_latest, mask)
        cellmaskfile = Path(self.savepath) / 'DIC_masks'
        cellmaskfile.unlink(missing_ok=True)
        cellmaskfile.symlink_to(maskfile_latest.parts[-1])
        logging.info(f"writing cell mask to {maskfile_latest}")


    def segment_dapi(self):
        logging.info(f'segmenting DAPI image')

        # model for nuclei
        model_type = "nuclei"
        model = models.CellposeModel(gpu=use_GPU, model_type=model_type)

        # pass parameters to model
        params = {
            'channels': [0, 0], # always define this with the model
            'rescale': None, # upscale or downscale your images, None = no rescaling
            'mask_threshold': 0.0, # erode or dilate masks with higher or lower values between -5 and 5
            'flow_threshold': 0.0,
            'min_size': 10,
            'diameter': 0.0,
            'invert': False,
            'transparency': True, # transparency in flow output
            'omni': True, # we can turn off Omnipose mask reconstruction, not advised
            'cluster': True, # use DBSCAN clustering
            'resample': True, # whether or not to run dynamics on rescaled grid or original grid
            'verbose': False, # turn on if you want to see more output
            'tile': False, # average the outputs from flipped (augmented) images; slower, usually not needed
            'niter': None, # default None lets Omnipose calculate # of Euler iterations (usually <20) but you can tune it for over/under segmentation
            'augment': False, # Can optionally rotate the image and average network outputs, usually not needed
            'affinity_seg': False, # new feature, stay tuned...
        }

        mask, flow, style = model.eval(self.grgb[2, ...], **params)
        self.dapi_masks = mask
        self.segmentation['dapi'] = params

        maskfile_latest = Path(self.savepath) / f'DAPI_masks_model={model_type}_chan={str(params["channels"]).replace(" ", "")}_diameter={params["diameter"]}_minsize={params["min_size"]}_mask={params["mask_threshold"]}_flow={params["flow_threshold"]}.tif'
        io.imwrite(maskfile_latest, mask)
        cellmaskfile = Path(self.savepath) / 'DAPI_masks'
        cellmaskfile.unlink(missing_ok=True)
        cellmaskfile.symlink_to(maskfile_latest.parts[-1])
        logging.info(f"writing DAPI mask to {maskfile_latest}")


    def postprocess_masks(self):
        logging.info(f'postprocessing masks')

        logging.info(f'..removing cell masks with no DAPI')
        self.cell_masks_with_dapi = remove_cell_masks_with_no_dapi(self.cell_masks, self.dapi_masks)[0]

        logging.info(f'..discarding masks of size < 200')
        self.regionprops = pd.DataFrame(regionprops_table(self.cell_masks, properties=['label', 'bbox', 'area', 'eccentricity', 'centroid', 'axis_major_length', 'axis_minor_length', 'orientation']))
        self.cell_masks_by_area = filter_cells_by_area(self.cell_masks_with_dapi, self.regionprops, min_cell_area=200)[0]

        logging.info(f'..discarding large, round masks (min size 1000, min excentricity 0.8')
        self.cell_masks_by_shape = filter_cells_by_shape(self.cell_masks_by_area, self.regionprops,
                                                         min_clump_area=1000, max_clump_eccentricity=0.8)[0]

    def expand_masks(self, nr_of_pixels):
        logging.info(f'expanding masks by {nr_of_pixels} pixels')
        self.cell_masks_expanded = expand_labels(masks, distance=distance)

    def find_focus(self):
        self.focus = ''


    def filter(self, sigma):
        self.filtered = ''


    def detect_spots(self):
        self.spots = ''


    def decompose_spots(self):
        self.decomposed_spots = ''


    def count_spots(self):
        self.counts = ''
        self.histograms = ''


    def view(self, mode):
        self.viewer = ''





    def delete(self):
        return True


