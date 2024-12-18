import logging
from pathlib import Path
from bioio import BioImage
from bioio.writers import OmeTiffWriter
from skimage import io
import numpy as np

from tools.utils import translate_image
import re
import json



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

    try:
      with open(Path(cfg.cfg['inputdir']) / vsi_file, 'r') as file:
        cfg.vsi_file = vsi_file
        logging.info(f"found vsi file {cfg.vsi_file}")
    except FileNotFoundError:
      logging.warning(f"could not find vsi file {cfg.vsi_file}")

    try:
      with open(Path(cfg.cfg['inputdir']) / cell_file, 'r') as file:
        cfg.cell_file = cell_file
        logging.info(f"found cell file {cfg.cell_file}")
    except FileNotFoundError:
      logging.warning(f"could not find cell file {cfg.cell_file}")

    try:
      cfg.outputdir = Path(cfg.cfg['outputdir'])
      cfg.stem = Path(vsi_file).stem
      cfg.stem = re.sub(r'_CY[A-Z0-9\s,._]+DAPI', '', cfg.stem)
      cfg.savepath = cfg.outputdir / cfg.stem
      Path(cfg.savepath).mkdir(parents=True, exist_ok=True)
      logging.info(f"created output dir {cfg.savepath}")
    except:
      logging.warning(f"failed to create output dir {cfg.savepath}")

    self.metadata = cfg


  def read_image(self):
    img = {}
    logging.info(f'reading fluorescence image: {self.metadata.vsi_file}')
    image = BioImage(Path(self.metadata.cfg['inputdir']) / self.metadata.vsi_file)
    logging.info(f'..found scenes: {image.scenes}')
    for s in image.scenes:
      if s == 'macro image':
        logging.info(f"....ignoring scene '{s}'")
      else:
        image.set_scene(s)
        logging.info(f"....reading scene '{s}' (shape {image.shape})")
        logging.info(f"....found channels {image.channel_names}")

        for ch in enumerate(s.replace(r'001 ', '').split(", ")):
          if ch[1] in self.metadata.filter2mRNA.keys():
            filter = ch[1]
            mrna = self.metadata.filter2mRNA[ch[1]]
            img[mrna] = {}
            logging.info(f"......reading in channel: C={ch[0]} filter={filter} mrna={mrna}")
            img[mrna]['data'] = image.get_image_data("ZYX", C=ch[0])
            self.mrna = img

          else:
            logging.warning(f"......ignoring unknown channel {filter}")


  def read_cells(self):
    logging.info(f'reading cell image: {self.metadata.cell_file}')
    cells = io.imread(Path(self.metadata.cfg['inputdir']) / self.metadata.cell_file)
    self.cells = {}
    self.cells['data'] = cells


  def crop(self, y1, x1, y2, x2):
      self.cropped = ''


  def align(self):
    logging.info(f'aligning DIC image by: {self.metadata.cfg['translate_dic_by']}')
    dy = self.metadata.cfg['translate_dic_by'][0]
    dx = self.metadata.cfg['translate_dic_by'][1]

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
    logging.info(f'saving layers to: {self.metadata.savepath}')
    # fluorescent channels
    for ch in self.mrna.keys():
      data = self.mrna[ch]['aligned']
      savepath = self.metadata.savepath / f'{ch}.tif'
      OmeTiffWriter.save(data, savepath, dim_order="ZYX", channel_names=[ch])
      logging.info(f'..saving {ch} channel to {savepath}')

    # cell channel
    data = self.cells['aligned']
    savepath = self.metadata.savepath / 'DIC.tif'
    io.imsave(savepath, data)
    logging.info(f'..saving DIC channel to {savepath}')

    # grgb file for segmentation
    np.save(self.metadata.savepath / 'grgb.npy', self.grgb)


  def save_metadata(self):
    # serialize config
    self.metadata.savepath = str(self.metadata.savepath)
    self.metadata.outputdir = str(self.metadata.outputdir)
    del self.metadata.cfg
    del self.metadata.filter2mRNA
    del self.metadata.mRNA2filter
    json.dump(vars(self.metadata), open(Path(self.metadata.savepath) / "img.json", "w"), indent=4)


  def segment(self, **params):
    self.cell_masks = ''
    self.dapi_masks = ''
    self.params = params
    self.regionprops = ''


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


