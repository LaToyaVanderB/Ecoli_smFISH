import logging
import json
from pathlib import Path
import re
from fish.image import Image


class Experiment:
    """
    An experiment is a list of images in the same directory,
    that share processing parameters
    """
    def __init__(self, cfg_file):
        try:
            with open(cfg_file, 'r') as file:
                logging.info(f'found cfg file {cfg_file}')
                cfg = json.load(file)
                for key, value in cfg.items():
                    setattr(self, key, value)
                self.cfg_file = cfg_file
                self.filter2mRNA = { v['filter']: v['mrna'] for v in cfg['channels'].values() }
                # deal with threshold 'None'
                for key in cfg['channels'].keys():
                    if 'threshold' not in cfg['channels'][key]:
                        self.channels[key]['threshold'] = None

        except FileNotFoundError:
            logging.warning(f'failed to find cfg file {cfg_file}')


    def create_image_list(self):
        self.images = {}

        # get the vsi files:
        vsi_files = list(Path(self.inputdir).glob('*.vsi'))

        for f in vsi_files:
            stem = re.sub(r'_CY[A-Z0-9\s,._]+DAPI', '', f.stem)
            params = { 'vsi_file': f.parts[-1], 'cell_file': None, 'valid': True }

            # we need a VSI directory
            if (Path(self.inputdir) / f'_{f.stem}_').is_dir():
                logging.info(f'found VSI directory _{f.stem}_')
            else:
                logging.warning(f'failed to find VSI directory _{f.stem}_')
                params['valid'] = False

            # we need a DIC file
            if (Path(self.inputdir) / f'{stem}_DIC.tif').is_file():
                params['cell_file'] = (Path(self.inputdir) / f'{stem}.tif').parts[-1]
                logging.info(f'found DIC file {params["cell_file"]}')
            else:
                logging.warning(f'failed to find DIC file {Path(self.inputdir)} / {stem}.tif')
                params['valid'] = False

            if params['valid'] == True:
                logging.info(f'{stem}: ok, adding image to image list')
                self.images[stem] = params
                # self.images[stem] = Image.from_dict(self, params)


    def read_image_list_from_jsons(self):
        # get the img.json files:
        self.json_files = sorted(Path(self.outputdir).glob('*/img.json'))

        # for f in json_files:
        #     self.images[f.parts[-1]] = Image.from_json(f, self)






