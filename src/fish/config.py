import logging
import json

class Config:
  def __init__(self, cfg_file):
    try:
      with open(cfg_file, 'r') as file:
        logging.info(f'found cfg file {cfg_file}')
        cfg = json.load(file)
        for key, value in cfg.items():
          setattr(self, key, value)
        # self.cfg = cfg
        self.cfg_file = cfg_file
        self.filter2mRNA = { v['filter']: v['mrna'] for v in cfg['channels'].values() }
        # deal with threshold 'None'

    except FileNotFoundError:
      logging.warning(f'failed to find cfg file {cfg_file}')





