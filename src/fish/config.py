import logging
import json

class Config:
  def __init__(self, cfg_file):
    try:
      with open(cfg_file, 'r') as file:
        logging.info("cfg file exists")
        cfg = json.load(file)
        for key, value in cfg.items():
          setattr(self, key, value)
        # self.cfg = cfg
        self.cfg_file = cfg_file
        self.filter2mRNA = { v['filter']: v['mrna'] for v in cfg['channels'].values() }

        # self.mRNA2filter = {c['mrna']: c['filter'] for c in cfg['channels']}

    except FileNotFoundError:
      logging.warning("cfg file does not exist")





