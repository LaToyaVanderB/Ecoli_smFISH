from glob import glob
import os.path
import sys
import re
import shutil
from pathlib import Path
import json
import time
import numpy as np
from bioio import BioImage
import bioio_bioformats
from bioio.writers import OmeTiffWriter
import argparse

import logging
logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s ', datefmt='%m/%d/%Y %I:%M:%S%p', level=logging.INFO)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file (process multiple images).")
    parser.add_argument("-only", "--only", help="VSI image file (process only one image. Requires a config file anyway).")
    args = parser.parse_args()
    logging.info(f'{args}')
    configfile = args.config
    only = args.only

    tic0 = time.time()
    logging.info(f'reading config file {configfile}.')

    # configfile contains:
    # - expected directory structure
    # - experiments to analyse
    # - expected channels (filter x mRNA)
    # all experiments are assumed to contain the same channels
    with open(configfile, 'r') as f:
        config = json.load(f)

    # channels
    filter2mRNA = {c['filter']: c['mrna'] for c in config['channels']}
    mRNA2filter = {c['mrna']: c['filter'] for c in config['channels']}

    # find source images for each experiment:
    for exp in config['experiments']:
        exp['images'] = [{'sourcefile': f} for f in glob(f"{config['inputdir']}/{exp['strain']}_{exp['condition']}*") if f.find(r'_DIC') == -1]

    # how many images do we have
    config['nr_images'] = sum([len(exp['images']) for exp in config['experiments']])

    # parse source file name
    # look for DIC image
    n = 0
    for exp in config['experiments']:
        pattern = re.compile(rf'({exp["strain"]}_{exp["condition"]}.*)' + '_' + r'(\d+)')
        for img in exp['images']:
            n = n + 1
            logging.info(f"parsing filename: {img['sourcefile']} [{n}/{config['nr_images']}]")
            img['rootdir'], img['basename'] = os.path.split(img['sourcefile'])
            img['basename'], img['format'] = os.path.splitext(img['basename'])
            img['format'] = img['format'][1:]
            logging.info(f"basename: {img['basename']}, format: {img['format']}")

            img['time'] = {}

            img['stem'] = re.sub(r'_CY[A-Z0-9\s,._]+DAPI', '', img['basename'])
            _, img['seqnr'] = re.match(pattern, img['stem']).groups()
            logging.info(f"stem: {img['stem']}, seqnr: {img['seqnr']}")
            os.makedirs(os.path.join(config['outputdir'], img['stem']), exist_ok=True)

            img['dicbasename'] = re.sub(r'_CY[A-Z0-9\s,.]+DAPI', '_DIC', img['basename']) + ".tif"
            logging.info(f"dicbasename: {img['dicbasename']}")
            img['dicfile'] = os.path.join(config['outputdir'], img['stem'], 'DIC.tif')

            # find DIC file
            logging.info(f"looking for DIC file {img['dicbasename']}")
            if os.path.isfile(os.path.join(config['inputdir'], img['dicbasename'])) is True:
                img['inputdicfile'] = os.path.join(config['inputdir'], img['dicbasename'])
                logging.info(f'..found DIC file {img["inputdicfile"]}, copying to output DIC dir for convenience')
                shutil.copy(img['inputdicfile'], img['dicfile'])
            elif os.path.isfile(os.path.join(config['inputdir'], 'DIC', img['dicbasename'])) is True:
                img['inputdicfile'] = os.path.join(config['inputdir'], 'DIC', img['dicbasename'])
                logging.info(f'..found DIC file {img["inputdicfile"]}, copying to output DIC dir for convenience')
                shutil.copy(img['inputdicfile'], img['dicfile'])
            else:
                logging.warning(f"..failed to find a DIC file, ignoring image {img['sourcefile']}.")

            img['cellmaskfile'] = os.path.join(config['outputdir'], img['stem'], 'DIC_masks.tif')
            img['nuclearmaskfile'] = os.path.join(config['outputdir'], img['stem'], 'DAPI_masks.tif')

            logging.info(f"....writing image parameters to {Path(config['outputdir']) / img['stem'] / 'img.json'}")
            with open(Path(config['outputdir']) / img['stem'] / 'img.json', 'w') as f:
                json.dump(img, f)


    # read images, convert to tiff,
    # split into DAPI and mRNA folders
    n = 0
    process = True
    for exp in config['experiments']:
        for img in exp['images']:
            if (only is None) or (img['sourcefile'] == only):

                n = n + 1
                tic = time.time()

                if process is True:

                    # check source image
                    image = BioImage(img['sourcefile'])
                    logging.info(f'opening {img["sourcefile"]} [{n}/{config['nr_images']}]')
                    logging.info(f'..found scenes: {image.scenes}')
                    for s in image.scenes:
                        if s == 'macro image':
                            logging.info(f"....ignoring scene '{s}'")
                        else:
                            image.set_scene(s)
                            logging.info(f"....reading scene '{s}' (shape {image.shape})")
                            logging.info(f"....found channels {image.channel_names}")
                            # save to tif
                            for ch in enumerate(s.replace(r'001 ', '').split(", ")):
                                if ch[1] in filter2mRNA.keys():
                                    filter = ch[1]
                                    mrna = filter2mRNA[ch[1]]
                                    img[mrna] = {}
                                    logging.info(f"......reading in channel: C={ch[0]} filter={filter} mrna={mrna}")

                                    # save layer data in tiff format
                                    layerdata = image.get_image_data("ZYX", C=ch[0])
                                    layerfile = os.path.join(config['outputdir'], img['stem'], f'{mrna}.tif')
                                    OmeTiffWriter.save(layerdata, layerfile, dim_order="ZYX", channel_names=[filter])
                                    logging.info(f'......saving {mrna} channel to {layerfile}')
                                    img[mrna]['rnafile'] = layerfile

                                    # save DAPI max projection along Z because we need it as input to compute the nuclear mask
                                    if mrna == "DAPI":
                                        maxprojdata = np.max(layerdata, axis=0)

                                        img[mrna]['maxprojfile'] = os.path.join(config['outputdir'], img['stem'], 'DAPI_max_proj.tif')
                                        OmeTiffWriter.save(maxprojdata, img[mrna]['maxprojfile'], dim_order="YX", channel_names=[filter])
                                        logging.info(f"......saving {mrna} max projection data to {img[mrna]['maxprojfile']}")

                                        logging.info(f"......expecting nuclear mask file in {img['nuclearmaskfile']}")

                                else:
                                    logging.warning(f"......ignoring unknown channel {filter}")

                img['time']['01-configure'] = time.time() - tic

    logging.info(f"output config file: {os.path.join(config['outputdir'], os.path.basename(configfile))}")
    with open(os.path.join(config['outputdir'], os.path.basename(configfile)), "w") as f:
        json.dump(config, f)

    img['time']['01-configure'] = time.time() - tic0

    logging.info(f"processed {n} image{'s' if n > 1 else ''} in {time.time() - tic0: 0.2f}s")
    logging.info("done.")


