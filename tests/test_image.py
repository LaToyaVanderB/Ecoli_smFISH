import pytest

from fish.experiment import Experiment
from fish.image import Image

import time


@pytest.fixture
def my_experiment():
    cfg_file = "exp16/config.json"
    return Experiment(cfg_file)


@pytest.fixture
def my_image(my_experiment):
    my_params = {
        'vsi_file': "MG1655_GLU_OD_0.3_left_CY5, CY3.5 NAR, CY3, DAPI_02.vsi",
        'cell_file': "MG1655_GLU_OD_0.3_left_DIC_02.tif"
    }
    return Image.from_dict(my_params, my_experiment)


@pytest.fixture
def my_image_from_json(my_experiment):
    img_json = "exp16/output_oo/MG1655_GLU_OD_0.3_left_02/img.json"
    return Image.from_json(img_json, my_experiment)


def test_from_json(my_image_from_json):
    assert isinstance(my_image_from_json, Image) == True
    assert isinstance(my_image_from_json.experiment, Experiment) == True
    pass


def test_from_params(my_image):
    assert isinstance(my_image, Image) == True
    assert isinstance(my_image.experiment, Experiment) == True
    pass


def test_pipeline(my_image):
    # load image (~ 01-configure)
    tic = time.time()
    my_image.read_image()
    my_image.read_cells()
    my_image.align()
    my_image.create_grgb()

    # crop image
    my_image.crop()

    # save image (write to dir)
    my_image.save_layers()
    my_image.time['01-configure'] = time.time() - tic
    my_image.save_metadata()

    # segment image (~ 02-segment)
    tic = time.time()
    my_image.segment_cells()
    my_image.time['02-segment-cells'] = time.time() - tic
    my_image.save_metadata()

    tic = time.time()
    my_image.segment_dapi()
    my_image.time['02-segment-dapi'] = time.time() - tic
    my_image.save_metadata()

    # postprocess masks
    tic = time.time()
    my_image.postprocess_masks()
    my_image.time['02-segment-pp'] = time.time() - tic
    my_image.save_metadata()

    # cropping: better process the whole picture and then crop
    # (i.e. select cells from good area)
    # duh...

    # detect spots (~ 03-detect-spots)
    tic = time.time()
    my_image.find_focus()
    my_image.filter()
    my_image.detect_spots()
    my_image.time['03-detect-spots'] = time.time() - tic
    my_image.save_metadata()

    # decompose spots (~ 04-decompose-spots)
    tic = time.time()
    my_image.decompose_spots()
    my_image.time['04-decompose-spots'] = time.time() - tic
    my_image.save_metadata()


    # assign spots (05-assign-spots)
    my_image.assign_spots()

    # save image (json pickle)
    tic = time.time()
    my_image.save("05")
    my_image.time['05-save'] = time.time() - tic
    my_image.save_metadata()

    pass
