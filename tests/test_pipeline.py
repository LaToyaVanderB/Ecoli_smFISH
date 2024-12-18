import pytest

from fish.image import Image
from fish.config import Config

@pytest.fixture
def my_config():
    cfg_file = "exp16/config.json"
    return Config(cfg_file)


@pytest.fixture
def my_image(my_config):
    vsi_file = "MG1655_GLU_OD_0.3_left_CY5, CY3.5 NAR, CY3, DAPI_02.vsi"
    cell_file = "MG1655_GLU_OD_0.3_left_DIC_02.tif"
    return Image(my_config, vsi_file, cell_file)

def test_image(my_image):
    assert my_image.metadata.vsi_file == "MG1655_GLU_OD_0.3_left_CY5, CY3.5 NAR, CY3, DAPI_02.vsi"
    assert my_image.metadata.cell_file == "MG1655_GLU_OD_0.3_left_DIC_02.tif"
    assert my_image.metadata.outputdir == "exp16/output_oo"
    assert my_image.metadata.savepath == "exp16/output_oo/MG1655_GLU_OD_0.3_left_02"

    # load image (~ 01-configure)
    my_image.read_image()
    my_image.read_cells()
    my_image.align()
    my_image.create_grgb()

    # save image (write to dir)
    my_image.save_layers()
    my_image.save_metadata()

    # segment image (~ 02-segment)
    my_image.segment()
    pass

    # detect spots (~ 03-detect-spots)

    # decompose spots (~ 04-decompose-spots)

    # assign spots (05-assign-spots)



def test_experiment(my_experiment, my_config):
    # iterate on all images in an experiment

    assert True