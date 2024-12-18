import pytest

from fish.image import Image
from fish.config import Config

@pytest.fixture
def my_config():
    cfg_file = "experiments/20240927-exp16/config.json"
    return Config(cfg_file)

@pytest.fixture
def my_image():
    vsi_file = "tests/exp16/input/MG1655_GLU_OD_0.3_left_CY5, CY3.5 NAR, CY3, DAPI_02.vsi"
    cell_file = "tests/exp16/input/MG1655_GLU_OD_0.3_left_DIC_02.tif"
    return Image(vsi_file, cell_file, my_config)

def test_init(my_image, my_config):
    assert my_config.cfg_file == "experiments/20240927-exp16/config.json"
    assert my_image.vsi_file == "tests/exp16/input/MG1655_GLU_OD_0.3_left_CY5, CY3.5 NAR, CY3, DAPI_02.vsi"
    assert my_image.cell_file == "tests/exp16/input/MG1655_GLU_OD_0.3_left_DIC_02.tif"

