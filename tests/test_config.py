import pytest

from fish.config import Config

@pytest.fixture
def my_config():
    my_config = Config("exp16/config.json")
    return my_config

def test_config(my_config):
    assert len(my_config.channels) == 4
