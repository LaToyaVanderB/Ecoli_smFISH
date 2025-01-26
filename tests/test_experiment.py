import pytest

from fish.experiment import Experiment

@pytest.fixture
def my_experiment():
    my_experiment = Experiment("exp24/config.json")
    return my_experiment


def test_experiment(my_experiment):
    assert len(my_experiment.channels) == 4
    my_experiment.create_image_list()
    pass


def test_experiment_from_jsons(my_experiment):
    assert len(my_experiment.channels) == 4
    my_experiment.read_image_list_from_jsons()
    pass
