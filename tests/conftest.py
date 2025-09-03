import matplotlib as mpl
from pytest import Config

def pytest_configure(config: Config):
    mpl.rcParams["text.usetex"] = False
