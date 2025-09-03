import matplotlib as mpl
from pytest import Session

def pytest_sessionstart(session: Session):
    mpl.rcParams["text.usetex"] = False
