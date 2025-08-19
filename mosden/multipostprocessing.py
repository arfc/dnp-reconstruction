from mosden.base import BaseClass

class MultiPostProcess(BaseClass):
    def __init__(self, input_paths: list[str]) -> None:
        print(input_paths)
        return None