import argparse
from . import runme
from . import preprocessor
from . import postprocessor  # you'll need to create this
from . import __version__

def main():
    parser = argparse.ArgumentParser(description="MoSDeN CLI")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-v", "--version", action="version", version=f"MoSDeN {__version__}")

    group.add_argument("-i", "--input", help="Input file for main run")
    group.add_argument("-p", "--preprocess", help="Input file for preprocessing")
    group.add_argument("-pp", "--postprocess", help="Input file for postprocessing")

    args = parser.parse_args()

    if args.input:
        runme.run(args.input)
    elif args.preprocess:
        preprocessor.preprocess(args.preprocess)
    elif args.postprocess:
        postprocessor.postprocess(args.postprocess)

if __name__ == "__main__":
    main()
