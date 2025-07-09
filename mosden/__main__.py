import argparse
from . import concentrations
from . import input_parser
from . import preprocessing
from . import postprocessing
from . import __version__

def main():
    parser = argparse.ArgumentParser(description="MoSDeN CLI")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-v", "--version", action="version", version=f"MoSDeN {__version__}")

    group.add_argument("-i", "--input", help="Input file for main run")
    group.add_argument("-pre", "--preprocess", help="Input file for preprocessing")
    group.add_argument("-post", "--postprocess", help="Input file for postprocessing")

    args = parser.parse_args()

    if args.input:
        print(args.input)
    elif args.preprocess:
        print(args.preprocess)
    elif args.postprocess:
        print(args.postprocess)
    else:
        print("No valid option provided. Use -h for help.")

if __name__ == "__main__":
    main()
