import argparse
from mosden.preprocessing import Preprocess
from mosden.concentrations import Concentrations
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
        concentrations = Concentrations(args.input)
        concentrations.generate_concentrations()
        raise NotImplementedError("Need to implement delayed neutron count rates and grouping")
    elif args.preprocess:
        preprocess = Preprocess(args.preprocess)
        preprocess.run()
    elif args.postprocess:
        print(args.postprocess)
    else:
        print("No valid option provided. Use -h for help.")

if __name__ == "__main__":
    main()
