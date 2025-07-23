import argparse
from mosden.preprocessing import Preprocess
from mosden.concentrations import Concentrations
from mosden.countrate import CountRate
from mosden.groupfit import Grouper
from . import __version__

def main():
    parser = argparse.ArgumentParser(description="MoSDeN CLI")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-v", "--version", action="version", version=f"MoSDeN {__version__}")

    group.add_argument("-main", "--input", help="Input file for main run")
    group.add_argument("-pre", "--preprocess", help="Input file for preprocessing")
    group.add_argument("-post", "--postprocess", help="Input file for postprocessing")
    group.add_argument("-a", "--all", help="Input file for all processes")

    args = parser.parse_args()

    if args.input:
        concentrations = Concentrations(args.input)
        concentrations.generate_concentrations()
        countrate = CountRate(args.input)
        countrate.calculate_count_rate()
        grouper = Grouper(args.input)
        grouper.generate_groups()
    elif args.preprocess:
        preprocess = Preprocess(args.preprocess)
        preprocess.run()
    elif args.postprocess:
        raise NotImplementedError("Postprocessing is not yet implemented")
    elif args.all:
        preprocess = Preprocess(args.all)
        preprocess.run()
        concentrations = Concentrations(args.all)
        concentrations.generate_concentrations()
        countrate = CountRate(args.all)
        countrate.calculate_count_rate()
        grouper = Grouper(args.all)
        grouper.generate_groups()
        print('Need postprocessing')
    else:
        print("No valid option provided. Use -h for help.")

if __name__ == "__main__":
    main()
