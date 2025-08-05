import argparse
from mosden.preprocessing import Preprocess
from mosden.concentrations import Concentrations
from mosden.countrate import CountRate
from mosden.groupfit import Grouper
from mosden.postprocessing import PostProcess
from mosden.base import BaseClass
from . import __version__


def main():
    parser = argparse.ArgumentParser(description="MoSDeN")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"MoSDeN {__version__}")

    group.add_argument("-main", "--input", help="Input file for main run")
    group.add_argument(
        "-pre",
        "--preprocess",
        help="Input file for preprocessing")
    group.add_argument(
        "-post",
        "--postprocess",
        help="Input file for postprocessing")
    group.add_argument("-a", "--all", help="Input file for all processes")

    args = parser.parse_args()

    if args.input:
        BaseClass(args.input).clear_post_data()
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
        postprocess = PostProcess(args.postprocess)
        postprocess.run()
    elif args.all:
        BaseClass(args.all).clear_post_data()
        preprocess = Preprocess(args.all)
        preprocess.run()
        concentrations = Concentrations(args.all)
        concentrations.generate_concentrations()
        countrate = CountRate(args.all)
        countrate.calculate_count_rate()
        grouper = Grouper(args.all)
        grouper.generate_groups()
        postprocess = PostProcess(args.all)
        postprocess.run()
    else:
        print("No valid option provided. Use -h for help.")


if __name__ == "__main__":
    main()
