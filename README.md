# MoSDeN
Molten Salt Delayed Neutron (MoSDeN) is a tool used for reconstruction of delayed neutron precursor groups in molten salt reactors.

## History
This tool had a previous version in this repository accessible with git hash `b56528a4`.
https://doi.org/10.5281/zenodo.14888551

## Understanding the tool

### Preprocessing
Preprocessing should be run to generate any data needed (if it does not already exist as processed data).
This data is dependent on the energy of the irradiating neutrons as well as the fissile nuclide target.

#### Getting data
The exact organization of raw, unprocessed data is flexible, with some notable exceptions:
- OpenMC chain files to be should all be in a subdirectory labeled with "omcchain" (see `preprocessing.py` for all keywords)
- ENDF NFY data should all be in a subdirectory labeled `nfy`
- ENDF NFY files should be named "nfy-<ZZZ>_<ID>_<AAA>.csv", so 235U would be `nfy-092_U_235.csv`.
- IAEA beta-delayed neutron emission data should be in a directory `iaea` and be called `eval.csv` (default when downloading data).

Data can be collected from different sources:
- [OpenMC depletion chains](https://openmc.org/depletion-chains/): these give half-lives and independent fission yields (linearly interpolated energy dependence)
- [ENDF data](https://www.nndc.bnl.gov/endf-releases/): these give (currently) cumulative fission yields (with energy dependence based on nearest energy)
- [IAEA data](https://www-nds.iaea.org/beta-delayed-neutron/database.html): these give emission probabilities and half-lives

### Processing
Processing consists of three steps:
1. Generate concentrations (or collect fission yield data).
2. Generate the delayed neutron count rate.
3. Fit a set of delayed neutron precursor group parameters that best fit the count rate.

### Postprocessing
Postprocessing handles plotting and data analysis from the processed results, including analysis of each step.

## Using the tool

`pip install .` or `pip install -e .` can be used to make the package available to use on the command line as `mosden` (or in an editable version for development purposes).
Use `mosden -i <input.json>` to do a normal run, `mosden -pre <input.json>` for preprocessing, or `mosden -post <input.json>` for post-processing.
Alternatively, `mosden --all <input.json>` can be used to run preprocessing, calculations, and post-processing.

### Input file

#### Log level
- [<10] is the debug level
- [<20] is the info level (This is the suggested level)
- [<30] is the warning level
- [<40] is the error level
- [<50] is the critical level