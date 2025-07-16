# MoSDeN
Molten Salt Delayed Neutron (MoSDeN) is a tool used for reconstruction of delayed neutron precursor groups in molten salt reactors.


## Understanding the tool

### Preprocessing
Preprocessing should be run to generate any data needed (if it does not already exist as processed data).
This data is dependent on the energy of the irradiating neutrons as well as the fissile nuclide target.

#### Getting data
The exact organization of raw, unprocessed data is flexible, with some notable exceptions:
- OpenMC chain files to be should all be in a subdirectory labeled "omcchain"
- OpenMC chain files should be named "chain_<data>_<energy>.csv", where `data` is "endfb71" or similar, and `energy` is "pwr" for thermal spectrum or "sfr" for fast spectrum neutrons.
- ENDF NFY data should all be in a subdirectory labeled "nfy"
- ENDF NFY files should be named "nfy-<ZZZ>_<ID>_<AAA>.csv", so 235U would be `nfy-092_U_235.csv`.
- IAEA beta-delayed neutron emission data should be called `eval.csv` (default when downloading data).

Data can be collected from different sources:
- [OpenMC depletion chains](https://openmc.org/depletion-chains/): these give half-lives and independent fission yields (linearly interpolated energy dependence)
- [ENDF data](https://www.nndc.bnl.gov/endf-releases/): these give (currently) cumulative fission yields (with energy dependence based on nearest energy)
- [IAEA data](https://www-nds.iaea.org/beta-delayed-neutron/database.html): these give emission probabilities and half-lives

### Processing
Processing consists of three steps:
1. Generate concentrations (or collect fission yield data).
2. Generate the delayed neutron count rate.
3. Fit a set of delayed neutron precursor group parameters that best fit the count rate.

The name of the processed data files is configured in the `mosden/utils/input_handler`, where different data `names` are aligned to specific file names.
This is also one of the locations where the tool will have to be modified to allow new datasets.

### Postprocessing
Postprocessing handles plotting and data analysis from the processed results, including analysis of each step.

## Using the tool

`pip install .` or `pip install -e .` can be used to make the package available to use on the command line as `mosden` (or in an editable version for development purposes).
Use `mosden -i <input.json>` to do a normal run, `mosden -pre <input.json>` for preprocessing, or `mosden -post <input.json>` for post-processing.
Alternatively, `mosden --all <input.json>` can be used to run preprocessing, calculations, and post-processing.