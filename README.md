# MoSDeN
Molten Salt Delayed Neutron (MoSDeN) is a tool used for reconstruction of delayed neutron precursor groups in molten salt reactors.


## Understanding the tool

### Preprocessing
Preprocessing should be run to generate any data needed (if it does not already exist as processed data).
This data is dependent on the energy of the irradiating neutrons as well as the fissile nuclide target.

#### Getting data
The exact organization of data is flexible, with some notable exceptions:
- OpenMC chain files to be preprocessed should all be in a subdirectory labeled "omcchain"

Data can be collected from different sources:
- [OpenMC depletion chains](https://openmc.org/depletion-chains/)
- [ENDF data](https://www.nndc.bnl.gov/endf-releases/)

### Processing
Processing consists of three steps:
1. Generate concentrations (or collect fission yield data).
2. Generate the delayed neutron count rate.
3. Fit a set of delayed neutron precursor group parameters that best fit the count rate.

### Postprocessing
Postprocessing handles plotting and data analysis from the processed results, including analysis of each step.

## Using the tool

`pip install .` or `pip install -e .` can be used to make the package available to use on the command line as `mosden` (or in an editable version for development purposes).
Use `mosden -i <input.json>` to do a normal run, `mosden -pre <preproc.json>` for preprocessing, or `mosden -post <postproc.json>` for post-processing.