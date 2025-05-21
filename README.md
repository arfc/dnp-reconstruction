# MoSDeN
Molten Salt Delayed Neutron (MoSDeN) is a tool used for reconstruction of delayed neutron precursor groups in molten salt reactors.


## Using the tool

### Preprocessing
Preprocessing should be run to generate any data needed (if it does not already exist as processed data).
This data is dependent on the energy of the irradiating neutrons as well as the fissile nuclide target.

### Processing
Processing consists of three steps:
1. Generate concentrations (or collect fission yield data).
2. Generate the delayed neutron count rate.
3. Fit a set of delayed neutron precursor group parameters that best fit the count rate.

### Postprocessing
Postprocessing handles plotting and data analysis from the processed results, including analysis of each step.