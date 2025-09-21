As of September 21, 2025, the independent fission yield is required to run the 
saturation simulation with custom nuclides (due to the NFY data requirement).
To get around this, simply run the saturation model, adjust the concentrations
to be 1000 for A, 100 for B, 10 for C, and 1 for D.
Then, run `mosden --main input_sat.json` followed by
`mosden --post input_sat.json` in order to use those values.
The pulse simulation can be fully run as-is with no modifications.