# specfem_beachball
Tooltip to plot moment tensor and body wave piercing points from specfem3d files. It requires the pygmt package.

The scrip takes specfem3d's CMTSOLUTION and STATIONS files as input to generate a full moment tensor plot.
The generated plot uses the CMTSOLUTION file to extract the full moment tensor parameters, used as input in pygmt's `meca` function. The STATIONS file information are used to compute the body wave piercing point location on the lower half sphere of the moment tensor.
