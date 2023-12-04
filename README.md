# Geostats

This repository provides a suite of geostatistical tools, including functions to compute the ground-motion correlation structure according to the relations by Bodenmann et al. (2023), and to do Kriging interpolation.


### Example
- [Example 1](https://github.com/RPretellD/geostats/blob/main/Examples/Example_1.ipynb): Ordinary Kriging using the $\rho_{E}$ correlation model by Bodenmann et al. (2023)
- [Example 2](https://github.com/RPretellD/geostats/blob/main/Examples/Example_2.ipynb): Ordinary Kriging using the $\rho_{EA}$ correlation model by Bodenmann et al. (2023)
- [Example 3](https://github.com/RPretellD/geostats/blob/main/Examples/Example_3.ipynb): Ordinary Kriging using any preferred correlation structure
- [Example 4](https://github.com/RPretellD/geostats/blob/main/Examples/Example_4.ipynb): Krige a map


### Acknowledgements
- Implementation of the Kriging code benefited from Scott Brandenberg's [random field](https://github.com/sjbrandenberg/ucla_geotech_tools/tree/main/random_field) python package.
- Some of the cython functions to compute ground-motion correlation are based on Lukas Bodenmann's [python functions](https://github.com/bodlukas/ground-motion-correlation-bayes).


### Citation
If you use these codes, please cite:<br>
Renmin Pretell. (2023). RPretellD/geostats: Initial release (0.1.0). Zenodo. https://doi.org/10.5281/zenodo.10253691 <br>

[![DOI](https://zenodo.org/badge/716446689.svg)](https://zenodo.org/doi/10.5281/zenodo.10253690)