# geostats: A suite of geostatistical tools

This repository provides a suite of geostatistical tools, including functions to compute the ground-motion correlation structure according to the relations by Bodenmann et al. (2023), and to do Kriging interpolation.


## Examples
- [Example 1 (Python backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example1_Python.ipynb): Ordinary Kriging using the $\rho_{E}$ correlation model. 
- [Example 1 (Cython backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example1_Cython.ipynb): Ordinary Kriging using the $\rho_{E}$ correlation model. 
- [Example 2 (Python backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example2_Python.ipynb): Ordinary Kriging using the $\rho_{EA}$ correlation model. 
- [Example 2 (Cython backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example2_Cython.ipynb): Ordinary Kriging using the $\rho_{EA}$ correlation model. 
- [Example 3 (Python backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example3_Python.ipynb): Ordinary Kriging using a user-defined correlation model. 
- [Example 3 (Cython backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example3_Cython.ipynb): Ordinary Kriging using a user-defined correlation model. 
- [Example 4 (Python backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example4_Python.ipynb): Krige a map using a user-defined correlation model. 
- [Example 4 (Cython backend)](https://github.com/RPretellD/geostats/blob/main/examples/Example4_Cython.ipynb): Krige a map using a user-defined correlation model. 


## Acknowledgements
- Implementation of the Kriging code benefited from Scott Brandenberg's [random field](https://github.com/sjbrandenberg/ucla_geotech_tools/tree/main/src/ucla_geotech_tools) python package.
- Some of the cython functions to compute ground-motion correlation are based on Lukas Bodenmann's [python functions](https://github.com/bodlukas/ground-motion-correlation-bayes).


## Citation
If you use these codes, please cite:<br>
Pretell, R. (2026). geostats: A suite of geostatistical tools (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.10253690 <br>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10253690.svg)](https://doi.org/10.5281/zenodo.10253690)


## Contact
For any questions or comments, contact Renmin Pretell (rpretell at unr dot edu).