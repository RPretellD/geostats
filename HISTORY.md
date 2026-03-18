# geostats version history

### V0.1.0
- Initial release. 

### V1.0.0
- Solved installation issues associated to Cython. Cython files are compiled during pip installation (if C/C++ tools are available).
- Functions that were previously available in Cython only are now also available in Python.
- The library functions are loaded by default using the Python implementation, but can also be used in Cython when desired. If the Cython implementation cannot be loaded or fails for some reason, then the library defaults back to the Python implementation.
- Other code changes to accommodate the above updates.

### V1.0.1
- Added "noise" to covariance matrix to prevent singular matrix issues. 

### V1.0.2
- Modified the backend selection to facilitate going from Python (default) to Cython and back to Python. This modification also allows for the backend to be selected from [gmKriger](https://github.com/RPretellD/gmKriger).
