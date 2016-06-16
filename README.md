# DESDM Standalone

Tools for running the DESDM single-epoch pipeline in standalone mode
at Fermilab. Supports both local and grid computing.


## Code Structure

There are several levels to the code:
* desdmlib.py - The base library for application control. Essentially builds python object wrappers for the individual command-line components of the DESDM pipeline (i.e., crosstalk, pixcorrect, immask, etc.). Parses the config file and builds the command-line arguments for the individual components.
* run_desdm.py - Executable script for running the DESDM pipeline for all CCDs in a single exposure. Essentially a loop through application wrappers in desdmlib.
* job_createFiles - Higher level control tool for setting up grid job submission. Creates shell scripts wrapping the environment setup and configuration for run_desdm. Takes as input a list of exposures for processing.

