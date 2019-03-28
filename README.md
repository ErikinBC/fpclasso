## FPC Lasso replication instructions

This repository contains the necessary code and files to reproduce the results shown in X.

In the ~/rcodes/ folder, run the shell script fpc_lasso_master.sh to generate the necessary datasets and files. The results will be save to an ~/output/ and ~/datasets/ folder. Ensure that "dir_base" and corresponding folder structure matches the shell script.

The ~/rcodes/fpc_wrapper.R contains the general FPC Lasso function that can be used.

An R Package for the FPC Lasso is under construction and will package the scripts accordingly.
