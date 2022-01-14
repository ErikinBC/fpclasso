## FPC Lasso replication instructions

This repository contains the necessary code and files to reproduce all results from [The False Positive Control Lasso](https://arxiv.org/abs/1903.12584) paper.

In the rcodes/ folder, run the shell script `fpc_lasso_master.sh` to generate the necessary datasets and files. The results will be save to the output/ and datasets/ folders. Ensure that `dir_base` and corresponding folder structure matches the shell script.

The `rcodes/fpc_wrapper.R` contains the general FPC Lasso function that can be used for other applications.
