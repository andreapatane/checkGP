# checkGP
Matlab implementation of the methods described in [Robustness Guarantees for Bayesian Inference with Gaussian Processes](https://arxiv.org/abs/1809.06452). Here we provide custom code for GP training in Matlab, and all the scripts and functions used to obtain the plots included in the paper.

## Requirements
The main requirement for running the code is of course a:
- Matlab installation (the code was tested on Matlab versions >= 2016a)

Additional toolbox may be required depending on the use:
- The optimization toolbox is used to solve quadratic programming problems. Those are used for ReLU kernel GP.
- The parallel toolbox may be used when solving quadratic programming problems, to solve the linear inequalities system in the constraints.
- The [vlfeat toolbox](http://www.vlfeat.org/install-matlab.html) is used for automatic feature extraction using SIFT. This is used in the experiments on the MNIST dataset.
- The [Random Field Simulation toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/27613-random-field-simulation) is used for efficient sampling from Gaussian Process. This is used to obtain empirical safety and invariance estimation. In the paper this is used only for comparison with the formal method proposed. 

Additionally, to run the experiments with the MNIST dataset it is necessary to first download the dataset. The code expect the training data to be stored into .csv files of the following form:
- "x_train.csv": [or similar name] a subsample of N images among all the ones included in the MNIST training dataset. Data are assumed to be organised as a Nx784 matrix.
- "y_train.csv": [or similar name] labels for the images included in the file above. That is a Nx1 matrix of integer values between 0 and 9.
- "x_test.csv": [or similar name] MNIST testing set (10000x784 matrix)
- "y_test.csv": [or similar name] MNIST testing set labels (10000x1 matrix)

## Settings
Some path and device specific parameters has to be set before running the code. Those are all directly included into the experiments scripts. Namely the variables to be updated are:
- gp_training_opts.data_folder : path to the directory in which the MNIST dataset is locally stored (used only for MNIST analysis).
- result_folder: path to the directory in which results of the analyses will be stored.
- var_comp_opts.vl_setup_script_path: path to the vl_feat toolbox. 