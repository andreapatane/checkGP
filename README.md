# checkGP
Matlab implementation of the methods described in [Robustness Guarantees for Bayesian Inference with Gaussian Processes](https://arxiv.org/abs/1809.06452).

## Requirements
The main requirement for running the code is of course a:
- Matlab installation (the code was tested on Matlab versions >= 2016a)

Additional toolbox may be required depending on the use:
- The optimization toolbox is used to solve quadratic programming problems. Those are used for ReLU kernel GP.
- The parallel toolbox may be used when solving quadratic programming problems, to solve the linear inequalities system in the constraints.
- The [vlfeat toolbox](http://www.vlfeat.org/install-matlab.html) is used for automatic feature extraction using SIFT. This is used in the experimnts on the MNIST dataset.

To run the experiments with the MNIST dataset it is necessary to first download the dataset. The code expect the training data to be stored into the following .csv files:
- "x_train2000.csv": a subsample of 2000 images among all the ones included in the MNIST training dataset. Data are assumed to be organised a 2000x784 matrix.
- "y_train2000.csv": labels for the images included in the file above. That is a 2000x1 matrix of integer values between 0 and 9.
- "x_test.csv": MNIST testing set (10000x784 matrix)
- "y_test.csv": MNIST testing set labels (10000x1 matrix)

## Settings
Some path and device specific parameters has to be set before running the code. Those are all included intp the m-files into utils/pc_utils/ subdirectory. Namely:
- get_data_folder(): this function has to return as output the directory in which the mnist dataset is locally stored.
- get_num_of_threads(): output of the function is the maximum number of concurrent threads for solving linear systems.
- get_result_folder(): the function output will be the directory in which results of the analyses will be stored
- setup_vlfeat(): runs the setup of the vl_feat toolbox. 