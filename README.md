# checkGP
Matlab implementation for providing "Robustness guarantees for Gaussian Process".

## Requirements
The main requirement for running the code is of course a:
- Matlab installation (the code was tested on Matlab versions >= 2016a)

Additional toolbox may be required depending on the use:
- The optimization toolbox is used to solve quadratic programming problems. Those are used for ReLU kernel GP.
- The parallel toolbox may be used when solving quadratic programming problems, to solve the linear inequalities system in the constraints.
- The [vlfeat toolbox](http://www.vlfeat.org/install-matlab.html) is used for automatic feature extraction using SIFT. This is used in the experimnts on the MNIST dataset.

To run the experiments with the MNIST dataset it is necessary to first download the dataset. The code expect the training data to be stored into the following .csv files:
- x_train2000.csv


## Settings
To run