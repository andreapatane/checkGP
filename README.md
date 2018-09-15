# checkGP
Matlab implementation for providing "Robustness guarantees for Gaussian Process".

## Requirements
The main requirement for running the code is of course:
- Matlab installation (the code was tested on Matlab versions >= 2016a)
Additional toolbox may be required depending on the use:
- The optimization toolbox is used to solve quadratic programming problems. Those are used for ReLU kernel GP.
- The parallel toolbox may be used when solving quadratic programming problems, to solve the linear inequalities system in the constraints.
- The [vlfeat toolbox](http://www.vlfeat.org/install-matlab.html) is used for automatic feature extraction using SIFT. This is used in the experimnts on the MNIST dataset.

##Settings