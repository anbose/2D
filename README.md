# Description

This is a code in C++ for simulation the diffusivity edge system as described in https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.010601 in 2-dimensional grid. This is a parallel program written using the MPI library.

## Author

Aritra Bose

### Details

We use a local Lax-Wendorf method to formulate a setup for the simulation. This setup is inspired from these two papers : https://epubs.siam.org/doi/10.1137/110853807 and https://epubs.siam.org/doi/10.1137/060673308. A version of this code, written in Julia, was used to simulate the same system in 1D for the paper https://iopscience.iop.org/article/10.1209/0295-5075/acdcb7/meta.

The directory contains the necessary files needed for the simulation. However, the file including functions (funcs.h) could not be shared due to privacy measures. Interested people can contact me.

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.



