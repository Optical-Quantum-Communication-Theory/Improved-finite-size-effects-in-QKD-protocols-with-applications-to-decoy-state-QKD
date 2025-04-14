# Improved finite-size effects in QKD protocols with applications to decoy-state QKD

This is a public version of the code used in *Improved finite-size effects in QKD protocols with applications to decoy-state QKD* \[[ArXiv](https://arxiv.org/abs/2502.05382)\]. This was built for [v2.0.2](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.2) of the Open QKD Security package.

The data in Figures 2 - 6 is created by running the following files:

| Figure  | Files |
| ------------- | ------------- |
| Fig. 2  | `mainFiniteLossyQubitBB84ExpecRate.m` |
| Fig. 3  | `mainFiniteLossyQubitBB84.m` |
| Fig. 4  | Key rates using 1-norm constraints run `mainFiniteLossyQubit1Norm.m`, <br> Key rates using entrywise constraints run `mainFiniteLossyQubitBB84.m` |
| Fig. 5  | `mainFiniteDecoyBB84.m` |
| Fig. 6  | `FourSixAdaptiveMain_DifInt.m`, <br> Switch between iid and ps needs to be done in preset `FourSixWCPDifInt_AdaptivePreset.m`. Both the `symbolic` and `statistics` toolbox are required. |

The folder '4-6Protocol' contains the fixed-length version of the variable-length 4-6 protocol, presented in section IX.A.

## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### As zip
1. Download the linked version of the code from above and follow all [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/v2.0.2).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/Improved-finite-size-effects-in-QKD-protocols-with-applications-to-decoy-state-QKD
```
2. Follow all further [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/v2.0.2).
3. Also follow the additional Mosek install instructions if you want an exact match.
