# KOSSTSEP
Progressive augmentation of $k-\omega$ SST turbulence model for separated flows
as proposed by Amarloo and Rincón (2023) for OpenFOAM.
Developed by Fluid Physics & Turbulence research group at Aarhus University.

Copyright Information
    Copyright © 2004-2018 OpenCFD Ltd (ESI Group)

## License
    This program is free software: you can redistribute and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Description
Progressive augmentation of the $k-\omega$ SST turbulence model in OpenFoam.
This correction will enhance the performance of $k-\omega$ SST turbulence model in 
capturing the separation phenomenon. This model has been tested on several cases
including periodic hills, curved backward-facing step, converging-diverging channel, and parametric bumps. 
Implementation of the Explicit algebraic Reynolds stress correction model (EARSCM)
The model has been developed for 2D flows but its applicability has been tested on 3D
flows, yielding similar results.
Five coefficients can be modified by the user to change the model's performance. 
Standard optimised values are given by default in the model.
More info is available in the publication listed at the end.

## Target platform
The code is known to work with OpenFOAM v2212.

## Authors
Ali Amarloo <amarloo@mpe.au.dk>
Mario Javier Rincón <mjrp@mpe.au.dk>

## Instructions

1. Download the source code using git:

         git clone git://github.com/AUfluids/KOSSTSEP.git

2 .Enter the directory where the source code has been extracted, and compile it by typing: 

         wmake

3. Add the following line to the _controlDict_ of your case:

         libs ( "libSEPIncompressibleTurbulenceModels" ) ;

4. Specify the following in _turbulentProperties_.

         RASModel kOmegaSSTSEP;

5. (Optional) Add the following subdictionary to _turbulentProperties_. Otherwise, the modelIV:4 is the default model to use. 

         separationMode  4; \\optional - default:4 - off:0 | ModelI:1 | ModelII:2 | ModelIII:3 | ModelIV:4

6. (Optional) You can also specify the 5 coefficients ($C_0, C_1, C_2, \lambda_1, \lambda_2$) corresponding to the separation factor in  _turbulentProperties_.
   Otherwise, the modelIV is the default model to use. 

           separationLambda1   20;             \\optional - default taken from separationMode 4
           separationLambda2   7.2513;         \\optional - default taken from separationMode 4
           C0                  -0.872209;      \\optional - default taken from separationMode 4
           C1                  0.0131861;      \\optional - default taken from separationMode 4
           C2                  -0.0766894;     \\optional - default taken from separationMode 4


You can also check the test case of CBFS13700 in the folder testCases.

## Test results

For more details, refer to the publication at: LINK Progressive augmentation of RANS models for separated flow prediction by CFD-driven surrogate optimisation

Results of using KOSSTSEP for curved backward-facing step with bulk Reynolds number of 13700 **I** and **III**:
![alt text](https://github.com/AUfluids/KOSSTSEP/blob/main/testCases/CBFS13700_KOSSTSEP/contours_comparisonCBFS.jpg)
![alt text](https://github.com/AUfluids/KOSSTSEP/blob/main/testCases/CBFS13700_KOSSTSEP/quantitative_comparison_CBFS.jpg)

## How to cite
Please, cite this library using the following DOI: DOI.

Amarloo and Rincón (2023)
@article{amarlooRincon2023progressive,
  title={Progressive augmentation of RANS models for separated flow prediction by CFD-driven surrogate optimisation},
  author={Amarloo, Ali and Rinc{\'o}n, Mario Javier and Abkar, Mahdi},
  journal={Internaltional Journal of Heat and Fluid Flow},
  volume={000},
  pages={000},
  year={0000},
  publisher={Elsevier}
}
(under review)

For release-specific DOIs, click on the badge and find the DOI corresponding to the desired version in the version list.

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM® and OpenCFD® trade marks.

Detailed information on the OpenFOAM trademark can be found at

http://www.openfoam.com/legal/trademark-policy.php
http://www.openfoam.com/legal/trademark-guidelines.php
For further information on OpenCFD and OpenFOAM, please refer to

http://www.openfoam.com
