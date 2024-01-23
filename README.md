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
Progressive augmentation of the kOmegaSST turbulence model in OpenFOAM.
This correction enhances the performance of kOmegaSST turbulence model in 
capturing the separation phenomenon. This model has been tested on several 2D canonical flow cases
including periodic hills, curved backward-facing step, converging-diverging channel, and parametric bumps. 
The implementation of the augmented model has been developed for 2D flows but 
its applicability has been tested on 3D flows, yielding similar results.
Five coefficients can be modified by the user to change the model's performance. 
Standard optimised values are given by default in the model.
More information is available in the publication listed at the of this file.

## Target platform
The code is known to work with OpenFOAM v2312 and previous ESI versions.

## Authors
Ali Amarloo <amarloo@mpe.au.dk>
Mario Javier Rincón <mjrp@mpe.au.dk>

## Instructions

1. Download the source code using git:

         git clone https://github.com/AUfluids/KOSSTSEP.git

2 .Enter the directory where the source code has been extracted, and compile it by typing: 

         wmake

3. Add the following line to the _controlDict_ of your case:

         libs ( "libSEPIncompressibleTurbulenceModels" ) ;

4. Specify the following in _turbulentProperties_.

         RASModel kOmegaSSTSEP;
   
    (We suggest first running your case with the standard kOmegaSST and then changing the model to kOmegaSSTSEP, to avoid possible errors)
   
6. (Optional) 4 different separation equations are obtained (More info inside the published manuscript).
   The model _IV_ is the default setting but you can change it by adding the following subdictionary to _turbulentProperties_.

         separationMode  4; //optional - default:4 - off:0 | ModelI:1 | ModelII:2 | ModelIII:3 | ModelIV:4
   
    If you use 0, the separation factor is deactivated, and the standard kOmegaSST is used.
   
8. (Optional) You can also specify the 5 coefficients (C_0, C_1, C_2, separationLambda1, separationLambda2) corresponding to the separation factor in  _turbulentProperties_.
   Otherwise, these coefficients are automatically assigned with values corresponding to the model specified by _separationMode_. 

           separationLambda1   20;             //optional 
           separationLambda2   7.2513;         //optional
           C0                  -0.872209;      //optional 
           C1                  0.0131861;      //optional 
           C2                  -0.0766894;     //optional 


You can also check the test case of CBFS13700 in the folder testCases.

## Test results

For more details, refer to the publication at: LINK Progressive augmentation of RANS models for separated flow prediction by CFD-driven surrogate optimisation

Results of using KOSSTSEP for curved backward-facing step with bulk Reynolds number of 13700 **I** and **III**:
![alt text](https://github.com/AUfluids/KOSSTSEP/blob/main/testCases/CBFS13700_KOSSTSEP/contours_comparisonCBFS.png)
![alt text](https://github.com/AUfluids/KOSSTSEP/blob/main/testCases/CBFS13700_KOSSTSEP/quantitative_comparison_CBFS.png)

## How to cite
Please, cite this library using the following DOI: [10.1063/5.0174470](https://doi.org/10.1063/5.0174470).

Amarloo and Rincón (2023)

    @article{amarloo2023progressive,
      title={Progressive augmentation of turbulence models for flow separation by multi-case computational fluid dynamics driven surrogate optimization},
      author={Amarloo, Ali and Rinc{\'o}n, Mario Javier and Reclari, Martino and Abkar, Mahdi},
      journal={Physics of Fluids},
      volume={35},
      number={12},
      year={2023},
      publisher={AIP Publishing}
    }

For release-specific DOIs, click on the badge and find the DOI corresponding to the desired version in the version list.

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM® and OpenCFD® trade marks.

Detailed information on the OpenFOAM trademark can be found at

http://www.openfoam.com/legal/trademark-policy.php
http://www.openfoam.com/legal/trademark-guidelines.php
For further information on OpenCFD and OpenFOAM, please refer to

http://www.openfoam.com
