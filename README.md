# EDoF-Miniscope

This is an open source repository of the EDoF-Miniscope project in the Computational Imaging Systems Lab (CISL) at Boston University. We aim to develop next-generation head-mounted fluorescence microscopes, known as miniscopes, that push the imaging capacity in scattering media by integrating customizable and optimized diffractive optics. In this repository, we provide: <br /> 
- The hardware design of the EDoF-Miniscope in its current tabletop iteration <br />  
- The Zemax models of the EDoF-Miniscope system including the ZMX files, spectra data and coating profile <br /> 
- Matlab scripts containing the genetic algorithm used to optimize the diffractive optics, optical simulations to predict the behavior, and a simple post-processing algorithm to extract the in-focus component

## Citation

If you find this project help, please consider citing our work:

[**Joseph Greene, Yujia Xue, Jeffrey Alido, Alex Matlock, Guorong Hu, Kivilcim Kilic, Ian Davison, and Lei Tian. "EDoF-Miniscope: pupil engineering for extended depth-of-field imaging in a fluorescence miniscope" BioRxiv (2022)**](https://doi.org/10.1101/2022.08.05.502947)

## EDoF-Miniscope Overview

Extended depth of field (EDoF) microscopy has emerged as a powerful solution to greatly increase the access into neuronal populations in table-top imaging platforms. Here, we present EDoF-Miniscope, which integrates an optimized thin and lightweight binary diffractive optical element (DOE) onto the gradient refractive index (GRIN) lens of a head-mounted fluorescence miniature microscope, i.e. “miniscope”. We achieve an alignment accuracy of 70 μm to allow a 2.8X depth-of-field extension between the twin foci. We optimize the phase profile across the whole back aperture through a genetic algorithm that considers the primary GRIN lens aberrations, optical property of the submersion media, and axial intensity loss from tissue scattering in a Fourier optics forward model. Compared to other computational miniscopes, our EDoF-Miniscope produces high-contrast signals that can be recovered by a simple algorithm and can successfully capture volumetrically distributed neuronal signals without significantly compromising the speed, signal-to-noise, signal-to-background, and maintain a comparable 0.9-μm lateral spatial resolution and the size and weight of the miniature platform. We demonstrate the robustness of EDoF-Miniscope against scattering by characterizing its performance in 5-μm and 10-μm beads embedded in scattering phantoms. We demonstrate that EDoF-Miniscope facilitates deeper interrogations of neuronal populations in a 100-μm thick mouse brain sample, as well as vessels in a mouse brain. Built from off-the-shelf components augmented by a customizable DOE, we expect that this low-cost EDoF-Miniscope may find utility in a wide range of neural recording applications.

<p align="center">
  <img src="/Images/overview.png">
</p>

## How to use
### 1) Hardware design

The directory 'CAD STL' contains the CAD files of the EDoF-Miniscope. This current iteration of the EDoF-Miniscope is designed to integrate on a tabletop setup for interrogating fixed samples. All CAD models are 3D printable on lab table-top 3D printers. The files are as follow:
- "DOEscope v1.4.2 - Opt Zemax DWL 670 GRIN Opt Filter Holes.stl" contains the body of the EDoF-Miniscope presented in this work with slots for off-the-shelf optics and a custom DOE.
- "PT3A Design v2.stl" is a custom adapter that holds a manufactured DOE and affixes to a PT3A stage so it may be aligned and glued to a GRIN lens to form a DOE-GRIN component (see Supplement).
- "GRIN Holder for EDOF-Miniscopev1.stl" contains just the holder for a glued DOE-GRIN component. This allows for easy insertion of a DOE-GRIN component and the holder maybe attached to the EDoF-Miniscope body with dental apoxy.

The part list of all optical and electronic components used in the EDoF-Miniscope prototype can be found [**here**](https://docs.google.com/spreadsheets/d/1PgIITukA03SGAjqEpHsR73N81aqUN8srO4x0Fl3sK8k/edit?usp=sharing).

All information on the tabletop setup and our procedure for integrating as well as aligning the EDoF-Miniscope can be found in our supplement.

<p align="center">
  <img src="/Images/assembly.png">
</p>

### 2) Zemax model

To use the Zemax models please follow these steps:  
1. Cloning this repository  
2. Unzip "Zemax Files.zip"
3. Copy the coating file "EDoF_Miniscope_Coating_just_Chroma.DAT" and "EDoF_Miniscope_Coating.DAT" to the directory "Zemax\Coatings\". The first file contains only the custom Chroma coatings. The second file contains the custom Chroma coatings as well as calls all the existing coatings available in the "Zemax\Coatings\" directory as of Zemax version 21.3.2.
4. Copy the spectra files "gfp_emission.spcd" and "led_spectrum_interp.spcd" to the directory "Zemax\Objects\Sources\Spectrum Files\"  
5. Open "EDoF_Miniscope_Characterization.zmx" to access the model of the emission model of the EDoF-Miniscope. By having the corresponding .CFG and .ZDA in the same directory, this script will autogenerate pregenerated universal plots detailing how adjusting the components within the miniscope affects the PSF.
6. Open "sequential_64520_miniscope_aberrations.zmx" to access our characterization of the aberrations within off-the-shelf GRIN lens used in the EDoF-Miniscope (Edmund Optics, #64520)

### 3) Matlab Genetic Algorithm, Simulation and Post-Porcessing

The Matlab folder contains three folders: Genetic Algorithm, Simulation and Post-Processing. Alls scripts were developed on Matlab R2020b and function on all versions between R2020b and the current tested version (R2022b).

#### 3a) Genetic Algorithm Deployment

The DOE optimization genetic algorithm is contained in the folder "Matlab/Genetic Algorithm". Ensure that the "Helpers", "Objects", "Supp_Analysis", "Export_Fig" are within the directory. Open "Optimize_DOE_GA.m" and adjust the user parameters to fit the desired optimization conditions (simulated enviornment and optimization parameters). Run the genetic algorithm. The algorithm will generate a new folder containing figures and parameters detailing the optimized DOE. We anaylzed the stability of the optimied solution as a function of genetic algorithm parameters in the supplement.

<p align="center">
  <img src="/Images/genetic algorithm.png">
</p>

#### 3b) 
