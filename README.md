# EDoF-Miniscope

This is an open source repository of the EDoF-Miniscope project in the Computational Imaging Systems Lab (CISL) at Boston University. We aim to develop next-generation head-mounted fluorescence microscope, known as miniscopes, that push the imaging capacity in scattering media by integrating customizable and optimized diffractive optics. In this repository, we provide: 1) the hardware design of the EDoF-Miniscope in its current tabletop iteration; 2) the Zemax models of the EDoF-Miniscope system including the ZMX files, spectra data and coating profile; 3) Matlab scripts containing the genetic algorithm used to optimize the diffractive optics, optical simulations to predict the behavior, and a simple post-processing algorithm to extract the in-focus component

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

The directory 'CAD_models' contains the CAD files of the EDoF-Miniscope. This current iteration of the EDoF-Miniscope is designed to integrate on a tabletop setup for interrogating fixed samples. All CAD models are 3D printable on lab table-top 3D printers. The subdirectory 'assembly' further provides further information and files on the EDoF-Minisope.

The part list of all optical and electronic components used in the EDoF-Miniscope prototypes can be found [**here**](https://docs.google.com/spreadsheets/d/1PgIITukA03SGAjqEpHsR73N81aqUN8srO4x0Fl3sK8k/edit?usp=sharing).

<p align="center">
  <img src="/Images/assembly.png">
</p>

### 2) Zemax model

To use the Zemax models please follow these steps:  
1.) cloning this repository  
2.) copy the coating file "cm2_coating_profiles_ver2.DAT" to the directory "Zemax\Coatings\"  
3.) copy the spectra files "gfp_emission.spcd" and "led_spectrum_interp.spcd" to the directory "Zemax\Objects\Sources\Spectrum Files\"  
4.) copy the CAD files (end with '.stl') to the directory "Zemax\Objects\CAD Files\"  
5.) then open "CM2_V1_opensource.zos" or "CM2_V2_opensource.zos" in Zemax to view the CM<sup>2</sup> design and ray tracing results.  
6.) Pre-rendered ray tracing data can be downloaded [**here**]     
    (https://drive.google.com/drive/folders/10xuGhUDethVntPEySfjDurRVyuz1iKvr?usp=sharing).  

### 3) Matlab Genetic Algorithm, Simulation and Post-Porcessing

The Matlab folder contains three folders: Genetic Algorithm, Simulation and Post-Processing. Alls scripts were developed on Matlab R2020b and function on all versions between R2020b and the current tested version (R2022b). 

The script "cm2_related_code.m" in the "Algorithm" folder provides a demo of CM<sup>2</sup> 3D reconstruction pipeline on a simulated measurement using [down-sampled PSFs](https://drive.google.com/drive/folders/10xuGhUDethVntPEySfjDurRVyuz1iKvr?usp=sharing). A full-scale experimental measurement is also provided under the "Algorithm" direcory, which requires large system memory to run. The GIF file below shows the flying-through of a reconstructed 3D object (a fluorescent fiber sample) from an experimental measurement.

<p align="center">
  <img src="/Images/genetic algorithm.png">
</p>

