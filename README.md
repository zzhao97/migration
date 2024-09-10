# Source Code for Image Analysis

This is a repository associated with manuscript "Chromatin Compaction During Confined Cell Migration Induces and Reshapes Nuclear Condensates" by Zhao et al, submitted to Nature Communications, and all subsequent version of the manuscript. It includes MATLAB scripts for all analyzing fluorescence intensity, software used to perform particle imaging velocimetry [PIVlab](https://pivlab.blogspot.com/), pairwise MSD analysis for microrheology experiments, and intracellular phase diagram mapping with SVM-based classifications of phase-separated vs non-phase separated data.

## 1. System Requirements

### Software Dependencies
- **Operating Systems**:  
  - Windows 10 or higher  
  - macOS Monterey or higher  
  - Linux (tested on Ubuntu 20.04 and 22.04)

- **MATLAB Version**:  
  - MATLAB R2021b or higher

- **MATLAB Toolboxes**:
  - Image Processing Toolbox (for handling and analyzing fluorescence images)
  - Statistics and Machine Learning Toolbox (optional, for advanced statistical analysis)
  - Signal Processing Toolbox (optional, for noise reduction and filtering in fluorescence data)
  - For PIVlab, follow instruction [here](https://pivlab.blogspot.com/)

### Required Non-Standard Hardware
- For large image processing (e.g., 2048x2048 timelapse fluorescence images):
  - At least 8GB RAM is recommended
  - A multi-core CPU for optimal performance
  - A GPU is optional but will improve performance for some image processing tasks (requires Parallel Computing Toolbox)

## 2. Installation Guide

### Instructions
- Install **MATLAB** from the official [MathWorks website](https://www.mathworks.com/).

## Reference
Thielicke, W. and Stamhuis, E.J. (2014): PIVlab – Towards User-friendly, Affordable and Accurate Digital Particle Image Velocimetry in MATLAB. Journal of Open Research Software 2(1):e30, DOI: http://dx.doi.org/10.5334/jors.bl
