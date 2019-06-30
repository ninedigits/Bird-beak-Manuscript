# Endograft mal-apposition algorithm
@author: Max Frohlich<br>
@date:   29 June 2019
# Project overview
## Objective

Create algorithm that can quantify endograft mal-apposition from segmented image data and relate it to pre-operative geometry.

## Scope

Models are generated from CT scans using SimVascular image segmentation software. 


## Code 

There were 6 main processes to complete this project.

1. Generate models using SimVascular. This step was completed by previous members of the Stanford VIBE lab.
2. `main_inner_outer.m`: Generates the inner and outer aortic curves. Determines aortic post- and pre-operative curve lines, which include an inner point for each aortic contour. Each contour is spaced roughly 5mm apart.
3. **Fourier smoothing**: (not published) This step involves a proprietary script that uses optimized Fourier smoothing to smooth and interpolate the inner and outer curves in 0.1mm increments. 
3. `main_birdbeak_v1.m`: This script measures endograft mal-apposition and relates it to pre-operative geometry.
4. `Aorta colormap_rev2.ipynb/py`: This script generates the figures used in the paper for publication. 
5. `Aorta Outcomes Area`: This script performs the statistical testing on the pre-operative geometry and post-operative bird-beak severity.

# Abstract / Executive Summary

## Introduction

Aortic geometry has been shown to influence the development of endograft mal-apposition (bird-beaking) in thoracic endovascular aortic repair (TEVAR), but the extent of this relationship lacks clarity. The aim of this study was to develop a reproducible method of measuring bird-beak severity and pre-operatively assess the risk of occurrence.
Methods: 

## Methods 

The study included 21 patients with thoracic aortic aneurysms or type-B dissections treated with TEVAR. Computed tomography scans were used to construct 3D geometric models of pre- and post-operative aortic and endograft geometry. Post-operative bird-beaking was quantified with length, height, and angle, categorized into a bird-beak group (BBG) and no bird-beak group (NBBG) using bird-beak height ≥ 5 mm as a threshold, and correlated to pre-operative aortic geometry using aortic area, diameter, and inner curvature, and graft diameter and oversizing at the proximal landing zone.

## Results

Aortic area (976±143 vs. 834±248 mm2), diameter (35.2±2.6 vs. 32.2±4.9 mm), and inner surface curvature (0.042±0.015 vs 0.031±0.012 mm-1) were not significantly different between BBG vs. NBBG, however, the dimensionless combination metric curvature × diameter, was significantly higher in BBG (1.5±0.5 vs. 1.0±0.3, P=0.016). Inner surface curvature and curvature × diameter were significantly correlated with bird-beak height (R=0.444, P=0.044; R=0.581, P=0.006) and bird-beak angle (R=0.710, P<0.001; R=0.736, P<0.001).


## Conclusion

TEVAR bird-beak severity can be automatically and robustly quantified with geometric modeling techniques, and the combination of high pre-operative aortic diameter and inner surface curvature increases the risk of developing TEVAR bird-beaking.

