# Deformable Mirror FEM Analysis for Zernike Mode Reproduction (Atmospheric correction)

This repository contains the numerical framework developed to evaluate the mechanical performance of a deformable (adaptive) mirror (DM) in reproducing optical Zernike modes using finite element analysis (FEA). The code implements the forward and inverse modeling of a pin-actuated deformable mirror and provides quantitative metrics of surface fitting accuracy.

The framework is based on a parametric, physics-driven approach and was developed as part of the numerical work supporting a peer-reviewed [study](https://opg.optica.org/abstract.cfm?URI=AOPT-2024-OTh1F.2) on deformable mirror design for adaptive optics applications.

---

## Overview

The deformable mirror is modeled as a thin circular plate actuated by an array of plungers arranged in configurable polar layouts. For a given actuator configuration and mirror geometry, the code:

-Generates a fully parametric 3D model of the mirror and actuators in FreeCAD
-Sets up and solves a mechanical FEM problem using CalculiX
-Computes Actuator Influence Functions (AIFs) via unit force excitation
-Constructs the AIF matrix and computes its Moore–Penrose pseudo-inverse
-Projects target Zernike modes (Noll indexing) onto the actuator basis
-Applies the resulting actuator forces and evaluates the mirror deformation
-Quantifies performance using RMS and peak residual surface errors
-Exports nodal displacement data for post-processing and visualization

This implementation focuses on the mechanical response and Zernike mode fitting accuracy of the deformable mirror. Optimization routines used in the associated publication are not included in this repository.

---
