# Deformable Mirror FEM Analysis for Zernike Mode Reproduction (atmospheric correction)

<p align="center">
  <img src="figures/FEM_data/trefoil.jpg" width="300">
</p>

This repository contains the numerical framework developed to evaluate the mechanical performance of a **deformable** (adaptive) **mirror** (DM) in reproducing optical **Zernike modes** using finite element analysis (**FEM**) implemented in **FreeCAD**. The code implements the forward and inverse modeling of a pin-actuated deformable mirror and provides quantitative metrics of surface fitting accuracy.

The framework is based on a parametric, physics-driven approach and was developed as part of the numerical work supporting a peer-reviewed [study](https://opg.optica.org/abstract.cfm?URI=AOPT-2024-OTh1F.2) on deformable mirror design for adaptive optics applications.

---

## Overview

The deformable mirror is modeled as a thin circular plate actuated by an array of plungers arranged in configurable polar layouts. For a given actuator configuration and mirror geometry, the code:

-Generates a fully parametric 3D model of the mirror and actuators in **FreeCAD**

-Sets up and solves a mechanical FEM problem using CalculiX

-Computes **Actuator Influence Functions** (AIFs) via unit force excitation

-Constructs the AIF matrix and computes its Moore–Penrose pseudo-inverse

-Projects target Zernike modes **(Noll indexing)** onto the actuator basis

-Applies the resulting actuator forces and evaluates the mirror deformation

-Quantifies performance using **RMS** and peak residual surface errors

-Exports nodal displacement data for post-processing and visualization

This implementation focuses on the mechanical response and Zernike mode fitting accuracy of the deformable mirror. Optimization routines used in the associated publication are not included in this repository.

---
## Requirements

FreeCAD 0.21.2 (please note that some features may differ in more recent versions)

CalculiX (via FreeCAD FEM tools)

Python 

All simulations are executed within the Macro FreeCAD (Python) environment.

---

## Running the FEM Analysis in FreeCAD


1.Place the next files in your working directory: 

FEM_main.FCMacro

fem_analysis.py

mirror_geometry.py

zernike_utils.py


2. Running the Macro

-Launch FreeCAD

-Open the Macro dialog:

-Macro → Macros…

-Select FEM_main.FCMacro

-Click Execute
