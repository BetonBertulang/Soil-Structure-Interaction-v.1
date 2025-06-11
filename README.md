# 2D Soil-Structure Interaction Analysis with OpenSeesPy

![Model Visualization](model_visualization.png) *Placeholder for model visualization*

## Overview
This repository contains an OpenSeesPy script (`OOP.py`) for performing 2D soil-structure interaction (SSI) analysis. The script models the dynamic response of a soil domain with a surface-founded structure and piles subjected to seismic loading, using an object-oriented approach for easy parameter management.

## Key Features
- **Automatic Mesh Generation**: Optimal soil mesh size determined by wave propagation requirements
- **Advanced Material Modeling**: 
  - `PressureIndependMultiYield` for nonlinear cyclic soil behavior
- **Boundary Conditions**:
  - Periodic boundaries (equalDOF) on sides
  - Viscous dashpot base for wave absorption
- **Structural Components**:
  - SDOF structure with pile foundation
  - Elastic beam-column elements
- **Two-Phase Analysis**:
  1. Gravity analysis (initial stresses)
  2. Dynamic analysis (seismic response)

## Quick Start

### Installation
```bash
pip install openseespy numpy matplotlib opsvis
