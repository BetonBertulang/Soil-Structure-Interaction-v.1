2D Soil-Structure Interaction Analysis with OpenSeesPy
This repository contains an OpenSeesPy script (OOP.py) for performing a 2D soil-structure interaction (SSI) analysis. The script models the dynamic response of a soil domain with a surface-founded structure and piles subjected to seismic loading.

The code is structured using an object-oriented approach to make parameters easy to manage and modify.


(Image is a placeholder for model visualization)

Features
Automatic Mesh Generation: The soil mesh size is automatically determined based on wave propagation requirements, ensuring the accurate resolution of seismic waves up to a specified maximum frequency.

Comprehensive Material Modeling: It uses the PressureIndependMultiYield material model to capture nonlinear cyclic soil behavior.

Advanced Boundary Conditions:

Periodic Boundaries: Simulates an infinite soil domain using equalDOF constraints on the lateral boundaries.

Viscous Dashpot: A non-reflective boundary at the base absorbs downward-propagating waves.

Structural Modeling: A single-degree-of-freedom (SDOF) structure with a pile foundation is modeled using elastic beam-column elements.

Two-Phase Analysis:

Gravity Analysis: Establishes the initial stress state in the soil under self-weight.

Dynamic Analysis: Simulates the response of the soil-structure system to an input ground motion.

How to Use
1. Prerequisites
Ensure you have Python and the following libraries installed:

openseespy

numpy

matplotlib (optional, for post-processing)

opsvis (optional, for visualization)

You can install them using pip:

pip install openseespy numpy matplotlib opsvis

2. Running the Analysis
Clone the Repository:

git clone <your-repository-url>
cd <your-repository-directory>

Input Ground Motion:
Place a ground motion file named velocityHistory.out in the root directory. This file should be a single column of velocity values recorded at each time step.
Note: If this file is not present, the script will create a dummy file to avoid errors, but the analysis results will not be meaningful.

Execute the Script:
Run the script from your terminal:

python OOP.py

Review Outputs:
The script will create an Outputs/ directory containing subdirectories for GravityAnalysis and PostGravityAnalysis, where recorder data (e.g., node accelerations, foundation reactions) will be saved.

3. Modifying Parameters
The script's object-oriented structure makes it easy to change key parameters. All primary inputs are defined within the MaterialParams, StructureParams, and AnalysisParams classes at the top of OOP.py.

Example: Changing the Element Size to 1m x 1m
To override the automatic mesh generation for a specific test case, directly modify the calculate_mesh_parameters method within the MeshParams class:

# Inside the MeshParams class in OOP.py
def calculate_mesh_parameters(self, mat_params):
    """Calculate optimal mesh parameters based on wave propagation requirements"""
    
    # --- MODIFICATION FOR 1x1 ELEMENT TEST CASE ---
    # Directly set the number of elements to match the domain dimensions
    self.numx = int(self.soil_width)
    self.numy = int(self.soil_depth)
    
    # Calculate actual element sizes (this will now result in 1.0)
    self.size_ele_y = self.soil_depth / self.numy
    self.size_ele_x = self.soil_width / self.numx
    # -------------------- END MODIFICATION --------------------

    # ... (rest of the function)

Code Structure
The script is organized into three main parts:

Parameter Classes (MaterialParams, StructureParams, etc.): These classes group related input parameters for materials, structure geometry, analysis settings, and mesh properties. This makes the code clean and easy to manage.

Core Functions (create_soil_mesh, run_gravity_analysis, etc.): Functional blocks that handle specific stages of the modeling and analysis process, such as:

Mesh and Boundary Condition setup

Material and Element Definition

Structure Definition and Connection

Analysis and Output Recording

Main Execution Block: The main() function orchestrates the entire simulation by calling the core functions in the correct sequence. The if __name__ == "__main__": construct ensures that main() is called only when the script is executed directly.

License
This project is licensed under the MIT License - see the LICENSE file for details.
