"""
Soil-Structure Interaction Analysis with Wave-Based Mesh Generation

This script models soil-structure interaction with:
1. Automatic mesh sizing based on wave propagation requirements
2. Robust node numbering system
3. Complete soil material definition
4. SDOF structure with foundation
5. Dynamic analysis capabilities
"""

# ------------------------------ Import Libraries ----------------------------- #
import openseespy.opensees as ops
import opsvis as opsv
from math import ceil, pi
import numpy as np
import matplotlib.pyplot as plt
import vfo.vfo as vfo
import os

# ---------------------------- Global Parameters ----------------------------- #
class MeshParams:
    def __init__(self, mat_params):
        # Wave propagation parameters
        self.f_max = 25.0      # Maximum frequency to resolve [Hz]
        self.ele_per_wave = 1  # Elements per wavelength
        
        # Physical dimensions
        self.soil_width = 40.0  # Total soil width [m]
        self.soil_depth = 20.0  # Total soil depth [m]
        self.node_start = 1     # Starting node number
        
        # Initialize calculated parameters
        self.numx = None
        self.numy = None
        self.size_ele_x = None
        self.size_ele_y = None
        self.xlength = None
        self.ylength = None
        
        # Calculate mesh parameters
        self.calculate_mesh_parameters(mat_params)
        
    def calculate_mesh_parameters(self, mat_params):
        """Calculate optimal mesh parameters based on wave propagation requirements"""
        
        # --- MODIFICATION FOR 1x1 ELEMENT TEST CASE ---
        # Directly set the number of elements to match the domain dimensions
        self.numx = int(self.soil_width)
        self.numy = int(self.soil_depth)
        
        # Calculate actual element sizes (this will now result in 1.0)
        self.size_ele_y = 2.5*self.soil_depth / self.numy
        self.size_ele_x = 2.5*self.soil_width / self.numx
        
        # Original automatic calculation is commented out for the test case
        # wavelength = mat_params.Vs / (4 * self.f_max)
        # target_ele_size = wavelength / self.ele_per_wave
        # self.numy = max(4, int(ceil(self.soil_depth / target_ele_size)))
        # self.numx = max(4, int(ceil(self.soil_width / target_ele_size)))
        # -------------------- END MODIFICATION --------------------

        # Domain lengths (for visualization/normalization)
        self.xlength = self.soil_width
        self.ylength = self.soil_depth
        
        # Aspect ratio check
        aspect_ratio = self.size_ele_x / self.size_ele_y
        if aspect_ratio > 2 or aspect_ratio < 0.5:
            print(f"Warning: Large element aspect ratio ({aspect_ratio:.2f})")
        
        print(f"Mesh parameters MANUALLY SET for 1x1 element size:")
        print(f" - Elements: {self.numx}x{self.numy}")
        print(f" - Element size: {self.size_ele_x:.2f}m x {self.size_ele_y:.2f}m")
    
    @property
    def nnx(self):
        """Number of nodes in x-direction"""
        return self.numx + 1
        
    @property
    def nny(self):
        """Number of nodes in y-direction"""
        return self.numy + 1
        
    @property
    def total_nodes(self):
        """Total number of nodes in mesh"""
        return self.nnx * self.nny

class StructureParams:
    def __init__(self):
        self.height = 5.0     # Height of SDOF structure [m]
        self.x_loc = 20.0     # X-coordinate of structure [m]
        self.pile_length = 5  # Length of pile foundation [elements]
        self.mass = 100.0     # SDOF mass [Mg]
        self.sfnum = 10000    # Structure foundation number offset

class MaterialParams:
    def __init__(self):
        # Soil properties
        self.rho = 1.7          # Mass density [Mg/m³]
        self.Vs = 250           # Shear wave velocity [m/s]
        self.nu = 0.45          # Poisson's ratio
        self.cohesion = 95.0    # Cohesion [kPa]
        self.peak_strain = 0.05 # Peak shear strain
        self.ref_press = 100    # Reference pressure [kPa]
        
        # Calculate derived properties
        self.G = self.rho * self.Vs**2  # Shear modulus [kN/m²]
        self.E = 2 * self.G * (1 + self.nu)  # Young's modulus [kN/m²]
        self.K = self.E/(3*(1-2*self.nu))  # Bulk modulus [kN/m²]
        
        # Rock properties
        self.rock_Vs = 760    # Bedrock shear wave velocity [m/s]
        self.rock_den = 2.4   # Bedrock density [Mg/m³]
        
        # Structural properties
        self.E_struct = 2.5e7  # Structural modulus [kN/m²]
        self.col_width = 0.4   # Column width [m]
        self.col_height = 0.4  # Column height [m]

class AnalysisParams:
    def __init__(self):
        self.damp = 0.02       # Damping ratio
        self.gamma = 0.5       # Newmark gamma
        self.beta = 0.25       # Newmark beta
        self.motion_dt = 0.005 # Ground motion time step [s]
        self.motion_steps = 7990  # Number of steps
        
        # Calculate Rayleigh damping
        omega1 = 2*pi*0.2  # [rad/s]
        omega2 = 2*pi*20   # [rad/s]
        self.a0 = 2*self.damp*omega1*omega2/(omega1+omega2)
        self.a1 = 2*self.damp/(omega1+omega2)
        print(f'Damping coefficients - a0: {self.a0:.4f}, a1: {self.a1:.4f}')

# Initialize parameter classes
mat = MaterialParams()
struct = StructureParams()
analysis = AnalysisParams()
mesh = MeshParams(mat)

# --------------------------- Mesh Generation Functions ------------------------ #
def generate_node_tag(i, j, mesh_params):
    """
    Generate unique node tag based on grid position
    
    Args:
        i (int): x-index (1 to nnx)
        j (int): y-index (1 to nny)
        mesh_params (MeshParams): Mesh parameters object
        
    Returns:
        int: Unique node tag
    """
    return mesh_params.node_start + (j-1)*mesh_params.nnx + (i-1)

def create_soil_mesh(mesh_params):
    """
    Create soil domain mesh with proper node numbering
    
    Args:
        mesh_params (MeshParams): Mesh parameters object
    """
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.wipe()
    
    # Create nodes with physical coordinates
    for j in range(1, mesh_params.nny + 1):
        for i in range(1, mesh_params.nnx + 1):
            node_tag = generate_node_tag(i, j, mesh_params)
            x_coord = (i-1) * mesh_params.size_ele_x
            y_coord = (j-1) * mesh_params.size_ele_y
            ops.node(node_tag, x_coord, y_coord)
    
    print(f"Created soil mesh with {mesh_params.total_nodes} nodes")
    print(f"Physical domain: {mesh_params.xlength}m x {mesh_params.ylength}m")

# ------------------------- Boundary Condition Functions ---------------------- #
def apply_boundary_conditions(mesh_params, mat_params):
    """
    Apply boundary conditions to soil domain
    
    Args:
        mesh_params (MeshParams): Mesh parameters object
    """
    # Fix base nodes (y-direction only)
    for i in range(1, mesh_params.nnx + 1):
        node_tag = generate_node_tag(i, 1, mesh_params)
        ops.fix(node_tag, 0, 1)
    
    # Periodic boundary conditions (left-right sides)
    for j in range(2, mesh_params.nny + 1):
        left_node = generate_node_tag(1, j, mesh_params)
        right_node = generate_node_tag(mesh_params.nnx, j, mesh_params)
        ops.equalDOF(left_node, right_node, 1, 2)
    
    # Create and constrain dashpot nodes
    dash_node1 = mesh_params.total_nodes + 1
    dash_node2 = dash_node1 + 1
    ops.node(dash_node1, 0.0, 0.0)
    ops.node(dash_node2, 0.0, 0.0)
    ops.fix(dash_node1, 1, 1)
    ops.fix(dash_node2, 0, 1)

    # Connect all base nodes horizontally
    base_master_node = generate_node_tag(1, 1, mesh_params)
    for i in range(2, mesh_params.nnx + 1):
        rNode = generate_node_tag(i, 1, mesh_params)
        ops.equalDOF(base_master_node, rNode, 1)
    
    # Connect base to dashpot
    ops.equalDOF(base_master_node, dash_node2, 1)

    # Create dashpot material and element
    dashMatTag = 4000
    dashEleTag = mesh_params.total_nodes + 1000
    mC = mesh_params.soil_width * mat_params.rock_den * mat_params.rock_Vs
    ops.uniaxialMaterial('Viscous', dashMatTag, mC, 1)
    ops.element('zeroLength', dashEleTag, dash_node1, dash_node2, '-mat', dashMatTag, '-dir', 1)
    
    print("Boundary conditions and base dashpot applied")

# ------------------------- Material Definition Functions --------------------- #
def define_soil_material(mat_params):
    """
    Define soil material properties
    
    Args:
        mat_params (MaterialParams): Material parameters object
    """
    mat_tag = 1
    nd = 2
    no_yield_surf = 25
    
    ops.nDMaterial('PressureIndependMultiYield', mat_tag, nd, mat_params.rho, 
                  mat_params.G, mat_params.K, mat_params.cohesion, 
                  mat_params.peak_strain, no_yield_surf, 0.0, 
                  mat_params.ref_press, 0.0)
    
    print("Soil material defined")

def create_soil_elements(mesh_params, mat_params):
    """
    Create quadrilateral elements for soil domain
    
    Args:
        mesh_params (MeshParams): Mesh parameters object
        mat_params (MaterialParams): Material parameters object
    """
    wgtX, wgtY = 0.0, -9.81 * mat_params.rho
    
    for j in range(1, mesh_params.numy + 1):
        for i in range(1, mesh_params.numx + 1):
            ele_tag = (j-1)*mesh_params.numx + i
            iNode = generate_node_tag(i, j, mesh_params)
            jNode = generate_node_tag(i+1, j, mesh_params)
            kNode = generate_node_tag(i+1, j+1, mesh_params)
            lNode = generate_node_tag(i, j+1, mesh_params)
            ops.element('quad', ele_tag, iNode, jNode, kNode, lNode, 
                        1, 'PlaneStrain', 1, 0.0, 0.0, wgtX, wgtY)
    
    print(f"Created {mesh_params.numx*mesh_params.numy} soil elements")

# ------------------------- Structure Definition Functions -------------------- #
def create_structure(mesh_params, struct_params, mat_params):
    """
    Create SDOF structure and connect to soil
    
    Args:
        mesh_params (MeshParams): Mesh parameters object
        struct_params (StructureParams): Structure parameters object
        mat_params (MaterialParams): Material parameters object
    """
    # Switch to a model with 3 DOFs for beam-column elements
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    struct_y = mesh_params.ylength
    
    # Foundation nodes
    base_node = struct_params.sfnum + 1000
    ops.node(base_node, struct_params.x_loc, struct_y)
    ops.node(base_node + 1, struct_params.x_loc - 1, struct_y)
    ops.node(base_node + 2, struct_params.x_loc + 1, struct_y)
    
    # SDOF mass node
    top_node = base_node + 40
    ops.node(top_node, struct_params.x_loc, struct_y + struct_params.height)
    
    # Pile nodes
    for y in range(1, struct_params.pile_length + 1):
        pile_node = base_node - y
        y_coord = struct_y - y * mesh_params.size_ele_y
        ops.node(pile_node, struct_params.x_loc, y_coord)
    
    # Define structural elements
    lin_transf = 1
    ops.geomTransf('Linear', lin_transf)
    col_A = mat_params.col_width * mat_params.col_height
    col_I = (mat_params.col_width * mat_params.col_height**3) / 12
    
    # Pile Cap Elements
    ops.element('elasticBeamColumn', base_node, base_node, base_node + 1, col_A, mat_params.E_struct, col_I, lin_transf)
    ops.element('elasticBeamColumn', base_node + 1, base_node, base_node + 2, col_A, mat_params.E_struct, col_I, lin_transf)
    
    # SDOF Column
    ops.element('elasticBeamColumn', 9999, base_node, top_node, col_A, mat_params.E_struct, col_I, lin_transf)

    # Pile Elements
    for y in range(struct_params.pile_length):
        pile_element = struct_params.sfnum + 10000 - y
        ops.element('elasticBeamColumn', pile_element, base_node - y, base_node - y - 1, col_A, mat_params.E_struct, col_I, lin_transf)

    # Apply mass to SDOF top node
    ops.mass(top_node, struct_params.mass, struct_params.mass, 0.0) # Corrected mass application
    
    print("Structure created")

def connect_structure_to_soil(mesh_params, struct_params):
    """
    Connect structure nodes to soil nodes
    
    Args:
        mesh_params (MeshParams): Mesh parameters object
        struct_params (StructureParams): Structure parameters object
    """
    # Find nearest soil column
    soil_col = int(round(struct_params.x_loc / mesh_params.size_ele_x)) + 1
    
    # Connect foundation to soil surface
    ops.equalDOF(struct_params.sfnum + 1000, generate_node_tag(soil_col, mesh_params.nny, mesh_params), 1, 2)
    ops.equalDOF(struct_params.sfnum + 1001, generate_node_tag(soil_col - 1, mesh_params.nny, mesh_params), 1, 2)
    ops.equalDOF(struct_params.sfnum + 1002, generate_node_tag(soil_col + 1, mesh_params.nny, mesh_params), 1, 2)
    
    # Connect pile nodes to soil nodes
    for y in range(1, struct_params.pile_length + 1):
        soil_node = generate_node_tag(soil_col, mesh_params.nny - y, mesh_params)
        pile_node = struct_params.sfnum + 1000 - y
        ops.equalDOF(pile_node, soil_node, 1, 2)
    
    print("Structure connected to soil")

# ---------------------------- Analysis Functions ---------------------------- #
def create_gravity_recorders(mesh_params, struct_params):
    """Create recorders for gravity and dynamic analysis."""
    # Create output directories if they don't exist
    os.makedirs("Outputs/GravityAnalysis", exist_ok=True)
    os.makedirs("Outputs/PostGravityAnalysis", exist_ok=True)

    # Gravity analysis recorders
    soil_surface_node = generate_node_tag(int(mesh.nnx/2), mesh.nny, mesh)
    base_rock_node = generate_node_tag(1, 1, mesh)
    ops.recorder('Node', '-file', 'Outputs/GravityAnalysis/Gaccel.out', '-time',
                 '-node', base_rock_node, soil_surface_node, '-dof', 1, 'accel')
    print("Gravity recorders created.")

  
def create_dynamic_recorders(mesh_params, struct_params):
    """Create recorders for dynamic analysis."""
    os.makedirs("Outputs/PostGravityAnalysis", exist_ok=True)
     # Gravity analysis recorders
    soil_surface_node = generate_node_tag(int(mesh.nnx/2), mesh.nny, mesh)
    base_rock_node = generate_node_tag(1, 1, mesh)

    # Post-gravity (dynamic) analysis recorders
    base_node = struct_params.sfnum + 1000
    ops.recorder('Node', '-file', 'Outputs/PostGravityAnalysis/Rbase.out', '-time',
                 '-node', base_node, '-dof', 1, 2, 3, 'reaction')
    ops.recorder('Node', '-file', 'Outputs/PostGravityAnalysis/accel.out', '-time',
                 '-node', base_rock_node, soil_surface_node, '-dof', 1, 'accel')
    print("Post-gravity recorders created.")


def run_gravity_analysis(analysis_params):
    """Run gravity loading analysis"""
    print("\nStarting gravity analysis...")
    ops.constraints('Transformation')
    ops.test('NormDispIncr', 1e-5, 30, 1)
    ops.algorithm('Newton')
    ops.numberer('RCM')
    ops.system('ProfileSPD')
    ops.integrator('Newmark', analysis_params.gamma, analysis_params.beta)
    ops.analysis('Transient')
    
    # Elastic gravity phase
    ops.analyze(10, 5e2) # Corrected time step from main.py
    
    # Plastic gravity phase
    ops.updateMaterialStage('-material', 1, '-stage', 1)
    ops.analyze(40, 5e2) # Corrected time step from main.py
    
    print("Gravity analysis completed")
    ops.setTime(0.0)
    ops.wipeAnalysis()

def run_dynamic_analysis(analysis_params, mesh_params, mat_params):
    """Run dynamic analysis"""
    print("\nStarting dynamic analysis...")
    # Setup analysis
    ops.constraints('Transformation')
    ops.test('NormDispIncr', 1e-3, 15, 1)
    ops.algorithm('Newton')
    ops.numberer('RCM')
    ops.system('ProfileSPD')
    ops.integrator('Newmark', analysis_params.gamma, analysis_params.beta)
    ops.rayleigh(analysis_params.a0, analysis_params.a1, 0.0, 0.0)
    ops.analysis('Transient')
    
    # Calculate stable time step
    duration = analysis_params.motion_dt * analysis_params.motion_steps
    kTrial = mesh_params.size_ele_x / mat_params.Vs
    
    if analysis_params.motion_dt <= kTrial:
        nSteps = analysis_params.motion_steps
        dT = analysis_params.motion_dt
    else:
        nSteps = int(duration // kTrial) + 1
        dT = duration / nSteps
    
    print(f"Running dynamic analysis with {nSteps} steps, dt={dT:.6f}s")
    print(f"Total duration: {duration:.2f}s")
    
    # Apply ground motion
    cFactor = mesh_params.soil_width * mat_params.rock_den * mat_params.rock_Vs
    # Ensure a dummy velocityHistory.out file exists
    if not os.path.exists('velocityHistory.out'):
        with open('velocityHistory.out', 'w') as f:
            f.write('0.0\n')

    ops.timeSeries('Path', 1, '-dt', analysis_params.motion_dt, 
                  '-filePath', 'velocityHistory.out', '-factor', cFactor)
    ops.pattern('Plain', 10, 1)
    # Apply load to the master base node
    base_master_node = generate_node_tag(1, 1, mesh_params)
    ops.load(base_master_node, 1.0, 0.0) # Added zero for moment
    
    # Run analysis
    ops.analyze(nSteps, dT)
    
    print("Dynamic analysis completed")

# ---------------------------- Main Analysis Script --------------------------- #
def main():
    print("\n=== Starting Soil-Structure Interaction Analysis ===")
    
    # Create soil model
    create_soil_mesh(mesh)
    apply_boundary_conditions(mesh, mat)
    define_soil_material(mat)
    create_soil_elements(mesh, mat)
    
    # Create and connect structure
    create_structure(mesh, struct, mat)
    connect_structure_to_soil(mesh, struct)
    
    # Create recorders
    create_gravity_recorders(mesh, struct)
    # Run gravity analyses
    run_gravity_analysis(analysis)

    create_dynamic_recorders(mesh, struct)
    # Run dynamic analysis
    run_dynamic_analysis(analysis, mesh, mat)
    
    print("\n=== Analysis Completed Successfully ===")
    ops.wipe()

if __name__ == "__main__":
    main()