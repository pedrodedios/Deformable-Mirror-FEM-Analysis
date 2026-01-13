"""
Perform FEM analysis of a deformable mirror and compute Zernike mode residuals.

This function combines geometry creation, FEM setup, actuator influence function
computation, and inverse wavefront reconstruction to evaluate how accurately
a deformable mirror can reproduce a set of target Zernike modes.

Steps
-----
1. Generate the mirror and plungers using `mirror_geometry.py`.
2. Set up the FEM analysis in FreeCAD using CalculiX.
3. Apply fixed constraints and assign material properties.
4. Generate the FEM mesh for the mirror assembly.
5. Compute the Actuator Influence Function (AIF) matrix by applying
test forces to each plunger and recording nodal displacements.
6. Compute the pseudo-inverse of the AIF matrix.
7. For each target Zernike mode:
    - Compute required actuator forces using the pseudo-inverse (using zernike_utils.py)
    - Apply forces and solve FEM
    - Compute RMS and peak residual errors
    - Save nodal displacement data to CSV

Parameters
----------
doc : FreeCAD.Document
    Active FreeCAD document for the analysis.
r_mech : float
    Outer mechanical radius of the mirror assembly.
r_mirror : float
    Radius of the mirror surface.
h_mirror : float
    Thickness of the mirror disk.
r1_plunger : float
    Major radius of the elliptical plungers.
r2_plunger : float
    Minor radius of the elliptical plungers.
h_plunger : float
    Height of the plunger extrusion.
mesh_size : float
    Maximum element size for FEM mesh generation.
set_modes_noll : list of [int, float]
    Target Zernike modes to analyze; each element contains:
    - Noll index of the mode
    - Weight (amplitude) of the mode
CA : float
    Radius of the clear aperture for evaluating nodal displacements.
force_test : float
    Test force magnitude used for computing the AIF matrix.
positions : list of [x, y]
    XY coordinates of the plungers relative to the mirror center.

Returns
-------
rms_mean : float
    Mean RMS residual error across all analyzed Zernike modes.
data_modes : list of [int, float, float]
    For each Zernike mode: [Noll index, RMS error, peak error]

Notes
-----
- Nodal displacements and residuals are saved as CSV files in:
  `nodes_displacements_<Noll index>.csv`
- The function removes temporary force constraints after each mode.
- It requires FreeCAD with FEM module and CalculiX (via `ObjectsFem` and `ccxtools`).

"""


import FreeCAD as App
import Part
import Fem
import ObjectsFem
import numpy as np
import math
from femtools import ccxtools

import numpy as np
import math

from mirror_geometry import draw_mirror
from zernike_utils import zernike_radial, noll_to_nm, zernike_noll, rms 
import os





def run_fem_analysis(
    doc,
    r_mech,
    r_mirror,
    h_mirror,
    r1_plunger,
    r2_plunger,
    h_plunger,
    mesh_size,
    set_modes_noll,
    CA,
    force_test,
    positions):
    



    BASE_DIR = os.path.dirname(os.path.abspath(__file__))





    # --------------------
    # Mirror geometry
    # --------------------
    draw_mirror(
        doc,
        r_mech,
        r_mirror,
        h_mirror,
        r1_plunger,
        r2_plunger,
        h_plunger,
        positions
    )
    


    num_plungers = len(positions)

    # --------------------
    # FEM setup
    # --------------------
    analysis = ObjectsFem.makeAnalysis(doc, "Analysis")

    solver = ObjectsFem.makeSolverCalculixCcxTools(doc, "CalculiX")
    solver.GeometricalNonlinearity = "linear"
    solver.ThermoMechSteadyState = True
    solver.MatrixSolverType = "default"
    solver.IterationsControlParameterTimeUse = False
    analysis.addObject(solver)

    material = ObjectsFem.makeMaterialSolid(doc, "Material")
    material.Material = {
        "Name": "Aluminum-6061",
        "YoungsModulus": "69000 MPa",
        "PoissonRatio": "0.33",
        "Density": "2700 kg/m^3",
    }
    analysis.addObject(material)

    #fix constraint on edge of the mirror 
    fixed_constraint = ObjectsFem.makeConstraintFixed(doc, "Fixed")
    fixed_constraint.References = [(doc.Fusion, "Face1")]
    analysis.addObject(fixed_constraint)

    # mesh construction
    mesh = doc.addObject("Fem::FemMeshShapeNetgenObject", "Mesh")
    mesh.Shape = doc.Fusion
    mesh.MaxSize = mesh_size
    mesh.Fineness = "Fine"
    mesh.Optimize = True
    mesh.SecondOrder = True
    mesh.NbSegsPerRadius = 250
    analysis.addObject(mesh)
    
    
    
    # --------------------------------
    # Influence Functions (AIF matrix)
    #--------------------------------

    aif_columns = []

    for i in range(num_plungers): #AIF columns


        plunger = doc.getObject(f"Plunger{i}")
        force = ObjectsFem.makeConstraintForce(doc, f"ForceIF{i}") # force on ith-pin
        force.References = [(plunger, "Face3")]  
        force.Force = force_test
        force.Reversed = True   # push pin
        analysis.addObject(force)

        # run FEM analysis all in one
        from femtools import ccxtools
        fea = ccxtools.FemToolsCcx()
        fea.purge_results()
        fea.run()

        displacement = np.array(doc.CCX_Results.DisplacementVectors)
        
        # mesh nodes on the active surface 
        nodes = np.array([
            [n.x, n.y, n.z]
            for n in doc.Mesh.FemMesh.Nodes.values()
        ])

        # mesh nodes within Clear Aperture
        inner = np.where(
            (nodes[:, 2] == 0.0)
            & ((nodes[:, 0] ** 2 + nodes[:, 1] ** 2) <= CA ** 2)
        )

        # displacements on the z direction
        dz = displacement[inner][:, 2]
        aif_columns.append(dz)


        # Remove current force over pin after used
        if doc.getObject(f"ForceIF{i}"):
            doc.removeObject(f"ForceIF{i}")

    # --------------------
    # AIF Matrix & Geometry
    # --------------------
    AIF = np.column_stack(aif_columns)
    AIF_pinv = np.linalg.pinv(AIF)  #AIF pseudoinverse

    nodes = nodes[inner]
    x = nodes[:, 0] / CA
    y = nodes[:, 1] / CA

    z_R = np.sqrt(x ** 2 + y ** 2)
    z_Phi = np.arctan2(y, x)





    rms_values = []
    all_forces = []
    data_modes = []



   # --------------------------------
   # Inverse wavefront reconstruction
   # --------------------------------
   # Solve for actuator forces by projecting target Zernike modes
   # onto the mirror actuator influence-function basis (AIF matrix)
    for mode in set_modes_noll:

        noll_index = int(mode[0])
        weight_mode = mode[1]
        displacements = []
        nodes_displacements = []

        #Evaluate inner nodes (clear aperture) on  target Zernike mode
        zernike_target = zernike_noll(
                         R      = z_R,
                         Phi    = z_Phi,
                         j      = noll_index,
                         weight = weight_mode
                          )	


        forces = np.matmul(AIF_pinv, zernike_target) #pins's forces
        all_forces.append(forces)


        # Apply forces on pins
        for i in range(0,num_plungers):
            plunger = doc.getObject(f"Plunger{i}")
            force_plunger = ObjectsFem.makeConstraintForce(doc, f"ForceMode{i}")
            force_plunger.References = [(plunger, "Face3")]
            force_plunger.Force = forces[i]
            force_plunger.Reversed = True
            analysis.addObject(force_plunger)

        # run FEM analysis
        from femtools import ccxtools
        fea = ccxtools.FemToolsCcx()
        fea.purge_results()
        fea.run()

        displacement = np.array(doc.CCX_Results.DisplacementVectors)
        dz = displacement[inner][:, 2]

        error = zernike_target - dz
        rms_mode = rms(error)
        peak_error = np.max(error)

        rms_values.append(rms_mode)
        data_modes.append([noll_index, rms_mode, peak_error])

        nodes_displacements = np.column_stack((x, y, zernike_target, dz))

        #save node displacements from Clear Aperture 
        filename = f"nodes_displacements_{noll_index}.csv"
        filepath = os.path.join(BASE_DIR, filename)

        np.savetxt(filepath, nodes_displacements, delimiter=',')
        
        

        for i in range(num_plungers):
            doc.removeObject(f"ForceMode{i}")

    rms_mean = float(np.mean(rms_values))


    for noll, rms, peak in data_modes:
    print(f"Noll {noll}: RMS = {rms:.3e}, Peak = {peak:.3e}")



    return rms_mean, data_modes