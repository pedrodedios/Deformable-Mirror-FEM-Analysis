import FreeCAD as App
import Part
import Fem
import ObjectsFem
import numpy as np
import math
from femtools import ccxtools

import numpy as np
import math

def draw_mirror(
    doc,
    r_mech,
    r_mirror,
    h_mirror,
    r1_plunger,
    r2_plunger,
    h_plunger,
    positions):
    """
    Generate a parametric deformable mirror with plungers in FreeCAD.

    This function creates a full mirror geometry including:
    1. The main mirror disk (deformable area)
    2. The mechanical outer ring (fixed ring)
    3. A set of circular (or elliptical) plungers at specified positions 
    4. Fusion of all components into a single solid

    It is designed to support FEM analysis and actuator influence function computation.

    Parameters
    ----------
    doc : FreeCAD.Document
        The active FreeCAD document where the geometry will be created.
    r_mech : float
        Outer mechanical radius of the mirror assembly.
    r_mirror : float
        Radius of the mirror surface.
    h_mirror : float
        Thickness of the mirror disk.
    r1_plunger : float
        Major radius of the elliptical plunger.
    r2_plunger : float
        Minor radius of the elliptical plunger.
    h_plunger : float
        Extrusion height of the plunger.
    positions : list of [x, y]
        XY coordinates of the plungers relative to the mirror center.

    Returns
    -------
    None
        The function directly modifies the FreeCAD document by adding
        the mirror and plunger geometry.

    Notes
    -----
    - Each plunger consists of a support cylinder and an elliptical extrusion
      that are fused together.
    - Fillets are applied to plunger edges for smoother geometry.
    - The final fusion object represents the complete mirror assembly and
      can be used as the base shape for FEM simulations.
    - This function does not perform FEM analysis itself; it only generates geometry.

    Example
    -------
    >>> import FreeCAD as App
    >>> doc = App.newDocument("Mirror")
    >>> positions = [[0,0],[32,32],[-32,32],[-32,-32],[32,-32]]
    >>> draw_mirror(doc, r_mech=100, r_mirror=95, h_mirror=2,
    ...             r1_plunger=1.2, r2_plunger=1.2, h_plunger=4, positions=positions)
    """


    radio_fillet = 1.2
    design_objects = []

    # --- Mirror disk ---
    mirror = doc.addObject("Part::Cylinder", "Mirror")
    mirror.Radius = r_mirror
    mirror.Height = h_mirror
    design_objects.append(mirror)

    # --- Mechanical ring ---
    line = doc.addObject("Part::Line", "RingProfile")
    line.X1, line.Y1, line.Z1 = r_mirror, 0, 0
    line.X2, line.Y2, line.Z2 = r_mech, 0, 0

    revolve = doc.addObject("Part::Revolution", "RingRevolve")
    revolve.Source = line
    revolve.Axis = (0, 0, 1)
    revolve.Angle = 360
    revolve.Solid = False
    revolve.Visibility = False

    ring = doc.addObject("Part::Extrusion", "Ring")
    ring.Base = revolve
    ring.LengthFwd = 2 * h_mirror
    ring.Solid = False
    design_objects.append(ring)

    # --- Plungers ---
    for i, (x, y) in enumerate(positions):

        angle = math.degrees(math.atan2(y, x))

        support_cyl = doc.addObject("Part::Cylinder", f"PlungerSupport{i}")
        support_cyl.Radius = 4 * r1_plunger
        support_cyl.Height = h_mirror / 2
        support_cyl.Placement = App.Placement(
            App.Vector(x, y, h_mirror / 2),
            App.Rotation(0, 0, 0)
        )

        ellipse = doc.addObject("Part::Ellipse", f"PlungerProfile{i}")
        ellipse.MajorRadius = r1_plunger
        ellipse.MinorRadius = r2_plunger
        ellipse.Placement = App.Placement(
            App.Vector(x, y, h_mirror),
            App.Rotation(App.Vector(0, 0, 1), angle + 90)
        )

        plunger = doc.addObject("Part::Extrusion", f"Plunger{i}")
        plunger.Base = ellipse
        plunger.LengthFwd = h_plunger
        plunger.Solid = True

        fusion = doc.addObject("Part::MultiFuse", f"PlungerFusion{i}")
        fusion.Shapes = [support_cyl, plunger]

        fillet = doc.addObject("Part::Fillet", f"PlungerFillet{i}")
        fillet.Base = fusion
        fillet.Edges = [(4, radio_fillet, 1)]

        design_objects.append(fillet)

    # --- Final fusion ---
    final_fusion = doc.addObject("Part::MultiFuse", "Fusion")
    final_fusion.Shapes = design_objects

    doc.recompute()

    return