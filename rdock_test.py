from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess, sys, glob, os, traceback
import string
import random
import tarfile
from scipy.optimize import leastsq
import numpy as np

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

import shapely, descartes
from descartes.patch import PolygonPatch
from shapely.geometry import Polygon, MultiPoint, Point



#----parameter & preparation
def dok2pdb(rdid, reconstructing_pdbtarget, n=0, base_path="outputs/"):
    """
    PRE: TAKES in a rdid with a corresponding mol2 and dok file for the guest
    POST: Will produce the n best complexes that correspond to the docking in pdb and mol2 formats
    """
	
    # ERASE OLD AND SPLITS IN INDIVIDUAL PDB FILES
    for f in glob.glob("{}{}-{}.pdb".format(base_path, rdid,"*")):
         with open(f, 'w') as w: pass
	
    i=0
    #SPLITS THE DOK INTO INDIVIDUAL PDBs
    guest_list = []
    #try:
    with open("{}{}.dok".format(base_path, rdid), 'r') as r:
        for line in r:
            if i>n: break
            curcomplex = "{}{}-{}.pdb".format(base_path, rdid, i)
            with open(curcomplex, "a") as a:
                #a.write(line)
                if "ATOM" in line:
                    pts = line.strip().split()
                    # EXAMPLE:  ATOM      2  C1  LIG     0      -6.550  -3.810  -2.641
                    rec_line = ("{0}{1: >7}  {2: <4}{3}     {4}{5: >12}{6: >8}{7: >8}\n".format(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5], pts[6], pts[7]))
                    #print rec_line
                else:
                    rec_line = line
                a.write(rec_line)
            if "END" in line:
                guest_list.append(curcomplex)
                i=i+1
	
    # CREATES THE COMPLEXES
    cb = Chem.MolFromPDBFile(reconstructing_pdbtarget, removeHs=False) #HETATM tags rather than ATOM are needed for safe reading
    complex_list = []
    for f in sorted(glob.glob("{}{}-{}.pdb".format(base_path, rdid, "*"))):
        print(f)
        try:
            guest = Chem.MolFromPDBFile(f, removeHs=False)
            cbguest = Chem.CombineMols(cb, guest)
            complex_pdb_fname = f.replace(rdid, rdid+'-CB')
            Chem.MolToPDBFile(cbguest, complex_pdb_fname)
            complex_list.append(complex_pdb_fname)
        except:
            print("Complex reassembling failed")
    return 

def make_docking_script(rdid, docking_targetPDB):
    """

    Parameters
    ----------
    docking_targetPDB : a pdb file that is the docking target for the docking procedure
        DESCRIPTION.
    rdid: A random string that defines the molecule and upon which later names will be based on 

    Returns
    -------
    None.
    """
    docking_script="""
Receptor
{0}
RMSD
1.0
Binding pocket
-20 20
-20 20
-20 20
Number of binding poses
20
Ligands list
outputs/{1}-ligands.list
    """.format(docking_targetPDB, rdid)
	
    with open("outputs/{}-ledock.config".format(rdid), "w") as w:
        w.write(docking_script)
    with open("outputs/{}-ligands.list".format(rdid), "w") as w:
        w.write("outputs/{}.mol2".format(rdid))
    return 

def tar_it_all(rdid):
    """
    rdid : a molecule id value
    will tar all related files to outputs/rdid-ALL-SAVED.tar
    """
    os.chdir("outputs")
    with tarfile.open("{}-ALL-SAVED.tar".format(rdid), "w:gz") as tar:
        for el in list(set(glob.glob("{}*".format(rdid)))-set(glob.glob("*tar"))):
            tar.add(el)
            os.remove(el)
    os.chdir("..")
    return 

def find_best_fitting_plane(XYZ):
    """
    Parameters
    ----------
    XYZ : array of xyz coordinates
        the xyz coordinates of the O forming the portal

    Returns
    -------
    the plane best approximating the portal as defined by the Oxygen atoms

    """
    # Inital guess of the plane, the plane should be close to perpendicular to x and pass near 3 Angstrom positive
    XYZ
    p0 = [1, 0, 0, -3]
    
    def f_min(X,p):
        plane_xyz = p[0:3]
        distance = (plane_xyz*X.T).sum(axis=1) + p[3]
        return distance / np.linalg.norm(plane_xyz)
    
    def residuals(params, signal, X):
        return f_min(X, params)
    
    
    sol = leastsq(residuals, p0, args=(None, XYZ))[0]
    
# =============================================================================
#     print("Solution: ", sol/sol[0])
#     print("Old Error: ", (f_min(XYZ, p0)**2).sum())
#     print("New Error: ", (f_min(XYZ, sol)**2).sum())
# =============================================================================
    return sol/sol[0]

def project_xyz_in_plane(xyz, plane_vec):
    """
    Parameters
    ----------
    xyz : list
        points to be projected in the plane
    plane_vec : TYPE
        four parameters defining the plane

    Returns
    -------
    The points projected in the plane using the formula p_proj = p - (n dot p + d) * n 

    """
    return  np.array([x-(np.dot(x, plane_vec[0:3])+plane_vec[3]) * plane_vec[0:3] for x in xyz.T]).T

def translate_rotate_points_to_XY_plane(points_to_transform, plane_vec):
    """
    Parameters
    ----------
    points_to_transform : list of points to convert to the XY plane, these points are ALREADY contained in the plane_vec plane
        DESCRIPTION.
    plane_vec : List
        Parameters describing the planes
    Returns
    -------
    None.
    """
    a, b, c, d = plane_vec
    translate_vector = [0,0,-d/c]
    
    cost = c/(a**2+b**2+c**2)**.5
    sint = ((a**2+b**2)/(a**2+b**2+c**2))**.5
    u1 = b/(a**2+b**2+c**2)**.5
    u2 = -a/(a**2+b**2+c**2)**.5

    rotate_matrix = [[cost + u1**2*(1-cost), u1*u2*(1-cost), u2*sint],
                     [u1*u2*(1-cost), cost+u2**2*(1-cost), -u1*sint],
                     [-u2*sint, u1*sint, cost]]
    
    rotate_matrix = np.array([np.array(x) for x in rotate_matrix])
    
    transformed_points = np.array([np.dot(rotate_matrix, x+translate_vector) for x in points_to_transform.T]).T
    return transformed_points

def get_coverage_score(rdid, base_path="outputs/"):
    """
    Parameters
    ----------
    rdid : str
        an rdid, there must be a file rdid-CB-0.pdb
        The complex is aligned so that the planes containing the oxygens from the CB portals are exclusively contained in the yz plane. One portal has O coordinates in the positive x domain and the other portal in negative x domain

    Returns
    -------
    the coverage number defined as the percentage of occulted portal. That is total area covered by cap atoms divided by portal surface
    """
    # Loads the complex molecule
    capped_complex = Chem.MolFromPDBFile("{}{}-CB-0.pdb".format(base_path, rdid), removeHs=False)
    components = Chem.GetMolFrags(capped_complex, asMols=True) # By default we will get them in the order of appearance in the pdb file, so CB, ADA, CAP
    O_upper_portal_coords = np.array([xyz for sym, xyz in zip([at.GetSymbol() for at in components[0].GetAtoms()], components[0].GetConformer(0).GetPositions()) if (sym=="O" and xyz[0]>0)]).T # get the xyz of the upper portal O
    cap_coords = components[2].GetConformer(0).GetPositions().T # get the cap coordinates
    portal_plane = find_best_fitting_plane(np.array(O_upper_portal_coords)) # finds the best fitting plane for the portal, should be close to (1,0,0,-3)
    
    O_upper_portal_on_plane = project_xyz_in_plane(O_upper_portal_coords,portal_plane) # Get the projection of O coordinates in the new plane
    cap_atoms_projected_on_plane = project_xyz_in_plane(cap_coords,portal_plane) # Get the projection of caps atoms on the portal plane
    
    xy_only_portal_points = translate_rotate_points_to_XY_plane(O_upper_portal_on_plane, portal_plane) # rotate translate portal planes to eliminate z contribution
    xy_only_cap_points = translate_rotate_points_to_XY_plane(cap_atoms_projected_on_plane, portal_plane) # rotate translate caps planes to eliminate z contribution
    
    #print(xy_only_portal_points, xy_only_cap_points, np.mean(xy_only_portal_points, axis=1))
    
    center = np.array(np.mean(xy_only_portal_points, axis=1)) # centering data based on the portal position
    xy_only_portal_points, xy_only_cap_points = np.array([x - c for x,c in zip(xy_only_portal_points, center)]), np.array([x - c for x,c in zip(xy_only_cap_points, center)]) # same, centering data based on the portal position
    
    portal_polygon = MultiPoint([tuple(x) for x in xy_only_portal_points[:2].T]).convex_hull # defines the portal area using the complex hull
    portal_area = portal_polygon.area
    
    occluded = Point(0,0).buffer(0)
    out = Point(0,0).buffer(0)

    for pts in [tuple(x) for x in xy_only_cap_points[:2].T]:
        a = Point(*pts)
        if a.within(portal_polygon):
            occluded = occluded.union(a.buffer(1.3))
        else:
            out = out.union(a.buffer(1.3))    
    portal_fraction_occluded = portal_polygon.intersection(occluded.union(out))
    
# =============================================================================
#     # VISULALIZE  THE AREAS
#     fig = plt.figure(1, dpi=90)
#     ax = fig.add_subplot(111)
#     ax.plot(*portal_polygon.exterior.xy)
#     ax.add_patch(PolygonPatch(occluded, alpha=0.2, zorder=2)) #ATOM IN ABOVE THE PORTAL
#     ax.add_patch(PolygonPatch(out, alpha=0.1, zorder=2)) # ATOMS OUT ABOVE THE PORTAL
#     ax.add_patch(PolygonPatch(portal_fraction_occluded, alpha=0.5, zorder=2)) # PORTAL intersection ATOMS ABOVE
#     plt.show()
#     
#     # PLOTTING TESTS 
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(*O_upper_portal_coords)
#     ax.scatter(*O_upper_portal_on_plane)
#     ax.scatter(*cap_coords)
#     ax.scatter(*cap_atoms_projected_on_plane)
#     ax.scatter(*xy_only_portal_points)
#     ax.scatter(*xy_only_cap_points)
#     plt.axis('tight')
#     ax.set_xlabel('X Label')
#     ax.set_ylabel('Y Label')
#     ax.set_zlabel('Z Label')
#     plt.show()
# =============================================================================

    
    # Find portal surface and plane
    return portal_fraction_occluded.area/portal_area
    
def rdock_score(compound, docking_target="adamantanone-docked-named-c.pdb"):

    input_smiles = str(compound)
    #num_docking = 3 # number of docking trials/ not an option here, hardcoded in the config file
	
    rdstring = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(32)])
    mol2_file = 'outputs/{}.mol2'.format(rdstring)
    sdf_file = 'outputs/{}.sdf'.format(rdstring)
    min_score = 10**10

    #----Translation from SMILES to sdf
    fw = Chem.SDWriter(sdf_file)
    m1 = Chem.MolFromSmiles(input_smiles)
    try:
        if m1!= None:
            m = Chem.AddHs(m1)

            cid = AllChem.EmbedMolecule(m)
            #fw.write(m)

            opt = AllChem.MMFFOptimizeMolecule(m,maxIters=2000)
            print(opt)

            fw.write(m)
            fw.close()
            
            #----convert to mol2
            convert_cmd = "obabel {} -O {}".format(sdf_file, mol2_file)
            proc = subprocess.call(convert_cmd.split() , shell=False)
            
            #----rdock calculation
            make_docking_script(rdstring, docking_target)
            cmd = "./ledock_linux_x86 outputs/{0}-ledock.config > outputs/{0}.log".format(rdstring)
            proc = subprocess.check_output(cmd.split() , shell=False )
            #proc = subprocess.call(cmd.split() , shell=False )
            #proc = subprocess.call(cmd , shell=True )
            
            dok2pdb(rdstring, docking_target, n=0)
            
            #----find the percentage of portal occluded
            percent_portal_covered = get_coverage_score(rdstring)
            
            
            #----find the minimum score of rdock from multiple docking results
            with open("outputs/{}.dok".format(rdstring), "r") as r:
                lines = r.readlines()
            


            min_score = float(lines[1].split()[-2])



            print('minimum rdock score', min_score)
            #----make backup
# =============================================================================
#             with open("outputs/{}-RES".format(rdstring), "w") as w:
#                 w.write("{}\t{}\t{}\n".format(rdstring, compound, min_score))
# =============================================================================
				
            tar_it_all(rdstring)
    except:
        print(traceback.format_exc())
        min_score=10**10
        percent_portal_covered=0



    return min_score, percent_portal_covered, rdstring

if __name__=="__main__":
# =============================================================================
#     import sys, glob
#     for f in glob.glob("/home/macenrola/Documents/ML/ChemTS/new_scoring_for_mcts/all_the_sdf_named_for_sale/sulfonamides_docked/*.dok"):
#         dok2pdb(f[:-4], reconstructing_pdbtarget="/home/macenrola/Documents/ML/ChemTS/new_scoring_for_mcts/all_the_sdf_named_for_sale/sulfonamides_docked/adamantanone-docked-named.pdb", base_path="", n=5)
#         print(f)
#         print(f[:-6]+"-CB-0.pdb")
# 
# =============================================================================

   # get_coverage_score("2BlgxNRGF5ko9PoiYXYO1VQCV4mUfgZw", base_path="/home/macenrola/Documents/ML/ChemTS/ledock_ligand_design_with_sa_scorer_adaincl_good_save/outputs/")

    pass
