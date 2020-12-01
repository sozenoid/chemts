#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 16:18:59 2020

@author: macenrola
"""
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, Descriptors 
import random, string
import subprocess, os, glob, sys
import numpy as np
import pickle 
import time
import tarfile
adamantanoneMOL2, CBMOL2, docking_targetPDB, cp2k_opti_file, apbs_inp = "docking_targets/adamantanone-GOOD.mol2", "docking_targets/CB7-GOOD.mol2", "docking_targets/adamantanone-docked-named-opti.pdb", "opti_vib_cp2k.inp", "apbs_inp"
Ha2kcal=627.5
binary_complex_values = {"E_tot":-0.6466416447*Ha2kcal, "S_tot": 367.730, "E_CPCM":-0.790020704373*Ha2kcal} #all kcal/mol except S_tot in cal/mol/K USING PM7 Int(Grid=SG1) and loose opti

#WORKSATION
wdir = "/media/macenrola/HUGO_+6590843639/Chemts/ledock_ligand_design_with_SMALL_MOLS_NEW_SCORING_RST"
obabel_path = "/usr/bin/obabel"
ledock_path = "" # Ledock should be in the same folder as wdir
#apbs_path = "/usr/bin/apbs"
#antechamber_path = "/home/macenrola/anaconda3/envs/chemts/bin"
#cp2k_path = "/home/macenrola/anaconda3/pkgs/cp2k-6.1.0-hc6cf775_3/bin/cp2k.sopt" # double check that cp2k is executing on a single core as it should
#g16_path = "/home/macenrola/Documents/Gaussian/g16/g16"
xtb_path = "/home/macenrola/xtb-6.3.3/bin/xtb"
os.chdir(wdir)
# =============================================================================
# #MYRIAD
#wdir = "/home/uccahcl/Scratch/FIND_CAP/ledock_ligand_design_with_SMALL_MOLS_NEW_SCORING_RST_XTB"
#obabel_path = "/home/uccahcl/anaconda2/envs/chemts/bin/obabel"
#ledock_path = "" # Ledock should be in the same folder as wdir
# apbs_path = "/home/uccahcl/apbs-pdb2pqr/bin/apbs"
# antechamber_path = "/home/uccahcl/anaconda2/envs/chemts/bin"
# cp2k_path = "/home/uccahcl/cp2k/exe/local/cp2k.sopt" # double check that cp2k is executing on a single core as it should
# g16_path = "/home/uccahcl/g16/g16"
# #g16_path = "/shared/ucl/apps/gaussian/g16-a03/pgi-2016.5/g16/g16"
#xtb_path = "/home/uccahcl/xtb-6.3.3/bin/xtb"
#os.chdir(wdir)
# os.environ["GAUSS_EXEDIR"]="/home/uccahcl/g16"
# =============================================================================
# -*- coding: utf-8 -*-
def build_mol_from_smiles(SMI=None, pdbfile=None, mol = None, NAME=None):
	"""

	Parameters
	----------
	SMI : TYPE a string representing a molecule, needs to be valid otherwise returns None
		DESCRIPTION.
		 trop is "C1=CC=C[CH+]C=C1"
	Returns a 3D version of the molecule
	-------
	Also produces a pdb file
	"""
        os.chdir(wdir)
	#BUILDS 3D structure
	if SMI is not None:
		mol = Chem.MolFromSmiles(SMI)
		mol = Chem.AddHs(mol)
		AllChem.EmbedMolecule(mol)
		
	elif pdbfile is not None:
		mol = Chem.MolFromPDBFile(pdbfile, removeHs=False)
	elif mol is not None:
		pass 
	else:
		print "No valid input. Provide valid SMI, pdbfile or mol"
		return 
	AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
	#GETS FORMAL CHARGE
	charge = Chem.GetFormalCharge(mol)
	#ASSIGN RANDOM NAME
	if NAME == None:
		rdstring = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(32)])
	else: rdstring = NAME
	#SAVE IN OUTPUTS
	Chem.MolToPDBFile(mol, "outputs/{}.pdb".format(rdstring), flavor=28)	

	shell("{0}  outputs/{1}.pdb -O outputs/{1}.mol2".format(obabel_path, rdstring).split())
	
	return mol, charge, rdstring



def make_docking_script(rdid, docking_targetPDB):
	"""

	Parameters
	----------
	docking_targetPDB : a pdb file that is the docking target for the docking procedure
		DESCRIPTION.

	Returns
	-------
	None.
	"""
        os.chdir(wdir)
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
	return 

def dock_mol_to_host(rdid, pdbtarget):
	"""
	Parameters
	----------
	mol : TYPE the mol id that is assumed to be in wdir/outputs/rdid.mol2
		DESCRIPTION.

	Returns 
	-------
	None.

	"""
        os.chdir(wdir)
	make_docking_script(rdid, pdbtarget)
	with open("outputs/{}-ligands.list".format(rdid), "w") as w:
		w.write("outputs/{}.mol2".format(rdid))
		
	ledock_cmd = "./ledock_linux_x86 outputs/{}-ledock.config".format(rdid).split()
	print ledock_cmd
	proc = subprocess.check_output(ledock_cmd, shell=False)

	return 

def MolToXYZFile(mol, fname):
        """
        PRE: Takes in a rdkitmol
        POST: Willwrite out a xyz file
        """
        conf = mol.GetConformer(-1)
        xyzs = conf.GetPositions()
        at = [a.GetSymbol() for a in mol.GetAtoms()]
        with open(fname, "w") as w:
                w.write("{}\n\n".format(len(at)))
                w.writelines(["{0:<3}{1: 3.8f}    {2: 3.8f}    {3: 3.8f}\n".format(
                a,x[0], x[1], x[2]) for a, x in zip(at,xyzs)])

def dok2pdb(rdid, charge, reconstructing_pdbtarget, n=0, makemol2=False):
	"""
	PRE: TAKES in a rdid with a corresponding mol2 and dok file for the guest
	POST: Will produce the n best complexes that correspond to the docking in pdb and mol2 formats
	"""
        os.chdir(wdir)

	summaryfile = "summary_file"
	with open(summaryfile, 'wb'): pass
	results = []
	
	# ERASE OLD AND SPLITS IN INDIVIDUAL PDB FILES
	for f in glob.glob("outputs/{}-{}.pdb".format(rdid,"*")):
		with open(f, 'wb') as w: pass
	
	i=0
	#SPLITS THE DOK INTO INDIVIDUAL PDBs
	guest_list = []
	#try:
	with open("outputs/{}.dok".format(rdid), 'rb') as r:
		for line in r:
			if i>n: break
			curcomplex = "outputs/{}-{}.pdb".format(rdid, i)
			with open(curcomplex, "ab") as a:
				#a.write(line)
				if "ATOM" in line:
					pts = line.strip().split()
					# EXAMPLE:  ATOM      2  C1  LIG     0      -6.550  -3.810  -2.641
					pt2lower = pts[2][0].upper() + pts[2][1:].lower()
					rec_line = ("{0: <8}{1: >3}  {2: <4}{3}     {4: >1}     {5: >7} {6: >7} {7: >7}\n".format(pts[0], pts[1], pt2lower, pts[3], pts[4], 
											"{0: 2.3f}".format(float(pts[5])), "{0: 2.3f}".format(float(pts[6])), "{0: 2.3f}".format(float(pts[7]))))
					print rec_line
				else:
					rec_line = line
				a.write(rec_line)
			if "END" in line:
				#make_pdb_with_named_residue(curcomplex, "GST", CAPS=True)
				guest_list.append(curcomplex)
				i=i+1
			elif "Cluster" in line:
				with open(summaryfile, 'ab') as a:
					#a.write('{}\t{}\t{}'.format(name, i, line))
					results.append((float(line.split()[-2]), rdid, i))
# =============================================================================
# 	except:
# 		pass
# =============================================================================
	
	# CREATES THE COMPLEXES
	cb = Chem.MolFromPDBFile(reconstructing_pdbtarget, removeHs=False) #HETATM tags rather than ATOM are needed for safe reading
	complex_list = []
	for f in sorted(glob.glob("outputs/{}-{}.pdb".format(rdid, "*"))):
		print(f)
		try:
			guest = Chem.MolFromPDBFile(f, removeHs=False)
			cbguest = Chem.CombineMols(cb, guest)
			complex_pdb_fname = f.replace(rdid, rdid+'-CB')
			guest_xyz_name, complex_xyz_name = f.replace(".pdb", ".xyz"), complex_pdb_fname.replace(".pdb", ".xyz")
			Chem.MolToPDBFile(guest, f)
			Chem.MolToPDBFile(cbguest, complex_pdb_fname)
                        MolToXYZFile(guest, guest_xyz_name) # here creating sdfs for xtb handling
                        MolToXYZFile(cbguest, complex_xyz_name) # creating sdfs for xtb handling
			complex_list.append(complex_pdb_fname)
		except Exception as e:
			print "Complex reassembling failed"
			print str(e)
	with open(summaryfile, 'ab') as a:
		for res in sorted(results):
			a.write("{}\t{}\t{} kcal/mol\n".format(res[1], res[2], res[0]))
	complex_list = [x.split("/")[-1][:-4] for x in complex_list]
	guest_list = [x.split("/")[-1][:-4] for x in guest_list]

	return complex_list[0], guest_list[0]



def get_xtb_eandg(molname, charge, test=False):
	"""
	PRE : There is a molname file that can be read in the output/ directory 
	POST: will perform the xtb optimization there
	"""
	os.chdir(wdir)
	os.chdir("outputs")
	os.environ["KMP_INIT_AT_FORK"]="FALSE"
	os.environ["OMP_NUM_THREADS"]="4"
	os.environ["MKL_NUM_THREADS"]="4"
	os.environ["OMP_STACKSIZE"]="4G"
	command = "{0} {1}.xyz -c {2} -u 0 -P 4 --alpb water --namespace {1} --ohess tight".format(
	xtb_path, molname, charge)
	E,G=0,0
	
	if test:
		command = "{0} {1}.xyz -c {2} -u 0 -P 4 --alpb water --namespace {1} --opt --cycles 20".format(
											xtb_path, molname, charge)
		G = 2.0

	res=shell(command.split())


	for line in res[::-1]:
		if "TOTAL ENERGY" in line:
			E=line.split()[3]
		if "TOTAL FREE ENERGY" in line:
			G=line.split()[4]
	return E, G


def shell(command):
    try:
        output = subprocess.check_output(command, shell=False, stderr=subprocess.STDOUT)
    except Exception as e:
        output = str(e)
    finished = str(output).split('\n')
#    for line in finished:
#        print line 
    return finished

def estimate_dG_xtb(mol_representation, ref_dic={"E_tot":-283.941373569825*627.5, "G_tot":-282.893658634987*627.5}, NAME=None):
	"""
	Parameters
	----------
	SMI : smiles
		smiles as an input for the computation of dG for a cap binding with a binary inclusion complex

	Returns
	-------
	the binding affinity in kcal/mol and a random id string
	the method will
	- Create a molecule and its structure in pdb and mol2
	- dock it using ledock 
	- build the complexes from the docked geometries
	- write the input files for geometry optimisation, frequency calculation and solvent corrected calculation
	at the time of writing, the method is based on PM7 and optimised in internal coordinates for AT MOST 100 stepts using a loose grid
	and loose stopping criteria
	- collect the values for these computations and write out an estimate for the binding affinity of the smiles 
	"""
	# print build_the_reference_dictionary().__repr__()
	t0 = time.time()
	os.chdir(wdir)
	print "BUILDING THE MOL, DOCKING IT AND RECONSTRUCTING THE RESULTING PDBs ({0:4.4f}s) SMI: {1}".format(time.time()-t0, mol_representation)
	mol, charge, rdstring = build_mol_from_smiles(SMI=mol_representation, NAME=NAME) # creates the molecule from smiles
	dock_mol_to_host(rdstring, docking_targetPDB) # will dock the molecule
	best_complex, best_guest = dok2pdb(rdstring, charge, docking_targetPDB) # will convert the dok file to a pdb again
	
	print "OPTIMISE THE RESULTING COMPLEX USING G16 ({0:4.4f}s)".format(time.time()-t0)
	complex_E, complex_G = get_xtb_eandg(best_complex, charge) # produces an rdkit opti file
	os.chdir(wdir)
	shell("{0}  outputs/{1}-CB-0.xtbopt.xyz -O outputs/{1}-CB-0.xtbopt.pdb".format(obabel_path, rdstring).split())

	print "OPTIMISE THE BEST GUEST USING G16 ({0:4.4f}s)".format(time.time()-t0)
	guest_E, guest_G = get_xtb_eandg(best_guest, charge) # produces an rdkit opti file

	print "="*30
	print "CAP"
	guest_E, guest_G = float(guest_E)*627.5, float(guest_G)*627.5
	cap_summary = "E_tot: {0:4.4f}, G_tot: {1:4.4f} (all kcal/mol)".format(guest_E, guest_G)
	print cap_summary

	print "TERNARY COMPLEX"
	complex_E, complex_G = float(complex_E)*627.5, float(complex_G)*627.5
	complex_summary = "E_tot: {0:4.4f}, G_tot: {1:4.4f}".format(complex_E, complex_G)
	print complex_summary

	print "BINARY COMPLEX (target)"
	binary_summary = "E_tot: {0:4.4f}, G_tot: {1:4.4f}".format(
		ref_dic["E_tot"],
		ref_dic["G_tot"]
		)
	print binary_summary
	print "difference"
	dE = complex_E-guest_E-ref_dic["E_tot"]
	dG = complex_G-guest_G-ref_dic["G_tot"]
	if not all([complex_E,guest_E, complex_G,guest_G]):
		dG=10**10
		
	difference_summary = "dE_tot: {0:4.4f}, dG_tot: {1:4.4f}".format(dE, dG)
	print difference_summary
	print "="*30
	#tar_it_all(rdstring)
	os.chdir(wdir)
	return rdstring, dG, {"CAP":cap_summary, "TER_CMP": complex_summary, "diff":difference_summary}

		
#def tar_it_all(rdid):
#	"""
#	rdid : a molecule id value
#	will tar all related files to outputs/rdid-ALL-SAVED.tar
#	"""
#       os.chdir(wdir)
#	with tarfile.open("outputs/{}-ALL-SAVED.tar".format(rdid), "w:tar") as tar:
#		for el in list(set(glob.glob("outputs/{}*".format(rdid)))-set(glob.glob("outputs/*tar"))):
#			tar.add(el)
#			os.remove(el)
#		
#	return 

def tar_it_all(rdid, keep=True):
    """
    rdid : a molecule id value
    will tar all related files to outputs/rdid-ALL-SAVED.tar
    """
    os.chdir(wdir) 
    os.chdir("outputs")
    GauNums = []
    for logs in glob.glob("{}*.log".format(rdid)): # DITCH THE GAU TMP FILES
        with open(logs,"r") as r:
            for line in r:
                if "PID" in line:
                    GauNums.append(line.split()[-1][:-1])
                    break
    for nums in GauNums:
        for Gau in glob.glob("Gau-{}.*".format(nums)):
            os.remove(Gau)

    with tarfile.open("{}-ALL-SAVED.tar".format(rdid), "w:gz") as tar:
        for el in list(set(glob.glob("{}*".format(rdid)))-set(glob.glob("*tar"))):
            tar.add(el)
            os.remove(el)
    if not keep:
        os.remove("{}-ALL-SAVED.tar".format(rdid))
    for core in glob.glob("core.*"): # DITCH THE CORES THAT POP UP 
        os.remove(core)

    os.chdir("..")
    return 

if __name__=="__main__":
	"""
	ACCESS to the RDKIT, ledock, antechamber, apbs and cp2k are required. Access to a docked 
	The workflow goes as follows.
	From a molecule represented as smiles, a 3D configuration is created with the RDKIT
	The 3D molecule is docked using ledock
	antechamber is used to create topology files 
	
	additional lines
	/home/macenrola/anaconda3/envs/chemts/bin/parmchk2 -i $guest.mol2 -f mol2 -o $guest.frcmod -f frcmod

	additional notes, the names of the atoms like "C14" need to be strictly conserved in the pdb, mol2 and complex pdb
	ideally you should have all the pdbs in the correct position, correct labelling and then get the mol2, frcmod, lib, 
	prmtops and finally call all for the complex
	"""


	
	## TEST
	#rdstring, dG, detailed_res_dic = estimate_dG_g16("ClC(Cl)(Cl)Cl", binary_complex_values, "Ctetrachloride")
	#print os.getcwd()
	#
	#os.chdir(wdir)
	#os.chdir("outputs")
	#print os.getcwd()
	#print get_xtb_eandg("CULfF8RkSpXMxlaBdZbo8I3l3sl1Fm2x", 0)

	#dok2pdb("Lhnf1Lh72n6hGaHKxhiJyC3GK24vOs1P", 0, docking_targetPDB)
