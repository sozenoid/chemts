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
obabel_path = "/home/macenrola/anaconda2/envs/obabel/bin/obabel"
ledock_path = "" # Ledock should be in the same folder as wdir
apbs_path = "/usr/bin/apbs"
antechamber_path = "/home/macenrola/anaconda3/envs/chemts/bin"
cp2k_path = "/home/macenrola/anaconda3/pkgs/cp2k-6.1.0-hc6cf775_3/bin/cp2k.sopt" # double check that cp2k is executing on a single core as it should
g16_path = "/home/macenrola/Documents/Gaussian/g16/g16"
xrb_path = "/home/macenrola/xtb-6.3.3/bin/xtb"
os.chdir(wdir)
# =============================================================================
# #MYRIAD
# wdir = "/home/uccahcl/Scratch/FIND_CAP/ledock_ligand_design_with_SMALL_MOLS_NEW_SCORING_RST"
# obabel_path = "/home/uccahcl/anaconda2/envs/chemts/bin/obabel"
# ledock_path = "" # Ledock should be in the same folder as wdir
# apbs_path = "/home/uccahcl/apbs-pdb2pqr/bin/apbs"
# antechamber_path = "/home/uccahcl/anaconda2/envs/chemts/bin"
# cp2k_path = "/home/uccahcl/cp2k/exe/local/cp2k.sopt" # double check that cp2k is executing on a single core as it should
# g16_path = "/home/uccahcl/g16/g16"
# #g16_path = "/shared/ucl/apps/gaussian/g16-a03/pgi-2016.5/g16/g16"
# os.chdir(wdir)
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
					rec_line = ("{0: <8}{1: >3}  {2: <4}{3}     {4: >1}      {5: >6}  {6: >6}  {7: >6}\n".format(pts[0], pts[1], pt2lower, pts[3], pts[4], pts[5], pts[6], pts[7]))
					#print rec_line
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
			Chem.MolToPDBFile(cbguest, complex_pdb_fname)
			complex_list.append(complex_pdb_fname)
		except:
			print "Complex reassembling failed"
	with open(summaryfile, 'ab') as a:
		for res in sorted(results):
			a.write("{}\t{}\t{} kcal/mol\n".format(res[1], res[2], res[0]))
	complex_list = [x.split("/")[-1][:-4] for x in complex_list]
	guest_list = [x.split("/")[-1][:-4] for x in guest_list]

	return complex_list[0], guest_list[0]


def get_g16_optimised_mol(rdid):
	"""
	PRE  : a molecule outputs/ZINC000003874807-opti-G16.log exists an contains an
	POST : Will return a file PDB containing the optimised version
	"""
        os.chdir(wdir)

	opti_traj = "outputs/{}-opti-G16.log".format(rdid)
	opti_best = "outputs/{}-OPTI-G16.pdb".format(rdid)

	convert_to_pdb_cmd = "{} {} -O {}".format(obabel_path, opti_traj, opti_best)
	shell(convert_to_pdb_cmd.split())
	return opti_best.split("/")[-1][:-4]

	
def run_g16_file(g16_input_file):
	"""
	Parameters
	----------
	g16_input_file : string
		DESCRIPTION. The name should contained the labels opti, solv or freq, ONCE 
		otherwise it gets confusing

	Returns None: will just run the g16_input_file and return the energy value for opt and the vib value for VIBRATIONAL_ANALYSIS
	Note that for the energies, here we specifically look for restricted HF and PM7
	-------
	None.
	"""
	#scratch_dir = "{}/{}-GAUSS_SCRATCH".format("$TMPDIR", g16_input_file)
	#run_cmd = "module load gaussian/g16-a03/pgi-2016.5 & mkdir -p $GAUSS_SCRDIR & g16 {1}".format(g16_path, g16_input_file).split()
        os.chdir(wdir)
	run_cmd = "{0} {1}".format(g16_path, g16_input_file).split()
	os.chdir("outputs")

	print run_cmd
	shell(run_cmd)

	outf = g16_input_file[:-4]+".log"
	with open(outf, 'r') as r:
		result = r.readlines()
	if "opti" in g16_input_file:
		for l in result[::-1]:
			"""
			Looking for the last occurence of this block 
			 SCF Done:  E(RPM7) = -0.961086179292     A.U. after   15 cycles
			"""
			if "SCF Done:  E(RPM7) =" in l:
				energy = float(l.split()[-5])
				os.chdir("..")
				return energy*627.5 #returns kcal/mol
	elif "freq" in g16_input_file:
		""" Looking for this block of text 
			            E (Thermal)             CV                S
                      KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
 Total                  854.563            325.764            406.001
 Electronic               0.000              0.000              0.000
 Translational            0.889              2.981             47.615
 Rotational               0.889              2.981             41.800
 Vibrational            852.786            319.803            316.586
		"""
		for l in result[::-1]:
			if " Total              " in l:
				s_tot = float(l.strip().split()[-1])
				os.chdir("..")
				return s_tot # Returns the total entropy as cal/mol
	if "solv" in g16_input_file:
		for l in result[::-1]:
			"""
			Looking for the last occurence of this block 
			 SCF Done:  E(RPM7) = -0.961086179292     A.U. after   15 cycles
			"""
			if "SCF Done:  E(RPM7) =" in l:
				energy = float(l.split()[-5])
				os.chdir("..")
				return energy*627.5 #returns kcal/mol
	else:
		os.chdir("..")
		return

def shell(command):
    try:
        output = subprocess.check_output(command, shell=False, stderr=subprocess.STDOUT)
    except Exception as e:
        output = str(e)
    finished = output.split('\n')
    for line in finished:
        print line
    time.sleep(5)
    return

def estimate_dG_g16(mol_representation, ref_dic={"E_tot":-0.6466416447*627.5, "S_tot": 367.730, "E_CPCM":-0.790020704373*627.5}, NAME=None):
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
	complex_opt_g16_file = make_g16_input_file(best_complex, charge, "opti") # produces an rdkit opti file
	complex_energy = run_g16_file(complex_opt_g16_file)
	complex_opti_pdb = get_g16_optimised_mol(best_complex)

	print "OPTIMISE THE BEST GUEST USING G16 ({0:4.4f}s)".format(time.time()-t0)
	guest_opt_g16_file = make_g16_input_file(best_guest, charge, "opti") # produces an rdkit opti file
	guest_energy = run_g16_file(guest_opt_g16_file)
	guest_opti_pdb = get_g16_optimised_mol(best_guest)
	
 	print "PERFORM THE NMA FOR THE COMPLEX AND EXTRACT ENTROPY VALUES ({0:4.4f}s)".format(time.time()-t0)
 	complex_vib_g16_file = make_g16_input_file(complex_opti_pdb, charge, "freq") # produces an rdkit vibrational file
 	S_tot_complex = run_g16_file(complex_vib_g16_file)
	 
  	print "PERFORM THE NMA FOR THE GUEST AND EXTRACT ENTROPY VALUES ({0:4.4f}s)".format(time.time()-t0)
 	guest_vib_g16_file = make_g16_input_file(guest_opti_pdb, charge, "freq") # produces an rdkit vibrational file
 	S_tot_guest = run_g16_file(guest_vib_g16_file)
	 
	 
	print "OBTAINS SRCF ENERGY FOR SOLVATION OF THE COMPLEX ({0:4.4f})".format(time.time()-t0)
	complex_solv_g16_file = make_g16_input_file(complex_opti_pdb, charge, "solv")
	E_CPCM_complex = run_g16_file(complex_solv_g16_file)
	
	print "OBTAINS SRCF ENERGY FOR SOLVATION OF THE GUEST ({0:4.4f})".format(time.time()-t0)
	guest_solv_g16_file = make_g16_input_file(guest_opti_pdb, charge, "solv")
	E_CPCM_guest = run_g16_file(guest_solv_g16_file)
	
	print "="*30
	print "CAP"
	cap_summary = "E_tot: {0:4.4f}, S_tot: {1:4.4f}, E_CPCM: {2:4.4f} (all kcal/mol except S_tot cal/mol/K)".format(
		guest_energy,
		S_tot_guest,
		E_CPCM_guest
		)
	print cap_summary
	print "TERNARY COMPLEX"
	complex_summary = "E_tot: {0:4.4f}, S_tot: {1:4.4f}, E_CPCM: {2:4.4f}".format(
		complex_energy,
		S_tot_complex,
		E_CPCM_complex
		)
	print complex_summary
	print "BINARY COMPLEX (target)"
	binary_summary = "E_tot: {0:4.4f}, S_tot: {1:4.4f}, E_CPCM: {2:4.4f}".format(
		ref_dic["E_tot"],
		ref_dic["S_tot"],
		ref_dic["E_CPCM"]
		)
	print binary_summary
	print "difference"
	dG = E_CPCM_complex-E_CPCM_guest-ref_dic["E_CPCM"]  - 298.0/1000*(S_tot_complex-S_tot_guest-ref_dic["S_tot"])
	difference_summary = "dE_tot: {0:4.4f}, dS_tot: {1:4.4f}, dE_CPCM: {2:4.4f} |||| dG= {3:4.4f} kcal/mol".format(
		complex_energy-guest_energy-ref_dic["E_tot"],
		S_tot_complex-S_tot_guest-ref_dic["S_tot"],
		(E_CPCM_complex-E_CPCM_guest-ref_dic["E_CPCM"]),
		dG
		)
	print difference_summary
	print "="*30
	tar_it_all(rdstring)
	return rdstring, dG, {"CAP":cap_summary, "TER_CMP": complex_summary, "diff":difference_summary}


def make_g16_input_file(rdid, charge, way):
	"""
	Returns
	-------
	there should be a file like outputs/rdid.pdb or outputs/rdid-CB-0.pdb
	inputs files outputs/rdid.com to run a Gaussian computation
	"""
        os.chdir(wdir)
	mem = 16
	proc = 1
	maxdisk = 30
	if way == "opti":
		method = "#n PM7 maxdisk={}GB Int(Grid=SG1) Opt=(maxcycles=100, loose) SCF=(YQC, maxcycle=1000)".format(maxdisk)  # cartesian to avoid crash with big rotations
	if way == "freq":
		method = "#n PM7 maxdisk={}GB Int(Grid=SG1) freq SCF=(YQC, maxcycle=1000)".format(maxdisk) # cartesian to avoid crash with big rotations
	if way == "solv":
		method = "#n PM7 maxdisk={}GB Int(Grid=SG1) SCRF=CPCM SP SCF=(YQC, maxcycle=1000)".format(maxdisk)
	route = """%NProcShared={3}
%Chk={0}.chk
%mem={4}gb
{2}

{0}

{1} 1
""".format(rdid, charge, method, proc, mem, maxdisk)

	mol = Chem.MolFromPDBFile("outputs/{}.pdb".format(rdid), removeHs=False)
	atoms = [a.GetSymbol() for a in mol.GetAtoms()]
	postring = ["{0:<3}    {1:6> 2.4f}        {2:6> 2.4f}        {3:6> 2.4f}\n".format(a,x[0], x[1], x[2]) for a, x in zip(atoms, mol.GetConformer(-1).GetPositions())]
	outf = "{}-{}-G16.com".format(rdid, way)
	with open("outputs/"+outf, "w") as w:
		w.write(route)
		w.write(''.join(postring))
		w.write("\n\n")
	
	return outf
		
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

def tar_it_all(rdid, keep=False):
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
	ideally you should have all the pdbs in the correct position, correct labelling and then get the mol2, frcmod, lib, prmtops and finally call all for the complex
	"""

	
	## TEST
	#rdstring, dG, detailed_res_dic = estimate_dG_g16("ClC(Cl)(Cl)Cl", binary_complex_values, "Ctetrachloride")
