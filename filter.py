def keep_common_atoms(fname, atom_list):
	"""
	Parameters
	----------
	fname : string
		list of smiles
	atom_list : list
		list of atoms that are to be kept
	Returns
	-------
	None.
	"""
	atom_numbers = [x[0] for x in atom_list]
	atoms_ref_mass = [x[1] for x in atom_list]
	#s = dict()
	suppl = Chem.SmilesMolSupplier(fname, delimiter='\t')
	w = Chem.SmilesWriter(fname+'_purifed', delimiter="\t")
	for mol in suppl:
		try:
			atoms_nbs = [a.GetAtomicNum() for a in mol.GetAtoms()]
			atoms_mass = [a.GetMass() for a in mol.GetAtoms()]
			#print atoms_nbs, atoms_mass

			if all([a in atom_numbers for a in atoms_nbs]):
				#for k in atoms_mass: 
				#	if k in s: s[k]+=1
				#	else: s[k]=1
				if all([a in atoms_ref_mass for a in atoms_mass]):
					w.write(mol)
			else:
				print "excluding mol:\t{}\t{}".format(Chem.MolToSmiles(mol), mol.GetProp("_Name"))
				pass
		except:
			pass
	#for k in sorted(s):
	#	print s, s[k]  
	w.close()
	return 


def common_atom_list():
	return [(1, 1.008), (6, 12.011), (7, 14.007), (8, 15.999), (9, 18.998), (15, 30.974), (16, 32.067), (17, 35.453), (34, 78.96),
		  (35, 79.904), (53, 126.904)]

def remove_irrelevant_entries(fname, irrelevant_entries):
	with open(fname+"_no_irrelevant", "wb") as w:
		with open(fname, "rb") as r:
			for line in r:
				smi = line.split()[0]
				if is_good(smi, irrelevant_entries):
					w.write(line)
				else:
					print "BAD LINE: {}".format(line)
			



def is_good(smiles, irrelevant_entries):
	for i in irrelevant_entries:
		if i in smiles:
			return False
	return True

def irrelevant_entries():
	all_of_them = ['\n', '&', 'S', 'M', 'I', 'L', 'E', 'C', '(', 'O', ')', 'N', '=', 'P', 'n', '1', 'c', '2', 'Cl', '[O-]', '#', '[PH]', '[N+]', 
		 '[nH]', '[NH4+]', '[Br-]', 'Br', '[C-]', '[O+]', '[Cl-]', 'o', '[n+]', '[S+]', '[SeH2]', '[H]', '[SeH]', '[CH3-]', '/', '[NH3+]', 
		 '[OH-]', '[H+]', '[Se]', 's', '-', '[S-]', 'F', '3', '\\', '[NH-]', '[I-]', '[n-]', '[CH-]', '[NH2+]', '[C@@H]', '[C@H]', '[CH2-]', 
		 '[c-]', '[cH-]', '[nH+]', '[N-]', '[P-]', '[C@]', '[P+]', '[NH+]', '[C]', '[SH-]', '[SH2-2]', 'p', '[Cl+3]', '[se]', '[Cl+2]', '[I+2]', 
		 '[Se+]', '[Br+2]', '[Cl+]', '[NH2-]', '[F-]', '[CH+]', '[CH]', '[O]', '[I+]', '[S-2]', '[N]', '[PH3-3]', '[SeH2-2]', '[PH2-]', '[I]', 
		 '[IH]', '[I+3]', '[C@@]', '[Se-]', '[OH2+]', '[Se+4]', '[Se+6]', '[Se-2]', '[c]', '[CH2]', '[S]', '[OH+]', '[NH2]', '[OH3+]', '4', '[C+]', 
		 '[o+]', '5', '[NH]', '[PH2]', '[SH2+]', '[P]', '[SH]', '[cH+]', '[OH]', '[O-2]', '[Br+]', '[PH4+]', '[H-]', '[Br+3]', '[N@+]', '[s+]', '6', 
		 '[SH3+]', '[SH+]', '[CH2+]', '[pH]', '[CH3+]', '[p-]', '[PH+]', '[CH3]', '[PH-]', '[ClH2+]', '[sH+]', '[SeH2+]', '[IH2+]', '[P-3]', '[SeH-]',
		  '[S@@]', '[Cl]', '[F]', '[Br]', '[SH3]', '[SeH3+]', '[SeH3]', '[PH3+]', '[PH-2]', '[PH2+]', '[SeH+]', '[IH+]', '[SeH2-]', '[ClH+]', '[S@]',
		  '[P@@]', '[P@]', '[c+]', '[se+]', '[p+]', '[P-2]', '[P@@H]', '[P@H]', '[SeH2+4]', '[N@@+]', '[oH+]', '[SeH2+6]']
	
	good_ones = ['\n', '&', 'C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', 'S', 'Cl', '[O-]', '[C@H]', '[NH+]', '[C@]', 's', 'Br',
			   '/', '[nH]', '[NH3+]', '4', '[NH2+]', '[C@@]', '[N+]', '[nH+]', '\\', '[S@]', '5', '[N-]', '[n+]', '[S@@]', '[S-]', '6', '7', 'I', '[n-]', 'P', '[OH+]', 
			    '[NH-]', '[P@@H]', '[P@@]', '[PH2]', '[P@]', '[P+]', '[S+]', '[o+]', '[CH2-]', '[CH-]', '[SH+]', '[O+]', '[s+]', '[PH+]', '[PH]', '8', '[S@@+]']
	
	
	return [x for x in all_of_them  if x not in good_ones ]

if __name__ == "__main__":
	import sys
	import rdkit
	from rdkit import Chem
	
	#keep_common_atoms(sys.argv[1], common_atom_list())
	remove_irrelevant_entries(sys.argv[1], irrelevant_entries())
	