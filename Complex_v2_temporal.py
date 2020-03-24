# import our functions I guess we will eventually change read_files name
from read_files import *
import os
import sys
import copy
import Bio

class Complex():
	"""Stores a complex that can be build from pairs of structures.
	functions contained:"""
	def __init__(self, structure = False, comunicator = print, stop_if_warnings = True):
		"""  """
		self.complex_models = structure
		if not self.complex_models:
			self.complex_models = [] # will contain all posible ways to build a complex.
		self.comunicator = comunicator # to be used instead of print, to redirect the messages if needed.
		self.stop_if_warnings = stop_if_warnings # if True, the program will stop when some minor error
		comunicator(self.get_pdbs("."))          # as an unreadable pdb and ask the user how to continue.


	def get_complex(self, pdb_folder = ".", stoichiometry = False, correct_chains = True, stop_if_warnings = True):
		"""Builds the complex using pdbs in pdb_folder wich should contain
		pairs of interacting structures. Stoichiometry can be given as a
		dictionary with chain names as key, if not given, all fitting structures
		will be included with an upper limit of 50 repeats per structure.
		If the chain names in the files are not meaningful use correct_chains = False"""

		pdbs = self.get_pdbs(pdb_folder)
		self.pairs = []
		self.chain_dict = {} # stores in what pairs can each chain be found
		for pdb in pdbs:
			try:
				self.pairs.append(read_pdb(pdb))
			except:
				self.comunicator("WARNING: Unable to read file "+pdb)
				if stop_if_warnings:
					u=input("""How do you want to proced?
					write "exit" to stop the program
					write "continue" to continue without the file
					write a path to a pdb folder to substitute the one giving errors""")
					if u=="exit":
						print("See u later") # have to check if this really stops the execution
						sys.exit(0)
					if u!="continue":
						pdbs.append(u) # the new pdb will be tried to open at the end
					# if the user wanted to continue, nothing have to be done

		if not correct_chains:
			self.get_chains_corrected()
		for pair_index in range(len(self.pairs)):
			for chain_name in self.pairs[pair_index][0].child_dict:
				try:
					self.chain_dict[chain_name].append(pair_index)
				except: # the first time a chain_name used, except will create the coresponding list
					self.chain_dict[chain_name] = [pair_index,]

		self.bulid_complex_no_dna_strange_thinghs_please(self.pairs[0], stoichiometry)


	def bulid_complex_no_dna_strange_thinghs_please(self, base, stoichiometry, chain_types = "missing", missing_tries = "all"):
		"""Builds the complex without taking into acount posible dificulties associated with dna
		is a semiprivate function, the users should use get_complex, that ends calling this function"""
		model_number = len(self.complex_models)
		self.complex_models.append(copy.deepcopy(base)) # we have to copy it as we want to have several models
		# if we pass it by reference we would end modifing all models the same time. A more memory eficient
		# method would be to store only the transformation matrices for the models but we do not expect memory
		# problems
		miss = False
		current_model = self.complex_models[model_number][0]

		if chain_types == "missing":
			miss = True
			chain_types ={}
		if missing_tries=="all":
			missing_tries = {}
			for chain_name in current_model.child_dict:
				missing_tries[chain_name] = self.chain_dict[chain_name]
				if miss:
					chain_types[chain_name]=chain_name # As some chains will be repeated, we will have to rename them. This dictionary will help to
					# keep track of what type it is each chain.
		if not stoichiometry:
			stoichiometry = {x:50 for x in self.chain_dict}

		while len(missing_tries)>0:
			for model_chain_id in missing_tries:
				chain_in_model = current_model.child_dict[model_chain_id]
				for pair_id in missing_tries[model_chain_id]:
					other_chain_id = [x for x in self.pairs[pair_id][0].child_dict.keys() if x != chain_types[model_chain_id]][0]
					print("-"*100)
					print(pair_id)
					other_chain = self.pairs[pair_id][0].child_dict[other_chain_id]

					#seq_in_model = get_sequence(chain_in_model)
					#other_seq = get_sequence(other_chain)
					#alignment = align_sequences(seq_in_model, other_seq)

					#start = alignment[3]+1
					#end = alignment[4]+1

					rotated_chain = superimpose_chain(current_model, self.pairs[pair_id][0], chain_in_model, other_chain) #start, end)
					if rotated_chain.id in current_model.child_dict.keys():
						# repeated chain name, need to be changes
						rotated_chain.id = [x for x in "QWERTYUIOPASDFGHJKLZXCVBNM" if x not in current_model.child_dict.keys()][0]

					if structure_clashes(current_model, rotated_chain):
						troble_makers = self.get_clash_responsibles(current_model, rotated_chain)
						if len(troble_makers)==1: # cheking if we are trying to add a chain that already is in the model
							try: # that may happen if such chain is in contact with multiple chains
								if get_rmsd(current_model[troble_makers[0]], rotated_chain) <= rmsd_threshold:
									continue # If the chain is already in the model we do not have to do anything
							except:
								pass # try and pass structure are to ignore possible rmsd errors due to diferent number of atoms
								# if they have different number of atoms they are not the same chain

						new_model = copy.deepcopy(current_model) # as there is a conflict(clashes) we are going to generate two models
						new_missing = {x: missing_tries[x] for x in missing_tries if x not in troble_makers}
						for chain in troble_makers:
							new_model.detach_child(chain)


						new_model.add(rotated_chain)
						new_missing[rotated_chain.id]=self.chain_dict[other_chain_id]
						new_chain_types = copy.deepcopy(chain_types)
						new_chain_types[rotated_chain.id] = other_chain_id
						self.bulid_complex_no_dna_strange_thinghs_please(new_model.get_parent(),stoichiometry, new_chain_types, new_missing)

					else:
						current_model.add(rotated_chain)
						chain_types[rotated_chain.id] = other_chain_id
						missing_tries[rotated_chain.id]=self.chain_dict[other_chain_id]


				missing_tries.pop(chain_missing)

		write_pdb(current_model)

		"""
		PSEUDOCODE -> missing functions
		superimpose (chain_in_model, other_chain)
		rotate and translate the other chain, add it to the current model
		Check for clashes.
		if clashes:
			search for the responsables of such clashes
			check if the clashes are becouse we are trying to add a chain that already was there -> rmsd
				if so, remove the chain from the model and continue, this is good.
				if not so:
					recursively call this function with the current model and the chash responsivles removed (without removing the new chain)
					remove the new chain and continue
		"""
	def get_clash_responsibles(self,current_model, rotated_chain):
		"""returns the chains in current_model coliding with the new rotated chain."""
		troublemakers = []
		for chain_id in current_model.child_dict:
			if structure_clashes(current_model[chain_id], rotated_chain):
				troublemakers.append(chain_id)
		return troublemakers

	def get_pdbs(self, folder):
		"""Gets all pdbs in the given folder or its subfolders."""
		files_and_subfolders = os.listdir(folder)
		pdb_files = []
		for f_or_d in files_and_subfolders:
			full_f_or_d = os.path.join(folder,f_or_d)

			if not os.path.isfile(full_f_or_d):
				try:
					# found some broken folders in my labtop and writen this to overcome them,
					# probably not needed for most users.
					pdb_files += self.get_pdbs(full_f_or_d)
				except:
					self.comunicator("Warning: tried to acces subdirectory " +
					 full_f_or_d + "and failed, files in such subdirectory will not be used.")
					if self.stop_if_warnings:
						u=input("Do you want to stop the execution?(y/n)")
						#print(u)
						if u =="y":
							print("See u later") # have to check if this really stops the execution
							sys.exit(0)
			else:
				if full_f_or_d.split(".")[-1]=="pdb":
					pdb_files.append(full_f_or_d)
		return pdb_files



	def get_chains_corrected(threshold = 1):
		"""Renames the chains so that equal or really similar sequences get the same name."""
		names = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		sequences = []
		for pair_id in range(len(self.pairs)):
			new_child_dict = {}
			for chain in self.pairs[pair_id].child_dict:
				seq = get_sequence(self.pairs[pair_id].child_dict[chain])

				for i in range(len(names)):
					if i==len(sequences):
						sequences.append(seq)
						new_child_dict[names[i]] = self.pairs[pair_id].child_dict[chain]
					else:
						alignment = align_sequences(seq, sequences[i])
						final_score = alignment[2]/max(len(seq), len(sequences[i]))
						if final_score >= threshold:
							new_child_dict[names[i]] = self.pairs[pair_id].child_dict[chain]

			self.pairs[pair_id].child_dict = new_child_dict


def get_rmsd(chain_a,chain_b):
	"""Uses QCPSuperimposer from Biopython to obtain the rmsd of two
	chains without superimposing."""
	Super = Bio.PDB.QCPSuperimposer()
	Super.set(chain_a, chain_b)
	return(Super.get_init_rms())
