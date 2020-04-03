# import our functions I guess we will eventually change read_files name
from read_files import *
from modeller_functions import *
from plots import *
import os
import sys
import copy
import Bio
import numpy as np

class Complex():
	"""Stores a complex that can be build from pairs of structures"""
	def __init__(self, structure = False, comunicator = "print", stop_if_warnings = True):
		self.complex_models = structure
		if not self.complex_models:
			self.complex_models = [] # will contain all posible ways to build a complex.
		self.comunicator = comunicator # to be used instead of print, to redirect the messages if needed.
		self.stop_if_warnings = stop_if_warnings # if True, the program will stop when some minor error
		#comunicator(self.get_pdbs("."))          # as an unreadable pdb and ask the user how to continue.


	def get_complex(self, pdb_folder = ".", stoichiometry = False, directory = "models", correct_chains = False, stop_if_warnings = True):
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
			print("corrected")
			print("corrected")
		for pair_index in range(len(self.pairs)):
			for chain_name in self.pairs[pair_index][0].child_dict:
				try:
					self.chain_dict[chain_name[0]].append(pair_index)
				except: # the first time a chain_name used, except will create the coresponding list
					self.chain_dict[chain_name[0]] = [pair_index,]

		if not os.path.exists(directory):
			os.makedirs(directory)

		print(self.pairs)
		print(self.chain_dict)
		self.bulid_complex_no_dna_strange_thinghs_please(self.pairs[0], stoichiometry, directory)


	def bulid_complex_no_dna_strange_thinghs_please(self, base, stoichiometry, directory, chain_types = "missing", missing_tries = "all", rmsd_threshold=10):
		"""Builds the complex without taking into acount posible dificulties associated with dna
		is a semiprivate function, the users should use get_complex, that ends calling this function"""
		model_number = len(self.complex_models)
		if model_number > 50:
			return # limiting number of models
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
			model_chain_id = [x for x in missing_tries.keys()][0]
			chain_in_model = current_model.child_dict[model_chain_id]
			while len(missing_tries[model_chain_id])>0:
				print(missing_tries)
				pair_id = missing_tries[model_chain_id].pop()
				print("-"*100)
				print(pair_id, self.pairs[pair_id][0].child_dict)
				other_chain_id = [x for x in self.pairs[pair_id][0].child_dict.keys() if x[0] == chain_types[model_chain_id]][0]
				other_chain = self.pairs[pair_id][0].child_dict[other_chain_id]

				#seq_in_model = get_sequence(chain_in_model)
				#other_seq = get_sequence(other_chain)
				#alignment = align_sequences(seq_in_model, other_seq)

				#start = alignment[3]+1
				#end = alignment[4]+1
				print("-"*100, "CHAIN_TYPES:", chain_types)
				print(self.pairs[pair_id][0].child_dict,chain_in_model, other_chain, chain_types[model_chain_id])
				rotated_chain = superimpose_chain(current_model, self.pairs[pair_id][0], chain_in_model, other_chain) #start, end)
				rotated_chain = copy.deepcopy(rotated_chain)# if the original pair rotates is not important but we do not want its chain
				if rotated_chain is None: # names to be chanched, and once the rotated chain added to the model, it should not move
					continue # Not sure why can it be None, but it happens
				rotated_id = copy.deepcopy(rotated_chain.id[0])
				if rotated_chain.id in current_model.child_dict.keys():
					# repeated chain name, need to be changes
					rotated_chain.id = [x for x in "QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm" if x not in current_model.child_dict.keys()][0]

				if structure_clashes(current_model, rotated_chain):
					troble_makers = self.get_clash_responsibles(current_model, rotated_chain)
					if len(troble_makers)==1: # cheking if we are trying to add a chain that already is in the model
						try: # that may happen if such chain is in contact with multiple chains
							if get_rmsd(current_model[troble_makers[0]], rotated_chain) <= rmsd_threshold:
								continue # If the chain is already in the model we do not have to do anything
						except:
							a = get_rmsd(current_model[troble_makers[0]], rotated_chain) <= rmsd_threshold
							pass # try and pass structure are to ignore possible rmsd errors due to diferent number of atoms
							# if they have different number of atoms they are not the same chain

					new_model = copy.deepcopy(current_model) # as there is a conflict(clashes) we are going to generate two models
					new_missing = {x: missing_tries[x] for x in missing_tries if x not in troble_makers}
					for chain in troble_makers:
						new_model.detach_child(chain)


					new_model.add(rotated_chain)
					new_missing[rotated_chain.id]=[x for x in self.chain_dict[rotated_id[0]] if x != model_chain_id]
					new_chain_types = copy.deepcopy(chain_types)
					new_chain_types[rotated_chain.id] = rotated_id
					self.bulid_complex_no_dna_strange_thinghs_please(new_model.get_parent(),stoichiometry, directory, new_chain_types, new_missing, rmsd_threshold = rmsd_threshold)

				else:
					current_model.add(rotated_chain)
					chain_types[rotated_chain.id] = rotated_id
					missing_tries[rotated_chain.id] = [x for x in self.chain_dict[rotated_id[0]] if x != model_chain_id]
				print(self.pairs[pair_id][0].child_dict)
			if len(self.chain_dict)> 100:
				break # maximum chains number


			missing_tries.pop(model_chain_id)
		try:
			write_pdb(current_model, model_number, directory)
		except:
			print(model_number,current_model.child_dict)

	def get_clash_responsibles(self,current_model, rotated_chain):
		"""returns the chains in current_model coliding with the new rotated chain."""
		print("GET CLASH RESPONSABLES, in-chain:",rotated_chain)
		troublemakers = []
		for chain_id in current_model.child_dict:
			print("   other chain_id:",chain_id)
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



	def get_chains_corrected(self,threshold = 1):
		"""Renames the chains so that equal or really similar sequences get the same name."""
		names = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		sequences = []
		original_names = [] # debugging use only
		for pair_id in range(len(self.pairs)):
			#new_child_dict = {}
			old_chain_names = [x for x in self.pairs[pair_id][0].child_dict.keys()]
			print(self.pairs[pair_id][0].child_dict, old_chain_names)
			for chain in old_chain_names:
				seq = get_sequence(self.pairs[pair_id][0][chain])
				print("chain", chain)

				for i in range(len(names)):
					if i>=len(sequences):
						print("i==len")
						sequences.append(seq)
						renamed_chain = self.pairs[pair_id][0][chain]
						original_names.append(copy.copy(renamed_chain.id))
						self.pairs[pair_id][0].detach_child(chain)
						try:
							renamed_chain.id = names[i]
							self.pairs[pair_id][0].add(renamed_chain)
						except:
							renamed_chain.id = names[i]+"1" # if a pair contains 2 equal structures this 1 will help diferenciate them
							self.pairs[pair_id][0].add(renamed_chain)
						break



						#new_child_dict[names[i]] = self.pairs[pair_id].child_dict[names[i]]
					else:
						print("i!=len")
						try:
							alignment = align_sequences(seq, sequences[i])
							final_score = alignment[2]/max(len(seq), len(sequences[i]))
						except:
							continue # if the alignment fails assume they are not the same molecule
						if final_score >= threshold:
							print("score, it's chain",original_names[i])
							renamed_chain = self.pairs[pair_id][0][chain]
							self.pairs[pair_id][0].detach_child(chain)
							try:
								renamed_chain.id = names[i]
								self.pairs[pair_id][0].add(renamed_chain)
							except:
								renamed_chain.id = names[i]+"1"
								self.pairs[pair_id][0].add(renamed_chain)
							break
			#self.pairs[pair_id].child_dict = new_child_dict

	def energy_profiles(self, directory):
		for i in range(len(self.complex_models)):
			get_profile(directory, i)
			plot_profile(directory, i)

	def optimize_models(self, directory):
		for i in range(len(self.complex_models)):
			optimize(directory, i)


def get_rmsd(chain_a,chain_b):
	"""Uses QCPSuperimposer from Biopython to obtain the rmsd of two
	chains without superimposing."""
	Super = Bio.SVDSuperimposer.SVDSuperimposer()
	str_atoms_sup = []
	other_atoms_sup = []

	for residue1, residue2 in zip(chain_a,chain_b):
		try:
			if residue1[('CA')].get_vector() and residue2[('CA')].get_vector():
				str_atoms_sup.append([x for x in residue1[('CA')].get_vector()])
				other_atoms_sup.append([x for x in residue2[('CA')].get_vector()])
		except:
			pass
	#print(len(str_atoms_sup), str_atoms_sup[0])
	if len(str_atoms_sup)==0:
		return(0)
	Super.set(np.array(str_atoms_sup), np.array(other_atoms_sup))
	return(Super.get_init_rms())
