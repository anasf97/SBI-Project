from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from Bio.PDB.PDBIO import PDBIO
import itertools
import os

def read_pdb(file):
    parser = PDBParser()
    filext = file.split(".pdb")
    pdb_id = filext[0]
    structure = parser.get_structure(pdb_id, file)
    return structure

def get_sequence(structure):
    ppb = PPBuilder()

    sequence = ""

    for polypeptide in ppb.build_peptides(structure):
        sequence += polypeptide.get_sequence()
    return sequence

def align_sequences(seq, other):
    blosum_matrix = MatrixInfo.blosum62
    alignment = pairwise2.align.localds(seq, other, blosum_matrix, -8, -1 )
    return alignment[0]


def compare_sequences(structure, other):

    similarity = {}

    chain_pairs = [ (str_chain, other_chain) for str_chain in structure[0] for other_chain in other[0]]

    for str_chain, other_chain in chain_pairs:

        str_chain_seq = get_sequence(str_chain)
        other_chain_seq = get_sequence(other_chain)

        alignment = align_sequences(str_chain_seq, other_chain_seq)

        final_score = alignment[2]/max(len(str_chain_seq), len(other_chain_seq))

        if final_score >= 1:
            start = alignment[3]
            end = alignment[4]
            similarity = ((str_chain, other_chain), start, end)
            return similarity

    return False

def superimpose_chain(structure, other, chain, other_chain): #(start, end):
    superimposer = Superimposer()

    CA_str_atoms_sup = [residue['CA'] for residue in chain.get_residues() if 'CA' in [atom.get_id() for atom in residue.get_atoms()]]
    CA_other_atoms_sup = [residue['CA'] for residue in other_chain.get_residues() if 'CA' in [atom.get_id() for atom in residue.get_atoms()]]

    P_str_atoms_sup = [residue['P'] for residue in chain.get_residues() if 'P' in [atom.get_id() for atom in residue.get_atoms()]]
    P_other_atoms_sup = [residue['P'] for residue in other_chain.get_residues() if 'P' in [atom.get_id() for atom in residue.get_atoms()]]

    if CA_str_atoms_sup:
        superimposer.set_atoms(CA_str_atoms_sup, CA_other_atoms_sup)
    else:
        superimposer.set_atoms(P_str_atoms_sup, P_other_atoms_sup)

    superimposer.apply(other.get_atoms())

    for chain in other:
        if chain != other_chain:
            chain.detach_parent()
            return chain

def complex_builder(filelist):

    structure_list = list(map(read_pdb, filelist))
    structure_combinations = list(itertools.combinations(structure_list, 2))

    complex = Structure.Structure("complex")

    start = False

    for structure, other in structure_combinations:

        similarity = compare_sequences(structure, other)

        if similarity:

            chain = similarity[0][0]

            other_chain = similarity[0][1]

            if not start:
                complex.add(structure[0])

            new_chain = superimpose_chain(structure[0], other[0], chain, other_chain)

            if not structure_clashes(complex, new_chain):

                complex[0].add(new_chain)

        start = True

    return complex


def structure_clashes(complex, chain):

    complex_atoms = [atom for atom in list(complex.get_atoms())]# if atom.id == "CA"]
    neighbor_search = NeighborSearch(complex_atoms)

    for chain_atom in list(chain.get_atoms()):
        for complex_atom in neighbor_search.search(chain_atom.get_coord(), 1.2, 'A'):
            clashing_chain = complex_atom.get_parent().get_parent().id
            return clashing_chain

    return False

def write_pdb(structure, number, directory):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(directory) + "/model_" + str(number) + ".pdb")




