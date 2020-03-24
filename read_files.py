from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from Bio.PDB.PDBIO import PDBIO
import itertools

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
    structure_matrix = MatrixInfo.blosum62
    alignment = pairwise2.align.localds(seq, other, structure_matrix, -8, -1 )
    return alignment[0]


def compare_sequences(structure, other):

    similarity = {}

    chain_pairs = [ (str_chain, other_chain) for str_chain in structure for other_chain in other]

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

    str_atoms_sup = []
    other_atoms_sup = []

    for residue in chain:
        try:
            str_atoms_sup.append(residue[('CA')])
            other_atoms_sup.append(residue[('CA')])
        except:
            pass

    superimposer.set_atoms(str_atoms_sup, other_atoms_sup)
    superimposer.apply(other.get_atoms())

    for chain in other:
        if chain != other_chain:
            chain.detach_parent()
            return chain


def complex_builder(filelist):

    structure_list = list(map(read_pdb, filelist))
    structure_combinations = list(itertools.combinations(structure_list, 2))

    complex = Structure.Structure("complex")

    for structure, other in structure_combinations:

        similarity = compare_sequences(structure, other)

        other_chain = similarity[1]

        if similarity:
            new_chain = superimpose_chain(structure, other, similarity)

            complex.add(structure)

            complex[0].add(new_chain)

    return complex


def structure_clashes(complex, chain):

    complex_atoms = [atom for atom in list(complex.get_atoms()) if atom.id == "CA"]
    neighbor_search = NeighborSearch(complex_atoms)

    chain_atoms = [atom for atom in list(chain.get_atoms()) if atom.id == "CA"]

    for chain_atom in chain_atoms:
        for complex_atom in neighbor_search.search(chain_atom.get_coord(), 1.2, 'A'):
            clashing_chain = complex_atom.get_parent().get_parent().id
            return clashing_chain

    return False

def write_pdb(structure, number):
    io = PDBIO()
    io.set_structure(structure)
    io.save("model_" + str(number) +".pdb")
    return "File was correctly saved"



if __name__ == "__main__":
    """
    parser = argparse.ArgumentParser(description="Take input and output files...")

    parser.add_argument("-i", "--input", dest="infile", action="store", default= os.getcwd(), help="Input FASTA formatted file")

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Print log in stderr")

    options = parser.parse_args()

    if os.path.isdir(options.infile):
        filelist = [f for f in os.listdir(options.infile) if f.endswith(".pdb")]
    """

    str1 = read_pdb("WyQ.pdb")

    str2 = read_pdb("QyH.pdb")

    str3 = read_pdb("otherH.pdb")

    seq1 = get_sequence(str1)

    seq2 = get_sequence(str2)

    seq3 = get_sequence(str3)

    sim = compare_sequences(str1, str2)

    filelist = ["WyQ.pdb", "QyH.pdb"]

    """

    chainQ1 = str1[0]["Q"]

    chainQ3 = str2[0]["Q"]

    seqQ1 = get_sequence(chainQ1)

    seqQ3 = get_sequence(chainQ3)

    align = align_sequences(seqQ1, seqQ3)
    """
