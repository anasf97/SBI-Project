from Bio.PDB import *
from Bio import pairwise2
from Bio import Align
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SubsMat import MatrixInfo

def read_pdb(file):
    parser = PDBParser()
    filext = file.split(".")
    pdb_id = filext[0]
    structure = parser.get_structure(pdb_id , file)
    return structure

def get_sequence(structure):
    ppb = PPBuilder()

    for polypeptide in ppb.build_peptides(structure):
        sequence = polypeptide.get_sequence()
    return sequence

def align_sequences(seq, other):
    structure_matrix = MatrixInfo.structure
    alignment = pairwise2.align.localds(seq, other, structure_matrix, -4, -1 )
    #alignment = pairwise2.align.localxx(seq, other)
    return alignment[0]


def compare_structures(structure, other):
    str_model = structure[0]
    other_model = other[0]

    similarity = {}

    chain_pairs = [ (str_chain, other_chain) for str_chain in str_model for other_chain in other_model]

    for str_chain, other_chain in chain_pairs:

        str_chain_seq = get_sequence(str_chain)

        other_chain_seq = get_sequence(other_chain)

        alignment = align_sequences(str_chain_seq, other_chain_seq)

        final_score = alignment[2]/max(len(str_chain_seq), len(other_chain_seq))

        if final_score >= 1:

            start = alignment[3]

            end = alignment[4]

            similarity[(str_chain, other_chain)] = (start, end)

    return similarity

def superimpose(structure, other, similarity):
    superimposer = Superimposer()

    str_model = structure[0]
    other_model = other[0]

    str_atoms = []
    other_atoms = []

    for chain_pair in similarity[1]:

        for str_res in str_model[chain_pair[0]]:

            str_atoms.append(str_res[('CA' or "P")])

        for other_res in other_model[chain_pair[1]]:

            other_atoms.append(str_res[('CA' or "P")])

    superimposer.set_atoms(str_atoms, other_atoms)

    superimposer.apply(other_model.get_atoms())

    for chain in other_model:
        for res in chain:
            if res.values() not in other_atoms:
                other_atoms_not_superimp.append(res['CA'])



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

    chainQ1 = str1[0]["W"]

    chainQ3 = str2[0]["H"]

    seqQ1 = get_sequence(chainQ1)

    seqQ3 = get_sequence(chainQ3)

    align = align_sequences(seqQ1, seqQ3)
