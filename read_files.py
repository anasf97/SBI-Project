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
    structure_matrix = MatrixInfo.blosum62
    alignment = pairwise2.align.localds(seq, other, structure_matrix, -8, -1 )
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
            similarity = ((str_chain, other_chain), start, end)

    return similarity

def superimpose(structure, other, similarity):
    superimposer = Superimposer()

    str_model = structure[0]
    other_model = other[0]

    start = similarity[1]
    end = similarity[2]

    str_atoms_sup = []
    other_atoms_sup = []

    str_chain_sup = similarity[0][0]
    other_chain_sup = similarity[0][1]

    for i in range(start+1, end+1):

        str_atoms_sup.append(str_chain_sup[i][('CA')])
        other_atoms_sup.append(other_chain_sup[i][('CA')])

    superimposer.set_atoms(str_atoms_sup, other_atoms_sup)
    superimposer.apply(other_model.get_atoms())

    other_atoms_not_sup = []
    all_str_atoms = []

    for chain in other_model:
        if chain != other_chain_sup:
            for res in chain:
                other_atoms_not_sup.append(res['CA'])

    for atom in str_model.get_atoms():
        all_str_atoms.append(res['CA'])

    final_structure = other_atoms_not_sup + all_str_atoms

    return final_structure


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

    chainQ1 = str1[0]["Q"]

    chainQ3 = str2[0]["Q"]

    seqQ1 = get_sequence(chainQ1)

    seqQ3 = get_sequence(chainQ3)

    align = align_sequences(seqQ1, seqQ3)

    sim = compare_structures(str1, str2)
