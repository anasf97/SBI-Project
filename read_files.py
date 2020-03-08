from Bio.PDB import *

def read_pdb(file):
    parser = PDBParser()
    filext = file.split(".")
    pdb_id = filext[0]
    structure = parser.get_structure(pdb_id , file)


def compare_structures(structure, other):
    superimposer = Superimposer()

    str_model = structure[0]
    other _model = other[0]

    str_atoms = []
    other_atoms = []

    similarity = {}

    chain_pairs = [ (x, y) for x in str_model.keys() for y in other_model.keys()]

    for pair in chain_pairs:

        for str_res in str_model[pair[0]]:

            str_atoms.append(str_res['CA'])

        for other_res in other_model[pair[1]]:

            other_atoms.append(other_res['CA'])

        superimposer.set_atoms(str_atoms, other_atoms)

        if superimposer.rms <= 2:

            similarity[1] = []

            similarity[1].append(pair)

        else:
            similarity[0] = []

            similarity[0].append(pair)
        #superimposer.apply(other_model.get_atoms())

    return similarity

def superimpose(struture, other):

    for pair in similarity[1]:



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Take input and output files...")

    parser.add_argument("-i", "--input", dest="infile", action="store", default= os.getcwd(), help="Input FASTA formatted file")

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Print log in stderr")

    options = parser.parse_args()

    if os.path.isdir(options.infile):
        filelist = [f for f in os.listdir(options.infile) if f.endswith((".fasta",".fa", ".fa.gz", ".fasta.gz"))]
