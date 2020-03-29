from Complex_v2_temporal import *
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description="Defining the input and output directories, stoichiometry, optimization and if chains have the same name across files")

    parser.add_argument('-i', '--input',
                        dest="indir",
                        action="store",
                        default = ".",
                        help="Input directory with all PDB files")

    parser.add_argument('-v', "--verbose",
                        dest="verbose",
                        default = False,
                        help="Print log in stderr",
                        action="store_true")

    parser.add_argument('-opt', "--optimize",
                        dest="optimize",
                        help="Optimize the macro-complex",
                        default = False,
                        action="store_true")

    parser.add_argument('-st', "--stoichiometry",
                        dest="stoichiometry",
                        help="Allows the user to pass the stoichiometry once the chains have been processed.",
                        type = json.loads,
                        default= None,
                        action='store')

    parser.add_argument('-o', '--outdir',
                        dest="outdir",
                        action='store',
                        default = "results",
                        help='Indicate a name to create a directory where the outputs will be written')

    parser.add_argument('-c', '--correct_chains',
                        dest= "cor",
                        action='store_true',
                        default = True,
                        help='Indicate if the name of the chains are correct')

    options = parser.parse_args()

    complex = Complex()

    complex.get_complex(options.indir, options.stoichiometry, options.cor, options.outdir)

if __name__ == "__main__":

    main()
