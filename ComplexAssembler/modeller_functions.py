from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
import os

def get_profile(directory, number):
    log.none()

    env = environ()

    env.libs.topology.read(file='$(LIB)/top_heav.lib') #read topology
    env.libs.parameters.read(file='$(LIB)/par.lib') #read parameters

    filename = directory + "/model_" + str(number) + ".pdb"
    code = directory + "/model_" + str(number)

    #read model files
    mdl = complete_pdb(env, filename)

    zscore = mdl.assess_normalized_dope()

    print("Z-score of model %d is %.4f" %(number, zscore))

    # Asess with DOPE:
    s = selection(mdl) # all atom selection

    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file = code + ".profile", normalize_profile=True, smoothing_window = 30)

def optimize(directory, number):
    log.none()
    env = environ()
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    filename = directory + "/model_" + str(number) + ".pdb"
    #read model files
    code = directory + "/model_" + str(number)
    mdl = complete_pdb(env, filename)


    # Select all atoms:
    atmsel = selection(mdl)

    # Generate the restraints:
    mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)


    # Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='NO_REPORT')
    md = molecular_dynamics(output='NO_REPORT')

    if not os.path.exists(directory + '/optimization_stats'):
        os.makedirs(directory + '/optimization_stats')

    # Open a file to get basic stats on each optimization
    trcfil = open(directory + '/optimization_stats/model_' + str(number) + '.D00000001', 'w')

    # Run CG on the all-atom selection; write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))
    # Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
    # 10 steps during the run, and write stats every 10 steps
    md.optimize(atmsel, temperature=300, max_iterations=50,
                actions=[actions.trace(10, trcfil)])
    # Finish off with some more CG, and write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20, actions=[actions.trace(5, trcfil)])

    #mpdf = atmsel.energy()

    mdl.write(file=code+'.pdb')
