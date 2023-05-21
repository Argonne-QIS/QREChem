# !pip install lmfit
import numpy as np
#import matplotlib.pyplot as plt
#from IPython.display import Image
import pandas as pd
import math
from pyscf import gto, scf, mcscf, lo, cc,ci
import itertools
import json
#from matplotlib import markers
#from lmfit.models import PolynomialModel

# read hardware specs from file
def read_hardware_specs(hardware_specs_filename):

    # device properties
#    cnot_fid = 1 - 0.0018 # Falcon R10 IBM device
#    av_gate_time = 450e-9 # CNOT gate time
#    success_prob = 0.99
#    precision = 1e-4 # target precision in Ha (decimal)
#    trot_error = 1.6e-3 # Chemical accuracy in Hartree

    with open(hardware_specs_filename, 'r') as file:
        json_data = json.load(file)

    # Access the data from the JSON
    cnot_fid = json_data['cnot_fid']
    av_gate_time = json_data['av_gate_time']
    success_prob = json_data['success_prob']
    precision = json_data['precision']
    trot_error = json_data['trot_error']

    failure_prob = 1 - success_prob
    n = round(-(np.log2(precision)))
    print("Converting to binary precision:\nNeed at least {} qubits to reach {} precision".format(n, precision))

    return cnot_fid, av_gate_time, success_prob, precision, trot_error, failure_prob, n

#In the following calculations that involve computation of error rates we assume that gate fidelity and error rate for a gate are equal. However, it is only true when only decoherent (depolarizing noise) is present. Coherent noise can significantly increase errors. More details can be found in AVS Quantum Sci. 4, 013801 (2022). All calculations are per 1 shot

def compute_t_known_eig(n, failure_prob):
    # for QPE when eigenvector is known
    t = (n + np.log2(2 + 1 / (2*failure_prob)))
    t = np.around(t, decimals=0)
    return t

def compute_t_unknown_eig(n, delta, dE, failure_prob):
    # for QPE when eigenvector is not known
    t = (n + np.log2(2 + delta**2 / (2*failure_prob*dE**2)))
    t = np.around(t, decimals=0)
    return t

def compute_t_trot(n, dE, failure_prob, trot_error):
    # for QPE when eigenvector Trotterization is used
    t = n + np.log2(2 + trot_error**2 / (2*failure_prob*dE**2))
    t = np.around(t, decimals=0)
    return t

def compute_depth(n_cnot, n_rot, n_orb):
    # very approximate formula, can be improved
    return round(n_rot/n_orb + n_cnot)

def compute_success_prob_from_fid(fid, depth):
    # simple success probability based on average gate fidelity and total depth
    return fid**depth

def compute_exec_time_from_gate_time(av_gate_time, depth):
    # total execution times based on provided average gate time
    return av_gate_time * depth

def compute_n_TS(n_orb):
    # Function to compute the number of Trotter steps required
    # N^2.5 upper limit from:
    # Reiher, M., Wiebe, N., Svore, K. M., Wecker, D. & Troyer,
    # M. . Proc. Nat. Acad. Sci. 114, 7555–7560 (2017).

    return round(n_orb**2.5)

def compute_req_err_rate(success_prob, depth):
    # based on desired success probability computes required average gate error rate
    return 1 - success_prob**(1 / depth)

def compute_d_surf_code(req_err_rate, phys_err_rate):
    # Computation of distance d for surface code
    # equation (2) from AVS Quantum Sci. 4, 013801 (2022)

    return round(math.log(10 * req_err_rate, 100 * phys_err_rate)*2 - 1)

def compute_dE(bases, atom):
    # Computing the energy difference dE between the trial state (HF) and the first excited state
    # Using CCSD for ground state and EOM-CCSD for the excited state (see Figure above)
    
    dE = []
    delta = []
    for basis in bases:
#         print("\nbasis:", basis)
        mol = gto.Mole(verbose=0, atom=atom, charge=0, spin=0, basis=basis)

        mol.build()
        mf = scf.RHF(mol)
        mf.kernel()
        mycc = cc.CCSD(mf).run()
        e_corr = abs(mycc.e_corr)
        e_ee, c_ee = mycc.eeccsd(nroots=1)

        dE.append(e_ee - e_corr)
        delta.append(abs(e_corr))
#         print("Vertical excitation energy:", e_ee)
        
    return dE, delta

def process_molecule(df, atom, filename, cnot_fid, av_gate_time, success_prob, precision, trot_error, failure_prob, n):
    # Doing analysis for a given molecule and putting it in Pandas dataframe

    print("Molecular geometry:\n", atom)
    bases = df.index.to_list()
    print("Read the data for these basis sets:", bases)
    mol_array = [filename.split('.')[0] for n in range(len(bases))]
    df['mol'] = mol_array

    dE, delta = compute_dE(bases, atom)
    df['dE'] = dE
    df['delta'] = delta
#     print(df)

    df['depth_1TS'] = df.apply(lambda row : compute_depth(row['CNOT'],
                         row['Rotation'], row['N_orbitals']), axis = 1)
    df['success_prob_1TS'] = df.apply(lambda row : compute_success_prob_from_fid(cnot_fid,
                         row['depth_1TS']), axis = 1)
    df['exec_time_1TS'] = df.apply(lambda row : compute_exec_time_from_gate_time(av_gate_time,
                         row['depth_1TS']), axis = 1)
    df["n_TS_upper"] = df.apply(lambda row : compute_n_TS(row['N_orbitals']), axis = 1)

    df['N_ancilla'] = df.apply(lambda row : compute_t_trot(n, row['dE'],
                                                           failure_prob, trot_error), axis = 1)

    df['Rotation_tot'] = round(df['Rotation'] * df['n_TS_upper'] * df['N_ancilla'])
    df['CNOT_tot'] = round(df['CNOT'] * df['n_TS_upper'] * df['N_ancilla'])
    df['depth_tot'] = round(df['depth_1TS'] * df['n_TS_upper'] * df['N_ancilla'])
    df['exec_time_tot'] = df['exec_time_1TS'] * df['n_TS_upper'] * df['N_ancilla']
    df['required_rate'] = df.apply(lambda row : compute_req_err_rate(success_prob,
                         row['depth_tot']), axis = 1)

    df['d_surf_code'] = df.apply(lambda row : compute_d_surf_code(row['required_rate'],
                         1-cnot_fid), axis = 1)
    df['n_phys_qubits'] = df['d_surf_code']**2 * df['N_orbitals']
# the number of T gates from rotation synthesis, 10+(12*log2(1/epsilon)), epsilon=10^(-9), arxiv:1212.6253
    df['T_tot'] = round(df['Rotation_tot'] * 368.7682342)

#     df

# Here we read the csv files with the precomputed gate count using Q#. Potentially Q# has python wrapper so everything can be done here. 

# Define main function
def main():

    # Read hardware specs from the file
    hardware_specs_filename = "input/hardware-specs.json"
    cnot_fid, av_gate_time, success_prob, precision, trot_error, failure_prob, n = read_hardware_specs(hardware_specs_filename)

    filename = 'output/data.csv'
#    r = 0.74
#    atom = 'H .0 .0 .0; H .0 .0 {}'.format(r)

    # Read the molecule string from the file
    with open('output/mol.xyz', 'r') as f:
        mol_string = f.read()
    atom = mol_string

    df = pd.read_csv(filename, index_col=0)
    process_molecule(df, atom, filename, cnot_fid, av_gate_time, success_prob, precision, trot_error, failure_prob, n)

    df.to_csv(r'output/QREChem-output.cvs')
    #print(df)

# Call main function
if __name__ == "__main__":
    main()

