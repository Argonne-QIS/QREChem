from pyscf import scf, gto
from pyscf.tools import fcidump

def set_mol():

    # Get geometry file name from user
    # print ()
    #file_name = input("Enter name of geometry file: ")
    #mol = gto.M(atom=file_name)
    mol = gto.M(atom='input/be2.xyz')

    # Get basis set input from user
    #basis_set = input("Enter the basis set: ")
    #mol.basis = basis_set
    mol.basis = '6-31g*'  #sto-3g, sto-6g, 6-31g*, cc-pvDz, cc-pvTz, cc-pvQz, cc-pv5z

    # Get the molecule charge from user
    #charge = input("Enter molecule charge: ")
    #mol.charge = charge
    mol.charge = 0

    # Build mol object
    mol.build()

    return mol

# Define main function
def main():

    # Read molecule input and return mol object
    mol = set_mol()

    # Read hardware parameters

    # Compute RHF
    rhf = scf.RHF(mol)
    e = rhf.kernel()

    # Generate the FCIDUMP file
    fcidump.from_scf(rhf, 'output/fcidump')

    # Convert the Mole object to a string representation
    mol_string = mol.tostring()

    # Save the string representation to a file
    with open('output/mol.xyz', 'w') as f:
        f.write(mol_string)

# Call main function
if __name__ == "__main__":
    main()
