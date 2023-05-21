
# Set the filename and text to search for
input_filename = "output/msqdk.output"
output_filename = "output/data.csv"
search_rotations = "# Rotations"
search_cnots = "# CNOT"
search_trotter = "Trotter"
basis_set = "sto-3g"
header = "H2,N_orbitals,Rotation,CNOT"

# Open files
input_file = open(input_filename, "r")
output_file = open(output_filename, "w")

# Loop through each line
for line in input_file:
    # Split the line into individual words
    words = line.strip().split(":", 1)
    words_space = line.strip().split()
#    print (words_space)

    # Strip commas from each word
    words = [word.strip(",") for word in words]

    #print (words)

    # Check if the search text is in any of the words
    if search_rotations in words:
        # Print the matching line
        #print(words[1])
        rotations = words[1]
    if search_cnots in words:
        # Print the matching line
        #print(words[1])
        cnots = words[1]

    if search_trotter in words_space[1]:
        # Print the matching line
        #print(words_space[3])
        orbitals = words_space[3]

# Write the matching text to the output file
print(header, file=output_file)
print(basis_set,orbitals,rotations, cnots, sep=',', file=output_file)
output_file.close()

