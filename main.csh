
# it will ask for molecule input
input.py
# prepare input for QSharp
qdk-chem convert --from fcidump --to broombridge fcidumpfilename --out data.yml
# run estimated in QSharp
dotnet run -- --path=./data.yml --format=Broombridge --skip-qubitization --skip-opt-qubitization
# prepare input for QREChem
parser.py
# Computes estimates and writes to output.csv
QREChem.py
