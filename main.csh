
# it will ask for molecule input
python input.py
# prepare input for QSharp
qdk-chem convert --from fcidump --to broombridge fcidumpfilename --out data.yml
#Change to GetGateCount directory
cd GetGateCount
# run estimated in QSharp
dotnet run -- --path=../data.yml --format=Broombridge --skip-qubitization --skip-opt-qubitization > out
cp out ../output/msqdk.out
# prepare input for QREChem
python parser.py
# Computes estimates and writes to output.csv
python QREChem.py

