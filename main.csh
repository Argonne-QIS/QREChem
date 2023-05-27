
# it will ask for molecule input
python input.py
# prepare input for QSharp
qdk-chem convert --from fcidump --to broombridge output/fcidump --out output/fcidump.yml
#Change to GetGateCount directory
cd /home/yuri/Downloads/Quantum/samples/chemistry/GetGateCount
# run estimated in QSharp
dotnet run -- --path=/home/yuri/QREChem-Otten-Tests/output/fcidump.yml --format=Broombridge --skip-qubitization --skip-opt-qubitization > /home/yuri/QREChem-Otten-Tests/output/msqdk.output
cd /home/yuri/QREChem-Otten-Tests
# prepare input for QREChem
python parser.py
# Computes estimates and writes to output.csv
python QREChem.py

