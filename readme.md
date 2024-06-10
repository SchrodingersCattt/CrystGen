Example:
# Perturb and scale structures
python main.py disturb -i example/Antipyrin_912076.cif -p 0.1 -s 0.9 1.0 1.1 -n 3

# Rotate molecules
python main.py rotate -i example/Antipyrin_912076.cif -a 10 -n 2 -r 0 0 1

# Discorder unraveller
python main.py disorder -i example/PAP-H5.cif -o pap-h5
