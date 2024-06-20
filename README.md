# CrystGen

## Example:
- Perturb and scale structures
```
python main.py disturb -i example/Antipyrin_912076.cif -p 0.05 -s 0.1 -t 0.05 -n 10
```

- Rotate molecules
```
python main.py rotate -i example/Antipyrin_912076.cif -a 10 -n 2 -r 0 0 1
```

- Supercell
```
python main.py supercell -i example/Antipyrin_912076.cif -l 2 2 2 
python main.py supercell -i example/Antipyrin_912076.cif -s 2
```

- Discorder unraveller
```
python main.py disorder -i example/PAP-H5.cif -o pap-h5
```
