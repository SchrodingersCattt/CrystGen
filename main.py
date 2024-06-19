import argparse
from mk_disturb import StructurePerturber
from mk_rotation import MoleculeRotator
from unravel_disorder import DisorderUnraveller


def parse_args():
    parser = argparse.ArgumentParser(description='Molecular Structure Manipulation')

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Arguments for mk_disturb.py
    disturb_parser = subparsers.add_parser('disturb', help='Perturb and scale structures')
    disturb_parser.add_argument('-i', '--input', type=str, required=True, help='Input CIF file')
    disturb_parser.add_argument('-p', '--perturb', type=float, default=0.02, required=True, help='Perturb distance')
    disturb_parser.add_argument('-s', '--scale', type=float,  default=0.1,
                                help='Scaling factors for the cell')
    disturb_parser.add_argument('-t', '--maxtilt', type=float, default=0.01,
                                help='Tilt factors for the cell')
    disturb_parser.add_argument('-n', '--number', type=int, default=5, help='Number of disturbed replicas')

    # Arguments for mk_rotation.py
    rotate_parser = subparsers.add_parser('rotate', help='Rotate molecules in the structure')
    rotate_parser.add_argument('-i', '--input', type=str, required=True, help='Input CIF file')
    rotate_parser.add_argument('-a', '--angle', type=float, default=5, help='Rotation angle')
    rotate_parser.add_argument('-r', '--rotate', type=float, nargs='+', default=[0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1],
                               help='Rotation axis directions')
    rotate_parser.add_argument('-n', '--number', type=int, default=10, help='Number of rotation replicas')

    # Arguments for disorder_unraveller.py
    disorder_parser = subparsers.add_parser('disorder', help='Unravel disorder in structures')
    disorder_parser.add_argument('-i', '--input', type=str, required=True, help='Input CIF file')
    disorder_parser.add_argument('-o', '--output', type=str, required=True, help='Output filename prefix')

    return parser.parse_args()

def main():
    args = parse_args()
    if args.command == 'disturb':
        perturber = StructurePerturber(args.input, args.perturb, args.scale, args.maxtilt, args.number)
        perturber.perturb_atoms()
        perturber.scale_cells()
        perturber.tilt_cells()
    elif args.command == 'rotate':
        rotate_axes = [list(map(int, args.rotate[i:i+3])) for i in range(0, len(args.rotate), 3)]
        rotator = MoleculeRotator(args.input, args.angle, rotate_axes, args.number)
        rotator.rotate_molecules()
    elif args.command == 'disorder':
        unraveller = DisorderUnraveller(args.input)
        unraveller.write_structure(args.output)

if __name__ == '__main__':
    main()