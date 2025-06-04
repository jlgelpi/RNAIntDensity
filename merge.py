''' Merge individual residue groups '''

import argparse
import os
import sys

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

import utils as ut


def load_pdb(pdb_id, group_id, input_folder, parser):
    ''' Load a PDB file and return its structure '''
    pdb_file_path = os.path.join(input_folder, f"{pdb_id}_{group_id}.pdb")
    if not os.path.isfile(pdb_file_path):
        print(f"Unexistent PDB file: {pdb_file_path}, skipping")
        return None
    structure = parser.get_structure(group_id, pdb_file_path)
    return structure


def save_structure(structure, file_path):
    ''' Save a PDB structure to a file '''
    io = PDBIO()
    io.set_structure(structure)
    io.save(file_path)
    print(f"Saved structure to {file_path}")


def main():
    ''' Main function to merge individual residue groups '''
    parser = argparse.ArgumentParser(description='Merge individual residue groups from PDB files')

    parser.add_argument(
        '-i', '--input_folder',
        type=str,
        help='Input Folder containing PDB files (Default current folder)',
        default='.'
    )
    parser.add_argument(
        '-o', '--output_folder',
        type=str,
        help='Output Folder for merged PDB file (Default input folder)',
        default=None
    )
    parser.add_argument(
        '--pdb_id_list', type=str, help='List of IDs for the merged structure', required=True
    )
    args = parser.parse_args()

    if args.output_folder is None:
        args.output_folder = args.input_folder

    # Ensure the output folder exists
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Prepare Master strucure from
    pdb_id_list = []
    with open(args.pdb_id_list, 'r') as pdb_id_file:
        for line in pdb_id_file:
            pdb_id = line.strip()
            pdb_id_list.append(pdb_id)

    print(f"Loaded {len(pdb_id_list)} PDB ids from {args.pdb_id_list}")
    # Initialize PDB parser
    parser = PDBParser(QUIET=True)
    master_structures = {}
    master_structure_id = {}
    for group_id in ut.RNA_GROUP_IDs:
        i = 0
        while group_id not in master_structures:
            structure = load_pdb(pdb_id_list[i], group_id, args.input_folder, parser)
            if structure is None:
                if i >= len(pdb_id_list):
                    print(f"No structures found for group {group_id}, skipping")
                    break
                i += 1
            else:
                master_structures[group_id] = structure
                master_structure_id[group_id]= pdb_id_list[i]

    for pdb_id in pdb_id_list:
        for group_id in ut.RNA_GROUP_IDs:
            if group_id not in master_structures:
                print(f"Group {group_id} not found in master structures, skipping")
                continue
            if pdb_id == master_structure_id[group_id]:
                print(f"Skipping {pdb_id}_{group_id} as it is already in master structures")
                continue
            structure = load_pdb(pdb_id, group_id, args.input_folder, parser)
            if structure is not None:
                for model in structure:
                    mod = model.copy()
                    mod.id = len(master_structures[group_id])
                    mod.serial_num = mod.id
                    master_structures[group_id].add(mod)
                    print(f"Added model from {pdb_id}_{group_id} for group {group_id}")
            else:
                print(f"Failed to load structure for {pdb_id}_{group_id}, skipping")

    for group_id, structure in master_structures.items():
        rotran, rmsd_sup, rmsd_all = ut.superimpose_models(structure)
        print(f"Superimposed group {group_id}: RMSD (sup): {max(rmsd_sup):.2f}, RMSD (all): {max(rmsd_all):.2f}")
        output_file = os.path.join(args.output_folder, f"Merged_{group_id}.pdb")
        save_structure(structure, output_file)
        print(f"Saved merged structure for {group_id} to {output_file}")


if __name__ == "__main__":
    main()
