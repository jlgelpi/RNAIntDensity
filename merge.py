''' Merge individual residue groups '''

import argparse
import os
import sys

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Superimposer import Superimposer

import biobb_structure_checking.modelling.utils as mu

RNA_GROUP_IDs = ['A', 'C', 'G', 'U']
DNA_GROUP_IDs = ['DA', 'DC', 'DG', 'DT']
NA_GROUP_IDs = RNA_GROUP_IDs + DNA_GROUP_IDs
NULL_ROTRAN = [
    [[1.0, 0.0, 0.0],
     [0.0,  1.0, 0.0],
     [0.0, 0.0, 1.0]],
    [0.0, 0.0, 0.0]
]
SUP_ATOMS = ['C4\'', 'O4\'', 'C1\'', 'C2\'', 'C3\'']  # Ribose ring atoms
             
def load_pdb(pdb_id, group_id, input_folder, parser):
    ''' Load a PDB file and return its structure '''
    pdb_file_path = os.path.join(input_folder, f"{pdb_id}_{group_id}.pdb")
    if not os.path.isfile(pdb_file_path):
        print(f"Unexistent PDB file: {pdb_file_path}, skipping")
    structure = parser.get_structure(group_id, pdb_file_path)
    return structure

def save_structure(structure, file_path):
    ''' Save a PDB structure to a file '''
    io = PDBIO()
    io.set_structure(structure)
    io.save(file_path)
    print(f"Saved structure to {file_path}")

def superimpose_models(st):
    ''' superimpose groups of models in a structure '''
    rmsd = [0.]
    rmsd_all = [0.]
    rotran = [NULL_ROTRAN]

    fix_atoms = [
        at
        for at in st[0].get_atoms()
        if at.id in SUP_ATOMS
    ]
    all_atoms = list(st[0].get_atoms())

    spimp = Superimposer()

    for mod in st:
        if mod.id == 0:
            continue
        mov_atoms = [
            at
            for at in st[mod.id].get_atoms()
            if at.id in SUP_ATOMS
        ]
        all_mov_atoms = list(st[mod.id].get_atoms())
        spimp.set_atoms(fix_atoms, mov_atoms)
        spimp.apply(st[mod.id].get_atoms())
        rmsd.append(mu.calc_RMSd_ats(fix_atoms, mov_atoms))
        rmsd_all.append(mu.calc_RMSd_ats(all_atoms, all_mov_atoms))
        rotran.append(spimp.rotran)
    return rotran, rmsd, rmsd_all


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
    for group_id in RNA_GROUP_IDs:
        master_structures[group_id] = load_pdb(pdb_id_list[0], group_id, args.input_folder, parser)
        print(f"Loaded master structure for {group_id} from {pdb_id_list[0]}")
    
    for pdb_id in pdb_id_list[1:]:
        for group_id in RNA_GROUP_IDs:
            structure = load_pdb(pdb_id, group_id, args.input_folder, parser)
            if structure is not None:
                for model in structure:
                    mod = model.copy()
                    mod.id = len(master_structures[group_id])
                    mod.serial_num = mod.id
                    master_structures[group_id].add(mod)
                    print(f"Added model from {pdb_id}_{group_id} for group {group_id}")
    for group_id, structure in master_structures.items():
        rotran, rmsd_sup, rmsd_all = superimpose_models(structure)
        print(f"Superimposed group {group_id}: RMSD (sup): {max(rmsd_sup):.2f}, RMSD (all): {max(rmsd_all):.2f}")
        output_file = os.path.join(args.output_folder, f"Merged_{group_id}.pdb")
        save_structure(structure, output_file)
        print(f"Saved merged structure for {group_id} to {output_file}")

if __name__ == "__main__":
    main()  