''' Script to combine X3DNA analysis on a 3D grid '''


import json
import argparse
import re
import os
import sys
import numpy as np

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBIO import PDBIO

import biobb_structure_checking.modelling.utils as mu

RIBOSE_RING_ATOMS = ['C4\'', 'O4\'', 'C1\'', 'C2\'', 'C3\'']

SUP_ATOMS = RIBOSE_RING_ATOMS

HB_QUALITY_ALL = ['standard', 'acceptable', 'questionable']

HB_QUALITY = ['standard', 'acceptable']

NULL_ROTRAN = [
    [[1.0, 0.0, 0.0],
     [0.0,  1.0, 0.0],
     [0.0, 0.0, 1.0]],
    [0.0, 0.0, 0.0]
]

STD_RESNAMES = ['A', 'C', 'G', 'U', 'DT', 'DA', 'DC', 'DG']

NEW_CHAIN_ID = 'A'  # New chain ID for the output structure

def get_atom_from_index(st, atom_index):
    ''' Get atom from structure by index '''
    for at in st.get_atoms():
        if at.serial_number == atom_index:
            return at
    return None


def get_atom_from_x3dna_id(st, atom_id):
    ''' Get atom from structure by ID '''
    at_id, res_id = atom_id.split('@')
    chain_id, resid = res_id.split('.')
    match = re.match(r"([A-Za-z]+)(\d+)", resid)
    if match:
        resname = match.group(1)  # 'GNG'
        resnum = match.group(2)  # '101'
        hetlb = ' '
        if resname not in STD_RESNAMES:
            hetlb = f"H_{resname}"
        res = st[0][chain_id][(hetlb, int(resnum), ' ')]
    else:
        return None
    return res.child_dict[at_id]


def atom_renumbering(st):
    """ Sets  Bio.PDB.Atom.serial_number for all atoms in the structure,
        overriding original if any.
    """
    i = 1
    for atm in st.get_atoms():
        atm.serial_number = i
        if hasattr(atm, 'selected_child'):
            atm.selected_child.serial_number = i
        i += 1


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

def process_hbonds_data(st, hbonds_data):
    ''' Process the X3DNA hydrogen bonds data '''
    hb_data = {}
    for hb in hbonds_data:
        at1 = get_atom_from_x3dna_id(st, hb['atom1_id'])
        at2 = get_atom_from_x3dna_id(st, hb['atom2_id'])

        res1 = at1.get_parent()
        if res1 not in hb_data:
            hb_data[res1] = []

        hb_data[res1].append({
            'atm': at1,
            'coords': at2.coord,
            'details': hb
        })

        res2 = at2.get_parent()
        if res2 not in hb_data:
            hb_data[res2] = []

        hb_data[res2].append({
            'atm': at2,
            'coords': at1.coord,
            'details': hb
        })
    return hb_data


def prepare_groups(st, hb_data):
    ''' Prepare groups of residues to accumulate HB data '''
    groups = {}
    nresidues = 0
    for res in st.get_residues():
        if res not in hb_data:
            continue
        if res.id[0] != ' ':
            continue
        if res.get_resname() not in groups:
            groups[res.get_resname()] = {
                'structure': Structure(res.get_resname()),
                'residues': [],
                'nmod': 0,
                'rmsd_sup': None,
                'rmsd_all': None,
                'rotran': None,
                'transformed_hbcoords': []
            }

        residue = res.copy()
        new_chain = Chain(NEW_CHAIN_ID)
        new_chain.add(residue)
        new_mod = Model(groups[res.get_resname()]['nmod'])
        new_mod.add(new_chain)
        groups[res.get_resname()]['structure'].add(new_mod)

        groups[res.get_resname()]['residues'].append({
            'residue' :res,
            'nmod': groups[res.get_resname()]['nmod'],
            'hbcoords': [hb['coords'] for hb in hb_data[res]]
        })

        groups[res.get_resname()]['nmod'] += 1
        nresidues += 1

    print(f"Found {len(hb_data)} hydrogen bonds on {nresidues} residues")
    return groups

def save_group(gr, gr_id, st_id, output_folder):
    ''' Save the group structure to a PDB file'''
    io = PDBIO()
    io.set_structure(gr['structure'])
    io.save(f"{output_folder}/{st_id}_{gr_id}.pdb")


def process_structure(st, hbonds_data, output_folder):
    ''' Process the structure and X3DNA data to create a PDB file with fake WAT molecules'''
    hb_data = process_hbonds_data(st, hbonds_data)
    if not hb_data:
        print("No hydrogen bonds data found.")
        sys.exit()

    # Create a dictionary to hold residue groups
    groups = prepare_groups(st, hb_data)

    for gr_id, gr in groups.items():
        print(f"Processing group: {gr_id}, number of models: {len(gr['structure'])}")
        gr['rotran'], gr['rmsd_sup'], gr['rmsd_all'] = superimpose_models(gr['structure'])
        max_rmsd_all = max(gr['rmsd_all'])
        max_rmsd_sup = max(gr['rmsd_sup'])
        print(
            f"Group: {gr_id}, RMSd (all): {max_rmsd_all:.2f}, RMSd (sup): {max_rmsd_sup:.2f}"
        )
        # transform hbcoords
        nwat = 1
        nmod = 0
        for i, res_data in enumerate(gr['residues']):
            rot, tran = gr['rotran'][i]
            for hb in res_data['hbcoords']:
                transformed_hbcoords = np.dot(hb, rot) + tran
                gr['transformed_hbcoords'].append(transformed_hbcoords)
                new_atom = Atom('O', transformed_hbcoords, 1.0, 1.0, ' ', 'O', 0, 'O')
                new_residue = Residue(('W', nwat, ' '), 'WAT', ' ')
                new_residue.add(new_atom)
                gr['structure'][nmod][NEW_CHAIN_ID].add(new_residue)
                nwat += 1
            nmod += 1
        save_group(gr, gr_id, st.id, output_folder)


def main():
    ''' Main function to process the tructure and X3DNA data '''
    parser = argparse.ArgumentParser(description='Combine X3DNA analysis on a PDB file')

    parser.add_argument(
        '-i','--input_folder', 
        type=str,
        help='Input Folder (Default current folder)',
        default='.'
    )
    parser.add_argument(
        '-o','--output_folder', 
        type=str,
        help='Output Folder (Default input folder)'
    )
    parser.add_argument(
        '--grid_output', help='Output grid file (Default: None)',
        type=str, default=None
    )
    parser.add_argument('--pdb_id', type=str, help='Input structure id', required=True)
    args = parser.parse_args()

    cif_file = f"{args.input_folder}/{args.pdb_id}.cif"
    json_file = f"{args.input_folder}/{args.pdb_id}.json"

    # Check if the input files exist
    if not os.path.isfile(cif_file):
        sys.exit(f"Input mmCIF file not found: {cif_file}")
    if not os.path.isfile(json_file):
        sys.exit(f"Input JSON file not found: {json_file}")

    if args.output_folder is None:
        args.output_folder = args.input_folder
    # Ensure the output folder exists
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Parse the PDB file
    parser = MMCIFParser()
    structure = parser.get_structure(args.pdb_id, cif_file)
    atom_renumbering(structure)
    print(f"Parsed {args.pdb_id} structure from {cif_file}")

    # Load the JSON data
    try:
        with open(json_file, 'r') as f:
            x3dna_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON file: {e}")
    # Validate the JSON data
    if 'hbonds' not in x3dna_data:
        print("JSON file does not contain 'hbonds'")
    print(f"Loaded X3DNA HB data from {json_file}")

    process_structure(structure, x3dna_data['hbonds'], args.output_folder)

if __name__ == '__main__':
    main()
