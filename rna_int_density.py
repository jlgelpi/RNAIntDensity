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

SUP_ATOMS = RIBOSE_RING_ATOMS + ['N1', 'N9', 'C1', 'C6']

HB_QUALITY_ALL = ['standard', 'acceptable', 'questionable']

HB_QUALITY = ['standard', 'acceptable']

NULL_ROTRAN = [
    [[1.0, 0.0, 0.0],
     [0.0,  1.0, 0.0],
     [0.0, 0.0, 1.0]],
    [0.0, 0.0, 0.0]
]

STD_RESNAMES = [
    'A', 'C', 'G', 'U',
    'DT', 'DA', 'DC', 'DG',
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLY', 'GLN', 'GLU', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
]

STD_NUCLEOTIDES = [
    'A', 'C', 'G', 'U',
    'DA', 'DC', 'DG', 'DT'
]

NEW_CHAIN_ID = 'A'  # New chain ID for the output structure
WAT_CHAIN_ID = 'W'  # Chain ID for water molecules in the output structure


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

        if hb['donAcc_type'] not in HB_QUALITY:
            continue
        at1 = get_atom_from_x3dna_id(st, hb['atom1_id'])
        at2 = get_atom_from_x3dna_id(st, hb['atom2_id'])
        res1 = at1.get_parent()
        res2 = at2.get_parent()
        if res1 == res2:
            continue

        if res1 not in hb_data:
            hb_data[res1] = []

        hb_data[res1].append({
            'atm': at1,
            'atm2': at2,
            'coords': at2.coord,
            'details': hb
        })

        if res2 not in hb_data:
            hb_data[res2] = []

        hb_data[res2].append({
            'atm': at2,
            'atm2': at1,
            'coords': at1.coord,
            'details': hb
        })
    return hb_data

class ResidueGroup():
    ''' Class to hold a group of residues '''
    def __init__(self, resname):
        self.id = resname
        self.structure = Structure(resname)
        self.residues = []
        self.nmod = 0
        self.rmsd_sup = None
        self.rmsd_all = None
        self.rotran = None

    def add_residue(self, res, hb_data):
        ''' Add a residue to the group '''

        residue = res.copy()
        new_chain = Chain(NEW_CHAIN_ID)
        new_chain.add(residue)
        new_mod = Model(self.nmod)
        new_mod.add(new_chain)
        new_wat_chain = Chain(WAT_CHAIN_ID)
        new_mod.add(new_wat_chain)
        self.structure.add(new_mod)
        self.residues.append({
            'residue': res,
            'nmod': self.nmod,
            'hbcoords': [hb['coords'] for hb in hb_data]
        })
        self.nmod += 1

    def superimpose(self):
        ''' Superimpose the models in the group '''
        self.rotran, self.rmsd_sup, self.rmsd_all = superimpose_models(self.structure)

    def transform_hbcoords(self):
        ''' Transform the hydrogen bond coordinates using groups's rotation and translation
            and add them as fake water molecules in the structure
        '''
        nmod = 0
        for i, res_data in enumerate(self.residues):
            nwat = 1
            for hb in res_data['hbcoords']:
                new_atom = Atom(
                    'O',
                    apply_transformation(hb, self.rotran[i]),
                    1.0, 1.0, ' ', 'O', 0, 'O'
                )
                new_residue = Residue(('W', nwat, ' '), 'WAT', ' ')
                new_residue.add(new_atom)
                self.structure[nmod][WAT_CHAIN_ID].add(new_residue)
                nwat += 1
            nmod += 1

    def save(self, output_file_name):
        ''' Save the group structure to a PDB file '''
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(output_file_name)

def prepare_groups(st, hb_data):
    ''' Prepare groups of residues to accumulate HB data '''
    groups = {}
    nresidues = 0
    for res in st.get_residues():
        if res not in hb_data or not hb_data[res] or res.get_resname() not in STD_NUCLEOTIDES:
            continue

        # Create a new group if it does not exist
        if res.get_resname() not in groups:
            groups[res.get_resname()] = ResidueGroup(res.get_resname())

        # Add the residue to the group

        groups[res.get_resname()].add_residue(res, hb_data[res])

        nresidues += 1

    print(f"Found {len(hb_data)} hydrogen bonds on {nresidues} residues")
    return groups

def apply_transformation(crds, rotran):
    ''' Apply the rotation and translation to the coordinates '''
    rot, tran = rotran
    transformed_crds = np.dot(crds, rot) + tran
    return transformed_crds

def process_structure(st, hbonds_data, output_folder, verbose=False):
    ''' Process the structure and X3DNA data to create a PDB file with fake WAT molecules'''
    hb_data = process_hbonds_data(st, hbonds_data)
    if not hb_data:
        print("No hydrogen bonds data found.")
        sys.exit()
    if verbose:
        for res, hbs in hb_data.items():
            for hb in hbs:
                print(
                    f"#DEBUG HB: {mu.atom_id(hb['atm']):<13}"
                    f"- {mu.atom_id(hb['atm2']):<13} {hb['details']['distance']:.3f}"
                )
    # Create a dictionary to hold residue groups and process them
    groups = prepare_groups(st, hb_data)

    for gr in groups.values():
        print(f"Processing group: {gr.id}, number of models: {len(gr.structure)}")
        gr.superimpose()
        print(
            f"Group: {gr.id}, RMSd (all): {max(gr.rmsd_all):.2f}, RMSd (sup): {max(gr.rmsd_sup):.2f}"
        )
        # transform hbcoords
        gr.transform_hbcoords()
        # Save the group structure to a PDB file
        gr.save(f"{output_folder}/{st.id}_{gr.id}.pdb")


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
    parser.add_argument('--no_filter', action='store_true', help='Filter for quality of hydrogen bonds')
    parser.add_argument('--pdb_id', type=str, help='Input structure id', required=True)
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Extra information during execution'
    )
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

    process_structure(structure, x3dna_data['hbonds'], args.output_folder, args.verbose)

if __name__ == '__main__':
    main()
