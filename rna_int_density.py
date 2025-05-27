''' Script to combine X3DNA analysis on a 3D grid '''


import json
import argparse
import re

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBIO import PDBIO

import biobb_structure_checking.modelling.utils as mu

SUP_ATOMS = ['C4\'', 'O4\'', 'C1\'', 'C2\'', 'C3\'']

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
        if resname not in ['A', 'C', 'G', 'U', 'T']:
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
    rmsd = []
    rmsd_all = []
    rotran = []

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

def main(st, x3dna_data, args):
    ''' Main function to process the tructure and X3DNA data '''
    groups = {}
    # Prep residue groups
    for res in st.get_residues():
        if res.id[0] != ' ':
            continue
        if res.get_resname() not in groups:
            groups[res.get_resname()] = {
                'residues': [],
                'structure': Structure(res.get_resname()),
                'nmod': 0,
                'rmsd_sup': None,
                'rmsd_all': None,
                'rotran': None
            }
        groups[res.get_resname()]['residues'].append(res)
        residue = res.copy()
        new_chain = Chain('A')
        new_chain.add(residue)
        new_mod = Model(groups[res.get_resname()]['nmod'])
        new_mod.add(new_chain)
        groups[res.get_resname()]['structure'].add(new_mod)
        groups[res.get_resname()]['nmod'] += 1

    for gr in groups:
        groups[gr]['rotran'], groups[gr]['rmsd_sup'], groups[gr]['rmsd_all'] = superimpose_models(groups[gr]['structure'])
        max_rmsd_all = max(groups[gr]['rmsd_all'])
        max_rmsd_sup = max(groups[gr]['rmsd_sup'])
        print(
            f"Group: {gr}, RMSd (all): {max_rmsd_all:.2f}, RMSd (sup): {max_rmsd_sup:.2f}"
    )
        io = PDBIO()
        io.set_structure(groups[gr]['structure'])
        io.save(f"{gr}.pdb")







    # for hb in x3dna_data['hbonds']:
    #     print(
    #         hb['atom1_id'], mu.atom_id(get_atom_from_x3dna_id(st, hb['atom1_id'])),
    #         hb['atom2_id'], mu.atom_id(get_atom_from_x3dna_id(st, hb['atom2_id']))
    #     )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine X3DNA analysis on a 3D grid')
    parser.add_argument('cif_file', type=str, help='Input mmCIF file')
    parser.add_argument('json_file', type=str, help='Input JSON file with X3DNA analysis')
    args = parser.parse_args()

    # Parse the PDB file
    parser = MMCIFParser()
    structure = parser.get_structure('X3DNA', args.cif_file)
    atom_renumbering(structure)
    # Load the JSON data
    with open(args.json_file, 'r') as f:
        x3dna_data = json.load(f)

    main(structure, x3dna_data, args)

