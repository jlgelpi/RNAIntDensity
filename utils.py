''' common utilities and constants'''


import re

from Bio.PDB.Superimposer import Superimposer
import biobb_structure_checking.modelling.utils as mu

RNA_GROUP_IDs = ['A', 'C', 'G', 'U']
DNA_GROUP_IDs = ['DA', 'DC', 'DG', 'DT']
NA_GROUP_IDs = RNA_GROUP_IDs + DNA_GROUP_IDs

X3DNA_ATOM_CONVERSION = {
    'OP1': 'O1P',
    'OP2': 'O2P'
}

OMEGA_DH_ATOMS = {
    'pyr': ['O4\'', 'C1\'', 'N1', 'C2'],
    'pur': ['O4\'', 'C1\'', 'N9', 'C4']
}

OMEGA_GROUPS = {

}

RIBOSE_RING_ATOMS = ['C4\'', 'O4\'', 'C1\'', 'C2\'', 'C3\'']
BASE_ORIENT_ATOMS = ['N1', 'N9', 'C1', 'C6']

SUP_ATOMS = RIBOSE_RING_ATOMS + BASE_ORIENT_ATOMS

HB_QUALITY_ALL = ['standard', 'acceptable', 'questionable']

NULL_ROTRAN = [
    [[1.0, 0.0, 0.0],
     [0.0, 1.0, 0.0],
     [0.0, 0.0, 1.0]],
    [0.0, 0.0, 0.0]
]

STD_RESNAMES = [
    'A', 'C', 'G', 'U', 'DT', 'DA', 'DC', 'DG',
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLY', 'GLN', 'GLU', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
]

STD_NUCLEOTIDES = [
    'A', 'C', 'G', 'U',
    'DA', 'DC', 'DG', 'DT'
]


def get_atom_from_index(st, atom_index):
    ''' Get atom from structure by index '''
    for at in st.get_atoms():
        if at.serial_number == atom_index:
            return at
    return None


def get_atom_from_x3dna_id(st, atom_id):
    ''' Get atom from structure by ID '''
    print(atom_id)
    at_id, res_id = atom_id.split('@')
    chain_id, resid = res_id.split('.')
    if not '/' in resid:
        match = re.match(r"(.*[A-Za-z]+)(\d+)", resid)
    else:
        match = re.match(r"(.+)\/(\d+)", resid)
    if match:
        resname = match.group(1)  # 'GNG'
        resnum = match.group(2)  # '101'
        hetlb = ' '
        if resname not in STD_RESNAMES:
            hetlb = f"H_{resname}"
        print(st[0][chain_id].child_dict.keys())
        try:
            res = st[0][chain_id][(hetlb, int(resnum), ' ')]
        except KeyError:
            print(f"Error: {chain_id} {resname} {resnum} in structure {st.id}, trying with HETATM label")
            hetlb = 'H_' + resname
            res = st[0][chain_id][(hetlb, int(resnum), ' ')]
    else:
        return None
    print(at_id, res_id, chain_id, resid, resname, resnum, hetlb)
    print(res.child_dict)
    alt_id_num = None
    if '.' in at_id:
        at_id, alt_id_num = at_id.split('.')
    if at_id not in res.child_dict:
        print(f"Warning: Atom {at_id} not found in residue {res_id}, trying conversion")
        at_id = X3DNA_ATOM_CONVERSION.get(at_id, at_id)
    if alt_id_num is not None:
        atm = res.child_dict[at_id].child_dict[alt_id_num]
    else:
        atm = res.child_dict[at_id]
    return atm

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

def calc_omega(res):
    """ Calculate omega dihedral angle for a residue in a structure """
    if mu.is_purine(res):
        at1 = res['O4\'']
        at2 = res['C1\'']
        at3 = res['N9']
        at4 = res['C4']
    elif mu.is_pyrimidine(res):
        at1 = res['O4\'']
        at2 = res['C1\'']
        at3 = res['N1']
        at4 = res['C2']
    else:
        raise ValueError(f"Unsupported residue type for omega calculation: {res.get_resname()}")
    return mu.calc_bond_dihedral(at1, at2, at3, at4)

