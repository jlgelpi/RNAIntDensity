## RNAIntDensity

Scripts to map X3DNA output to a collection of PDBs, with fake water molecules in the interaction sites
usage: rna_int_density.py [-h] [-i INPUT_FOLDER] [-o OUTPUT_FOLDER]
                          [--no_filter] --pdb_id PDB_ID [-v]

Combine X3DNA analysis on a single PDB file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Input Folder (Default current folder)
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output Folder (Default input folder)
  --no_filter           Filter for quality of hydrogen bonds
  --pdb_id PDB_ID       Input structure id
  -v, --verbose         Extra information during execution
usage: merge.py [-h] [-i INPUT_FOLDER] [-o OUTPUT_FOLDER] --pdb_id_list
                PDB_ID_LIST

Merge individual residue groups from PDB files

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Input Folder containing PDB files (Default current
                        folder)
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output Folder for merged PDB file (Default input
                        folder)
  --pdb_id_list PDB_ID_LIST
                        List of IDs for the merged structure
