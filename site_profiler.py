import sys
# UPDATED IMPORT: Added 'is_aa' directly to the import list
from Bio.PDB import PDBList, PDBParser, NeighborSearch, is_aa

# --- CONFIGURATION ---
TARGET_PDB_ID = "1IEP"  # Gleevec bound to Abl Kinase
LIGAND_NAME = "STI"     # The PDB code for Gleevec (Imatinib)
SEARCH_RADIUS = 5.0     # Distance in Angstroms to define "binding"
# ---------------------

def get_structure(pdb_id):
    """Downloads and parses the PDB file."""
    pdbl = PDBList()
    # Download the file to the current directory
    # 'obsolete=False' ensures we get the standard file
    filename = pdbl.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, filename)
    return structure

def find_binding_pocket(structure, ligand_code, radius):
    """Finds residues interacting with the ligand."""
    model = structure[0] # Get the first model
    
    # 1. Separate Atoms into 'Ligand' and 'Protein'
    ligand_atoms = []
    protein_atoms = []
    
    print(f"Scanning structure for ligand '{ligand_code}'...")
    
    for residue in model.get_residues():
        # Check if this residue is our drug
        if residue.resname == ligand_code:
            for atom in residue:
                ligand_atoms.append(atom)
        # UPDATED CHECK: Use is_aa(residue) directly
        elif is_aa(residue, standard=True):
            for atom in residue:
                protein_atoms.append(atom)
                
    if not ligand_atoms:
        print(f"Error: Ligand '{ligand_code}' not found in structure.")
        sys.exit(1)
        
    print(f"Found {len(ligand_atoms)} atoms in the drug molecule.")
    print(f"Found {len(protein_atoms)} atoms in the protein.")
    
    # 2. Use KD-Tree for Fast Neighbor Search
    print("Building neighbor search tree...")
    ns = NeighborSearch(protein_atoms)
    
    binding_residues = set()
    
    # For every atom in the drug, ask: "Who is near me?"
    for atom in ligand_atoms:
        # search returns a list of nearby atoms
        # level='R' means return the Residue object, not just the Atom
        neighbors = ns.search(atom.coord, radius, level='R') 
        for res in neighbors:
            binding_residues.add(res)
            
    # Sort by residue number (ID)
    return sorted(list(binding_residues), key=lambda x: x.id[1])

if __name__ == "__main__":
    # 1. Get Data
    print(f"Fetching PDB: {TARGET_PDB_ID}...")
    structure = get_structure(TARGET_PDB_ID)
    
    # 2. Analyze
    print(f"\nAnalyzing Binding Site (Radius: {SEARCH_RADIUS} A)...")
    pocket_residues = find_binding_pocket(structure, LIGAND_NAME, SEARCH_RADIUS)
    
    # 3. Output Results
    print("-" * 40)
    print(f"BINDING POCKET FOR {LIGAND_NAME} (Gleevec):")
    print("-" * 40)
    print(f"Total interacting residues: {len(pocket_residues)}")
    print("\nResidue\t ID\t Chain")
    
    with open("binding_site.txt", "w") as f:
        f.write(f"Binding Site for {TARGET_PDB_ID} / {LIGAND_NAME}\n")
        for res in pocket_residues:
            res_name = res.resname
            res_id = res.id[1]
            chain = res.get_parent().id
            
            line = f"{res_name}\t {res_id}\t {chain}"
            print(line)
            f.write(line + "\n")
            
    print("-" * 40)
    print("Success! List saved to 'binding_site.txt'")