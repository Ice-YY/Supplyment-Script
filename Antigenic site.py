from Bio import PDB
import numpy as np

def parse_pdb(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    print(f"Structure loaded: {structure.get_id()}")
    return structure

def get_atoms_by_chain(structure, chain_ids):
    atoms = []
    for model in structure:
        for chain in model:
            if chain.get_id() in chain_ids:
                for residue in chain:
                    for atom in residue:
                        atoms.append((chain.get_id(), atom))  # Add chain ID
    return atoms

def calculate_distance(atom1, atom2):
    coord1 = atom1.get_coord()
    coord2 = atom2.get_coord()
    return np.linalg.norm(coord1 - coord2)

def find_contacts(structure, antibody_chains, antigen_chains, distance_threshold=5.0):
    # Find contact points between antibody and antigen chains within a distance threshold.
    antibody_atoms = get_atoms_by_chain(structure, antibody_chains)
    antigen_atoms = get_atoms_by_chain(structure, antigen_chains)
    
    if not antibody_atoms or not antigen_atoms:
        print("No atoms found for the specified chains.")
    
    contacts = set()  # Use a set to deduplicate
    for chain1, atom1 in antibody_atoms:
        for chain2, atom2 in antigen_atoms:
            distance = calculate_distance(atom1, atom2)
            if distance < distance_threshold:
                contacts.add((
                    chain1,  # Add chain ID
                    atom1.get_parent().get_resname(),
                    atom1.get_parent().get_id(),
                    atom1.get_name(),
                    chain2,  # Add chain ID
                    atom2.get_parent().get_resname(),
                    atom2.get_parent().get_id(),
                    atom2.get_name(),
                    distance
                ))
    return contacts

def save_contacts_to_tsv(contacts, output_file):
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("Antibody Chain,Antibody Residue,Antibody Residue ID,Antibody Atom,Antigen Chain,Antigen Residue,Antigen Residue ID,Antigen Atom,Distance (Å)\n")
        for contact in contacts:
            f.write("\t".join(map(str, contact)) + "\n")

def main(pdb_file, antibody_chains, antigen_chains, output_file):
    structure = parse_pdb(pdb_file)
    contacts = find_contacts(structure, antibody_chains, antigen_chains)
    if not contacts:
        print("No contacts found under the threshold.")
    else:
        save_contacts_to_tsv(contacts, output_file)
        print(f"Contacts saved to {output_file}")

# Example: Given antibody and antigen chain identifiers
pdb_file = '4zyp.pdb'  # PDB file path
antibody_chains = ['L', 'M', 'O', 'J', 'K', 'N']  # Antibody chain identifiers
antigen_chains = ['A', 'B', 'C']   # Antigen chain identifiers
output_file = 'contacts.tsv'  # Output TSV file path

main(pdb_file, antibody_chains, antigen_chains, output_file)