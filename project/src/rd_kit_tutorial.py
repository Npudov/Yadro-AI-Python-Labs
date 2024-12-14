from rdkit import Chem
from rdkit.Chem import Draw

# Basic Molecular Representations
# Create a molecule from a SMILES string
smiles = "CCO"  # Ethanol
molecule = Chem.MolFromSmiles(smiles)
# Draw the molecule
img = Draw.MolToImage(molecule)
img.show()

# Substructure Search
# Define the molecule and the substructure to search for
benzene = Chem.MolFromSmiles("c1ccccc1")
ethanol = Chem.MolFromSmiles("CCO")
# Perform the substructure search
match = ethanol.HasSubstructMatch(benzene)
print("Benzene ring found in ethanol:", match)

# Molecular Visualization
# Create a list of molecules
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
# Draw the molecules in a grid
img = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(200, 200),
                           returnPNG=False)
img.show()
