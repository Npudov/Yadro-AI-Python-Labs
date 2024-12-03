from rdkit import Chem


def substructure_search(molecule_list, substructure_smiles):
    """
    Находит все молекулы из заданного списка, соответствующие определенной структуре

    Parameters:
        molecule_list (list): Список молекул в формате SMILES.
        substructure_smiles (str): Субструктура в SMILES формате.

    Returns:
        list: Список строк молекул в формате SMILES, содержащие субструктуру.
    """
    if not molecule_list:
        raise ValueError("Molecule list cannot be empty")

    if not substructure_smiles.strip():
        raise ValueError("Substructure SMILES cannot be empty")

    substructure = Chem.MolFromSmiles(substructure_smiles)

    if substructure is None:
        raise ValueError("Invalid substructure SMILES")

    matching_molecules = []

    for smiles in molecule_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule and molecule.HasSubstructMatch(substructure):
            matching_molecules.append(smiles)

    return matching_molecules


# Example usage
molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
substructure = "c1ccccc1"
result = substructure_search(molecules, substructure)
print("Molecules containing the substructure:", result)
