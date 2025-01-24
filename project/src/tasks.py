from celery_worker import celery
from rdkit import Chem


@celery.task
def substructure_search_task(molecule_list: list, substructure_smiles: str) -> list:
    """
    Задача Celery для поиска по подструктуре.
    """
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if not substructure:
        raise ValueError("Invalid substructure SMILES string.")

    matching_molecules = []
    for smiles in molecule_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule and molecule.HasSubstructMatch(substructure):
            matching_molecules.append(smiles)

    return matching_molecules
