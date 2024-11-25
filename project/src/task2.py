from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from os import getenv

import task1

app = FastAPI()

# Хранилище молекул (в оперативной памяти)
molecules = {}


# Модель данных для молекулы
class Molecule(BaseModel):
    id: int
    smiles: str


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


# API: Добавление молекулы
@app.post("/add")
def add_molecule(molecule: Molecule):
    if molecule.id in molecules:
        raise HTTPException(status_code=400, detail="Molecule with this ID already exists.")
    # Проверяем корректность SMILES
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    molecules[molecule.id] = molecule.smiles
    return {"message": "Molecule added successfully."}


# API: Получение молекулы по идентификатору
@app.get("/get/{molecule_id}")
def get_molecule(molecule_id: int):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    return {"id": molecule_id, "smiles": molecules[molecule_id]}


# API: Обновление молекулы по идентификатору
@app.put("/update/{molecule_id}")
def update_molecule(molecule_id: int, molecule: Molecule):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    # Проверяем корректность SMILES
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    molecules[molecule_id] = molecule.smiles
    return {"message": f"Molecule with ID {molecule_id} updated successfully."}


# API: Удаление молекулы по идентификатору
@app.delete("/delete/{molecule_id}")
def delete_molecule(molecule_id: int):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    del molecules[molecule_id]
    return {"message": f"Molecule with ID {molecule_id} deleted successfully."}


# API: Получение всех молекул
@app.get("/list")
def list_molecules():
    return [{"id": mol_id, "smiles": smiles} for mol_id, smiles in molecules.items()]


# API: Поиск по подструктуре
@app.get("/search")
def substructure_search(substructure_smiles: str):
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES string.")

    molecule_list = list(molecules.values())  # Получаем только SMILES из словаря
    matching_molecules = task1.substructure_search(molecule_list, substructure_smiles)  # Возвращает SMILES молекул

    # Генерируем список соответствующих молекул с их идентификаторами
    results = [
        {"id": mol_id, "smiles": smiles}
        for mol_id, smiles in molecules.items()
        if smiles in matching_molecules
    ]
    return results


# [Optional] API: Загрузка файла с молекулами
@app.post("/upload")
def upload_molecules(file: bytes):
    lines = file.decode().splitlines()
    for line in lines:
        try:
            mol_id, smiles = line.split(",")
            if mol_id in molecules:
                continue  # Пропускаем, если молекула уже существует
            if Chem.MolFromSmiles(smiles):
                molecules[mol_id] = smiles
        except ValueError:
            continue  # Пропускаем некорректные строки
    return {"message": "Molecules uploaded successfully."}
