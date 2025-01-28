import json
from typing import Iterator

from redis.asyncio import Redis
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from os import getenv

import task1
import logging

# Настройка логирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

redis_client = Redis(host='redis', port=6379, db=0)

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
        logger.error(f"Molecule with ID {molecule.id} already exists.")
        raise HTTPException(status_code=400, detail="Molecule with this ID already exists.")
    # Проверяем корректность SMILES
    if not Chem.MolFromSmiles(molecule.smiles):
        logger.error(f"Invalid SMILES string: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    molecules[molecule.id] = molecule.smiles
    logger.info(f"Molecule with ID {molecule.id} added successfully.")
    return {"message": "Molecule added successfully."}


# API: Получение молекулы по идентификатору
@app.get("/get/{molecule_id}")
def get_molecule(molecule_id: int):
    if molecule_id not in molecules:
        logger.error(f"Molecule with ID {molecule_id} not found.")
        raise HTTPException(status_code=404, detail="Molecule not found.")
    logger.info(f"Molecule with ID {molecule_id} retrieved successfully.")
    return {"id": molecule_id, "smiles": molecules[molecule_id]}


# API: Обновление молекулы по идентификатору
@app.put("/update/{molecule_id}")
def update_molecule(molecule_id: int, molecule: Molecule):
    if molecule_id not in molecules:
        logger.error(f"Molecule with ID {molecule_id} not found.")
        raise HTTPException(status_code=404, detail="Molecule not found.")
    # Проверяем корректность SMILES
    if not Chem.MolFromSmiles(molecule.smiles):
        logger.error(f"Invalid SMILES string: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    molecules[molecule_id] = molecule.smiles
    logger.info(f"Molecule with ID {molecule_id} updated successfully.")
    return {"message": f"Molecule with ID {molecule_id} updated successfully."}


# API: Удаление молекулы по идентификатору
@app.delete("/delete/{molecule_id}")
def delete_molecule(molecule_id: int):
    if molecule_id not in molecules:
        logger.error(f"Molecule with ID {molecule_id} not found.")
        raise HTTPException(status_code=404, detail="Molecule not found.")
    del molecules[molecule_id]
    logger.info(f"Molecule with ID {molecule_id} deleted successfully.")
    return {"message": f"Molecule with ID {molecule_id} deleted successfully."}


# API: Получение всех молекул
@app.get("/list")
def list_molecules(skip: int = 0, limit: int = 10) -> Iterator[dict]:
    logger.info(f"Listing molecules with skip={skip} and limit={limit}")

    for i, (mol_id, smiles) in enumerate(molecules.items()):
        if i < skip:
            continue
        if i >= skip + limit:
            break
        yield {"id": mol_id, "smiles": smiles}

    logger.info(f"Returned {min(limit, len(molecules) - skip)} molecules.")

    # return [{"id": mol_id, "smiles": smiles} for mol_id, smiles in molecules.items()]


# API: Поиск по подструктуре
@app.get("/search")
async def substructure_search(substructure_smiles: str):
    logger.info(f"Substructure search requested for: {substructure_smiles}")
    # Ключ для кэша
    cache_key = f"search:{substructure_smiles}"

    # Проверяем кэш
    cached_result = await redis_client.get(cache_key)  # Используем await
    if cached_result:
        logger.info(f"Cache hit for substructure: {substructure_smiles}")
        return json.loads(cached_result)
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if not substructure:
        logger.error(f"Invalid substructure SMILES string: {substructure_smiles}")
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES string.")

    molecule_list = list(molecules.values())  # Получаем только SMILES из словаря
    matching_molecules = task1.substructure_search(molecule_list, substructure_smiles)  # Возвращает SMILES молекул

    # Генерируем список соответствующих молекул с их идентификаторами
    results = [
        {"id": mol_id, "smiles": smiles}
        for mol_id, smiles in molecules.items()
        if smiles in matching_molecules
    ]
    logger.info(f"Found {len(results)} matching molecules.")
    # Сохраняем результат в кэш
    await redis_client.setex(cache_key, 60, json.dumps(results))  # Кэшируем на 60 секунд
    return results


# [Optional] API: Загрузка файла с молекулами
@app.post("/upload")
def upload_molecules(file: bytes):
    logger.info("Uploading molecules from file.")
    lines = file.decode().splitlines()
    for line in lines:
        try:
            mol_id, smiles = line.split(",")
            if mol_id in molecules:
                logger.warning(f"Molecule with ID {mol_id} already exists, skipping.")
                continue  # Пропускаем, если молекула уже существует
            if Chem.MolFromSmiles(smiles):
                molecules[mol_id] = smiles
                logger.info(f"Molecule with ID {mol_id} added from file.")
        except ValueError:
            logger.error(f"Invalid line format: {line}")
            continue  # Пропускаем некорректные строки
    logger.info("Molecules uploaded successfully.")
    return {"message": "Molecules uploaded successfully."}
