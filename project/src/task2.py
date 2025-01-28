import logging
from os import getenv
from typing import Iterator

from celery.result import AsyncResult
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from redis.asyncio import Redis

from tasks import substructure_search_task

# Настройка логирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

redis_client = Redis(host='redis', port=6379, db=0)

app = FastAPI()

# Хранилище молекул (в оперативной памяти)
molecules = {}


# Модель для запроса поиска
class SearchRequest(BaseModel):
    substructure_smiles: str


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
@app.get("/list", response_model=None)
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
'''@app.get("/search")
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
    return results'''


# Эндпоинт для запуска задачи поиска
@app.post("/search/start")
async def start_search(request: SearchRequest):
    # Получаем список SMILES из хранилища
    molecule_list = list(molecules.values())

    # Запускаем задачу Celery
    task = substructure_search_task.delay(molecule_list, request.substructure_smiles)
    logger.info(f"Started search task with ID: {task.id}")
    return {"task_id": task.id}


# Эндпоинт для проверки статуса задачи и получения результатов
@app.get("/search/status/{task_id}")
async def get_search_status(task_id: str):
    task_result = AsyncResult(task_id)

    if task_result.state == 'PENDING':
        logger.info(f"Task {task_id} is still processing")
        return {"status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        logger.info(f"Task {task_id} completed successfully")
        # Формируем результат
        matching_molecules = task_result.result
        results = [
            {"id": mol_id, "smiles": smiles}
            for mol_id, smiles in molecules.items()
            if smiles in matching_molecules
        ]
        return {"status": "Task completed", "result": results}
    else:
        logger.error(f"Task {task_id} failed with state: {task_result.state}")
        return {"status": task_result.state, "error": str(task_result.result)}


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
