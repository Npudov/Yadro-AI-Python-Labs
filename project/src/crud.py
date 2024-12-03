from sqlalchemy.orm import Session

from task6 import Molecule


# Функция для добавления молекулы
def add_molecule(db: Session, name: str, smiles: str):
    db_molecule = Molecule(name=name, smiles=smiles)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    return db_molecule


# Функция для поиска молекулы по SMILES
def get_molecule_by_smiles(db: Session, smiles: str):
    return db.query(Molecule).filter(Molecule.smiles == smiles).first()


# Функция для получения всех молекул
def get_molecules(db: Session, skip: int = 0, limit: int = 10):
    return db.query(Molecule).offset(skip).limit(limit).all()
