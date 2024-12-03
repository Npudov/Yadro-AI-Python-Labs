import pytest
from task1 import substructure_search


def test_substructure_found():
    # Проверяем, что молекулы с субструктурой находятся правильно
    molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "c1ccccc1"  # Бензольное кольцо
    result = substructure_search(molecules, substructure)
    assert result == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


def test_no_match():
    # Проверяем, что возвращается пустой список, если совпадений нет
    molecules = ["CCO", "CC(=O)O"]
    substructure = "c1ccccc1"  # Бензольное кольцо
    result = substructure_search(molecules, substructure)
    assert result == []


def test_invalid_smiles():
    # Проверяем, что некорректный SMILES игнорируется
    molecules = ["CCO", "INVALID_SMILES", "c1ccccc1"]
    substructure = "c1ccccc1"  # Бензольное кольцо
    result = substructure_search(molecules, substructure)
    assert result == ["c1ccccc1"]


def test_substructure_empty():
    # Проверяем случай, когда субструктура пуста
    molecules = ["CCO", "c1ccccc1"]
    substructure = ""
    with pytest.raises(ValueError, match="Substructure SMILES cannot be empty"):
        substructure_search(molecules, substructure)


def test_molecule_list_empty():
    # Проверяем случай, когда список молекул пуст
    molecules = []
    substructure = "c1ccccc1"
    with pytest.raises(ValueError, match="Molecule list cannot be empty"):
        substructure_search(molecules, substructure)


def test_full_structure_match():
    # Проверяем, что полное совпадение также считается
    molecules = ["c1ccccc1", "c1ccccc1C"]
    substructure = "c1ccccc1"  # Бензольное кольцо
    result = substructure_search(molecules, substructure)
    assert result == ["c1ccccc1", "c1ccccc1C"]
