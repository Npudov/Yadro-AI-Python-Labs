from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import os

# Получаем параметры подключения из окружения
user = os.getenv("POSTGRES_USER", "postgres")
password = os.getenv("POSTGRES_PASSWORD", "postgres")
host = os.getenv("POSTGRES_HOST", "localhost")
port = os.getenv("POSTGRES_PORT", "5432")
db_name = os.getenv("POSTGRES_DB", "molecules_db")

# Формируем строку подключения к базе данных
DATABASE_URL = os.getenv(
    "DATABASE_URL",
    f"postgresql+psycopg2://{user}:{password}@{host}:{port}/{db_name}"
)

# Создание соединения с базой данных
engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Настройка базы данных с SQLAlchemy
Base = declarative_base()


class Molecule(Base):
    __tablename__ = 'molecules'

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    smiles = Column(String, unique=True, index=True)


# Создание таблиц в базе данных
def init_db():
    Base.metadata.create_all(bind=engine)


# Зависимость для получения сессии
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
