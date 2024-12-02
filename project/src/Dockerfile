FROM python:3.11


# Установка рабочей директории

# Создаём директорию /app
RUN mkdir /app


WORKDIR /app/

# Копирование зависимостей
COPY requirements.txt .
RUN pip install -r requirements.txt

# Копирование проекта
COPY . .

# Открытие порта приложения
EXPOSE 8000

# Команда запуска приложения
CMD ["uvicorn", "task2:app", "--host", "0.0.0.0", "--port", "8000"]