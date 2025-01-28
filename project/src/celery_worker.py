from celery import Celery

# Создаем экземпляр Celery
celery = Celery(
    'tasks',  # Имя модуля с задачами
    broker='redis://redis:6379/0',  # Брокер сообщений (Redis)
    backend='redis://redis:6379/0'  # Хранение результатов (Redis)
)

# Настройки Celery
celery.conf.update(
    task_track_started=True,  # Отслеживать начало выполнения задач
    result_expires=3600       # Время жизни результатов (1 час)
)
