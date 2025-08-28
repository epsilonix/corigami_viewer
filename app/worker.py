# app/worker.py
import os, sys, platform
from redis import Redis
from rq import Queue
from rq.worker import Worker, SimpleWorker  # import both

# Make the project root importable when this launches in ECS/Fargate
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

redis_url = (
    os.getenv("RQ_REDIS_URL")
    or os.getenv("REDIS_URL")
    or "redis://localhost:6379/0"
)
redis_conn = Redis.from_url(redis_url)
listen = ["default"]

def pick_worker_class():
    # Override via env if you want: RQ_WORKER_CLASS=Worker or SimpleWorker
    forced = os.getenv("RQ_WORKER_CLASS")
    if forced:
        return {"Worker": Worker, "SimpleWorker": SimpleWorker}[forced]
    # Default: avoid forking only on macOS (Darwin)
    if platform.system() == "Darwin" and os.getenv("RQ_ALLOW_FORK", "0") != "1":
        return SimpleWorker
    return Worker

if __name__ == "__main__":
    queues = [Queue(name, connection=redis_conn) for name in listen]
    worker_cls = pick_worker_class()
    worker = worker_cls(queues, connection=redis_conn)
    worker.work()
