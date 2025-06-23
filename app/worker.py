# worker.py  – listen to the “default” queue and use the same Redis URL
import os, sys
from redis import Redis
from rq    import Worker, Queue

# Make the project root importable when this launches in ECS/Fargate
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

redis_url = (
    os.getenv("RQ_REDIS_URL")
    or os.getenv("REDIS_URL")
    or "redis://localhost:6379/0"
)
redis_conn = Redis.from_url(redis_url)

listen = ["default"]

if __name__ == "__main__":
    queues = [Queue(name, connection=redis_conn) for name in listen]
    Worker(queues, connection=redis_conn).work()
