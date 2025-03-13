# worker.py
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import os
from redis import Redis
from rq import Worker, Queue

listen = ['default']
redis_conn = Redis(host='localhost', port=6379, db=0)

if __name__ == '__main__':
    queues = [Queue(name, connection=redis_conn) for name in listen]
    worker = Worker(queues, connection=redis_conn)
    worker.work()
