# app/__init__.py  – minimal update, keeps ALL of your original logic

import os
import redis
from flask import Flask
from redis import Redis
from rq    import Queue

from app.routes import main as main_blueprint
from app.api    import api           # ← NEW  (JSON job‑queue API)

# ──────────────────────────────────────────────────────────────────────────────
def create_app():
    # Template / static folders (unchanged)
    base_dir      = os.path.dirname(os.path.abspath(__file__))
    template_path = os.path.join(base_dir, "..", "templates")
    static_path   = os.path.join(base_dir, "..", "static")

    app = Flask(
        __name__,
        template_folder=template_path,
        static_folder=static_path
    )
    app.secret_key = "YOUR_SECRET_KEY_CHANGE_THIS"

    # ---- Blueprints ----
    app.register_blueprint(main_blueprint)   # existing UI
    app.register_blueprint(api)              # ← NEW  /api/* endpoints

    # ---- Inject AWS flag into every Jinja template ----
    USE_AWS = os.getenv("USE_AWS") == "1"
    @app.context_processor
    def inject_aws_flag():
        return {"USE_AWS": USE_AWS}

    # ---- Redis / RQ setup ----
    rq_url    = os.environ.get("RQ_REDIS_URL")
    redis_url = rq_url or os.environ.get("REDIS_URL") or "redis://localhost:6379/0"

    app.logger.info("RQ_REDIS_URL = %s", rq_url)
    app.logger.info("REDIS_URL    = %s", os.environ.get("REDIS_URL"))

    # Ping Redis so startup fails loudly if it’s unreachable
    try:
        Redis.from_url(redis_url).ping()
        app.logger.info("Redis ping OK: %s", redis_url)
    except Exception as e:
        app.logger.error("Redis ping FAILED (%s): %s", redis_url, str(e))

    # Expose a shared RQ queue on the Flask app for easy access (NEW)
    app.q = Queue("default", connection=Redis.from_url(redis_url))

    return app
