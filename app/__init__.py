from flask import Flask
from app.routes import main as main_blueprint
import os
import redis

def create_app():
    # Set the template folder to the 'templates' directory at the project root
    template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'templates')
    static_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'static')
    app = Flask(__name__, template_folder=template_path, static_folder=static_path)
    app.secret_key = "YOUR_SECRET_KEY_CHANGE_THIS"
    app.register_blueprint(main_blueprint)

    # ← HERE: inject the flag globally
    USE_AWS = os.getenv("USE_AWS") == "1"
    @app.context_processor
    def inject_aws_flag():
        return {"USE_AWS": USE_AWS}

    # Debug Redis configuration
    rq_url = os.environ.get('RQ_REDIS_URL')
    redis_url = rq_url or os.environ.get('REDIS_URL')

    app.logger.info('RQ_REDIS_URL = %s', rq_url)
    app.logger.info('REDIS_URL    = %s', os.environ.get('REDIS_URL'))

    # Ping Redis to verify connectivity
    try:
        redis.Redis.from_url(redis_url).ping()
        app.logger.info('Redis ping OK: %s', redis_url)
    except Exception as e:
        app.logger.error('Redis ping FAILED (%s): %s', redis_url, str(e))

    return app
