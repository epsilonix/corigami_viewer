from flask import Flask
from app.routes import main as main_blueprint
import os

def create_app():
    # Set the template folder to the 'templates' directory at the project root
    template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'templates')
    static_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'static')
    app = Flask(__name__, template_folder=template_path, static_folder=static_path)
    app.secret_key = "YOUR_SECRET_KEY_CHANGE_THIS"
    app.register_blueprint(main_blueprint)
    return app
