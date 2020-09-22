
'''
when you import a directory, the __init__.py file got executed,
deciding which variable/symbols expose to outside world

'''

from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_login import LoginManager
from flask_bootstrap import Bootstrap

app = Flask(__name__)    # an Flask instance called app, pass __name__ variable which I guess is __main__? Flask use this as a starting point
app.config.from_object(Config)
db = SQLAlchemy(app)   # <SQLAlchemy engine=sqlite:////Users/migu7781/Documents/dev/flask/microblog2/app.db>
migrate = Migrate(app,db)
login = LoginManager(app)
login.login_view = 'login'   # which page handle login
bootstrap = Bootstrap(app)

from app import routes,models    # 1. this app means the app directory in microblog folder, which is a package because there's a __init__.py file under this folder
                             # 2. You have to import it after instantiate app because there are some dependency under the hood happening