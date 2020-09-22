'''
set up enviroment:
conda create -n flask python=3.6
pip install flask
pip install flask-wtf
pip install flask-sqlalchemy
pip install flask-migrate
pip install flask-login
pip install flask-bootstrap

pip install email-validator   # in order to have Email validator 



for flask-migrate:   # need to run migrate when you change db related attribute in app.model class 
    flask db init    # will generate a migration directory under top-level folder
    flask db migrate -m 'users table'  # generate migration, not database, -m is comment info, optional
    flask db upgrade  # generate app.db database, in SQLite, it's fine, in MySQL or postgreSQL, need to generate database before apply upgrade
    flask db downgrade # undo last upgrade operation
'''


'''
So this microblog.py is the entry point when running flask run,
then when import app, __init__.py will be invoked

'''

from app import app,db   # from app directory/package import app instance
from app.models import User,Post


@app.shell_context_processor    # when flask shell, this function will be called
def make_shell_context():
    return {'db':db, 'User':User,'Post':Post}  # type db or User, Post to reference cognate object

'''
always remember to register an environmental variable
    export FLASK_APP=microblog.py
after that you can run by:
    flask run   # have to make sure microblog.py is in current folder
server will run on localhost, which is own computer and listen for connection at port 5000,
this means in order to view, just type localhost:5000 and when you operate on this port like click or submit form, server will listen and respond accordiningly


flask shell



'''