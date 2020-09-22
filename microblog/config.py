
'''
Seperation of concerns,
a isolated config class to address all possible configuration

'''
import os
basedir = os.path.abspath(os.path.dirname(__file__))  # __file__ is the path of the current program, abspath and realpath will return full path and realpath will convert symbolic link(softlink files if any)

class Config(object):   # after connect to app, you can assess it by app.config['SECRET_KEY]
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'  # to prevent cross-sitereference forgery
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or 'sqlite:///' + os.path.join(basedir, 'app.db')  # app.db file at the same level as config.py
    SQLALCHEMY_TRACK_MODIFICATIONS = False   # turn off a function we don't want for flask-sqlalchemy