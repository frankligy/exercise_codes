from datetime import datetime
from app import db,login

from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin

from hashlib import md5



@login.user_loader  # pre-requisite for flask_login work well with Users, when call current_user, will invoke this
def load_user(id):
    return User.query.get(int(id))

followers = db.Table(
    'followers',
    db.Column('follower_id', db.Integer, db.ForeignKey('user.id')),
    db.Column('followed_id', db.Integer, db.ForeignKey('user.id'))
)

class User(UserMixin,db.Model):   # multiple inheritence, using UserMixin to inherit 4 required items for flask_login to work
    '''
    is_authenticated
    is_active
    is_anonymous
    get_id()
    '''
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    email = db.Column(db.String(120), index=True, unique=True)
    password_hash = db.Column(db.String(128))
    about_me = db.Column(db.String(140))
    last_seen = db.Column(db.DateTime, default=datetime.utcnow)
    posts = db.relationship('Post', backref='author', lazy='dynamic')  # first have to be model name, refer to 'many' side
    '''
    u.posts will return all posts written by this author
    post.author will return the author who writes this post

    '''

    followed = db.relationship(
        'User', secondary=followers,
        primaryjoin=(followers.c.follower_id == id),
        secondaryjoin=(followers.c.followed_id == id),
        backref=db.backref('followers', lazy='dynamic'), lazy='dynamic')  # dynamic means don't excecute unless explicitly stated

    def __repr__(self):
        return '<User {}>'.format(self.username)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)
    
    def avatar(self, size):
        digest = md5(self.email.lower().encode('utf-8')).hexdigest()
        return 'https://www.gravatar.com/avatar/{}?d=identicon&s={}'.format(digest, size)
    
    def follow(self, user):
        if not self.is_following(user):
            self.followed.append(user)

    def unfollow(self, user):
        if self.is_following(user):
            self.followed.remove(user)

    def is_following(self, user):
        return self.followed.filter(
            followers.c.followed_id == user.id).count() > 0

    def followed_posts(self):
        followed = Post.query.join(
            followers, (followers.c.followed_id == Post.user_id)).filter(
                followers.c.follower_id == self.id)
        own = Post.query.filter_by(user_id=self.id)
        return followed.union(own).order_by(Post.timestamp.desc())
    


class Post(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    body = db.Column(db.String(140))
    timestamp = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))   # have to use table name "users" instead of Users, if the model name is AddressAndPhone, the table name would be address_and_phone

    def __repr__(self):
        return '<Post {}>'.format(self.body) 


'''
Play with SQLalchemy:

from app import db
from app.models import User,Post

u = Users(username='John',email='john@example.com')
db.session.add(u)
db.session.commit()

u = User(username='susan', email='susan@example.com')
db.session.add(u)
db.session.commit()

users = User.query.all()  # return a list, each will be a User instance
for u in users:
    print(u.id,u.username)  # even though id,username looks like class variable, but instance is able to access

u = Users.query.get(1)   # qeury by primary key, which is automatically incremented when add instance

u = Users.query.get(1)
p = Post(body='my first post!',author=u)  # backref in Users confer author attribute to Post
db.session.add(p)
db.session.commit()

posts = u.posts.all()  # a list contain all Post instance written by this author
p = Post.query.get(1)
p.author.username   # p.author return the a Users instance and can access all its attribute

User.query.order_by(User.username.desc()).all()

users = Users.query.all()
for u in users:
    db.session.delete(u)
db.session.commit()


'''


