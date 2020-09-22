from flask import render_template,flash,redirect,url_for,request  # this function will take the html template and replace placeholder as python object, backended by jinja2
from app import app,db   # from app package/directory import the app instance
from app.forms import LoginForm,RegistrationForm,EditProfileForm,EmptyForm

from flask_login import current_user,login_user,logout_user,login_required
from app.models import User

from werkzeug.urls import url_parse



from datetime import datetime

@app.before_request
def before_request():   # execute before calling any view function
    if current_user.is_authenticated:
        current_user.last_seen = datetime.utcnow()
        db.session.commit()




@app.route('/')    # associate URL / to index function, meaning to say requesting http:localhost:5000/ will invoke function index
@app.route('/index')   # same for the above, associate /index to index function
@login_required  # if you visit index without being logged in, will redirect to login page(specified in __init__), and might trigger a flash automically?
def index():
    posts = [
        {'author':{'username':'Micheal'},'body':'Ok, you asked me to be the first post'},
        {'author':{'username':'Burnet'},'body':'Well, I will be the second'}
    ]
    return render_template('index.html',title='Home',posts=posts)


@app.route('/login',methods=['GET','POST'])   # you can use both get and post, default is only GET
def login():
    if current_user.is_authenticated:
        return redirect(url_for('index'))
    form = LoginForm()
    if form.validate_on_submit():   # if submit is click
        # check if using POST, then validator attached to each field
        user = User.query.filter_by(username=form.username.data).first()   # first is in contrast to all()
        if user is None or not user.check_password(form.password.data):
            flash('Invalid username or password')
            return redirect(url_for('login'))
        login_user(user,remember=form.remember_me.data)   # form.remember_me.data is a boolean variable
        next_page = request.args.get('next')   # 'next' string in URL
        if not next_page or url_parse(next_page).netloc != '':
            next_page = url_for('index')
        return redirect(next_page)
    return render_template('login.html',title='Sign In',form=form)


@app.route('/logout')
def logout():
    logout_user()
    return redirect(url_for('index'))   # will redirect to /index, then trigger @login_required decorator to bounce back to login


@app.route('/register', methods=['GET', 'POST'])
def register():
    if current_user.is_authenticated:
        return redirect(url_for('index'))
    form = RegistrationForm()
    if form.validate_on_submit():
        user = User(username=form.username.data, email=form.email.data)
        user.set_password(form.password.data)
        db.session.add(user)
        db.session.commit()
        flash('Congratulations, you are now a registered user!')
        return redirect(url_for('login'))
    return render_template('register.html', title='Register', form=form)


@app.route('/user/<username>')   # dynamic component, accept susan and then pass to function
@login_required
def user(username):
    user = User.query.filter_by(username=username).first_or_404()  # if not found, raise 404 exception
    posts = [
        {'author': user, 'body': 'Test post #1'},
        {'author': user, 'body': 'Test post #2'}
    ]
    form = EmptyForm()
    return render_template('user.html', user=user, posts=posts, form=form)




@app.route('/edit_profile', methods=['GET', 'POST'])
@login_required
def edit_profile():
    form = EditProfileForm()
    if form.validate_on_submit():
        current_user.username = form.username.data
        current_user.about_me = form.about_me.data
        db.session.commit()
        flash('Your changes have been saved.')
        return redirect(url_for('edit_profile'))
    elif request.method == 'GET':   # probably get to this page by pressing some link instead of pressing submit button
        form.username.data = current_user.username
        form.about_me.data = current_user.about_me
    return render_template('edit_profile.html', title='Edit Profile',
                           form=form)


@app.route('/follow/<username>', methods=['POST'])
@login_required
def follow(username):
    form = EmptyForm()
    if form.validate_on_submit():
        user = User.query.filter_by(username=username).first()
        if user is None:
            flash('User {} not found.'.format(username))
            return redirect(url_for('index'))
        if user == current_user:
            flash('You cannot follow yourself!')
            return redirect(url_for('user', username=username))  # if username are not referenced, it will be added as query string ?username=
        current_user.follow(user)
        db.session.commit()
        flash('You are following {}!'.format(username))
        return redirect(url_for('user', username=username))
    else:
        return redirect(url_for('index'))


@app.route('/unfollow/<username>', methods=['POST'])
@login_required
def unfollow(username):
    form = EmptyForm()
    if form.validate_on_submit():
        user = User.query.filter_by(username=username).first()
        if user is None:
            flash('User {} not found.'.format(username))
            return redirect(url_for('index'))
        if user == current_user:
            flash('You cannot unfollow yourself!')
            return redirect(url_for('user', username=username))
        current_user.unfollow(user)
        db.session.commit()
        flash('You are not following {}.'.format(username))
        return redirect(url_for('user', username=username))
    else:
        return redirect(url_for('index'))