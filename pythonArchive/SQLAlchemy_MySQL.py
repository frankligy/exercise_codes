
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String
engine = create_engine('mysql+pymysql://root:formysql1@localhost:3306/hw2',echo=True)

meta = MetaData()  # MetaData object

students = Table(   # Table object
   'students', meta,
   Column('id', Integer, primary_key = True),
   Column('name', String(50)),
   Column('lastname', String(50)),
)
meta.create_all(engine)

ins = students.insert().values(name='Frank')   # insert object
str(ins)
ins.compile().params

conn = engine.connect()
conn.execute(ins)  # can also do multiple insertion

'''multiple'''
conn.execute(students.insert(), [
   {'name':'Rajiv', 'lastname' : 'Khanna'},
   {'name':'Komal','lastname' : 'Bhandari'},
   {'name':'Abdul','lastname' : 'Sattar'},
   {'name':'Priya','lastname' : 'Rajhans'},
])

s = students.select().where(students.c.id>2)
str(s)
s.compile().params

result = conn.execute(s)  # ResultProxy
for row in result:
    print(row)
result.fetchone()  # generator, one after one
result.fetchall()  # a list with all tuple returned

# another way to initiate Select object
from sqlalchemy.sql import select
s = select([students])
result = conn.execute(s)

# textual SQL
from sqlalchemy.sql import text,bindparam
s = text("select students.name, students.lastname from students where students.name between :x and :y")
conn.execute(s, x = 'A', y = 'L').fetchall()

# anotehr way
s = s.bindparams(
   bindparam("x", type_= String),
   bindparam("y", type_= String)
)

result = conn.execute(s, {"x": "A", "y": "L"})

# alias, rename
from sqlalchemy.sql import alias, select
st = students.alias("a")
s = select([st]).where(st.c.id > 2)
conn.execute(s).fetchall()


# update and delete
stmt=students.update().where(students.c.lastname=='Khanna').values(lastname='Kapoor')
conn.execute(stmt)

stmt = students.delete().where(students.c.id > 2)
conn.execute(stmt)


'''Multiple table'''
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, ForeignKey
engine = create_engine('mysql+pymysql://root:formysql1@localhost:3306/hw2', echo=True)
meta = MetaData()

students = Table(
   'students', meta,
   Column('id', Integer, primary_key = True),
   Column('name', String(50)),
   Column('lastname', String(50)),
)

addresses = Table(
   'addresses', meta,
   Column('id', Integer, primary_key = True),
   Column('st_id', Integer, ForeignKey('students.id')),
   Column('postal_add', String(50)),
   Column('email_add', String(50)))

meta.create_all(engine)
conn = engine.connect()
conn.execute(students.insert(), [
   {'name':'Ravi', 'lastname':'Kapoor'},
   {'name':'Rajiv', 'lastname' : 'Khanna'},
   {'name':'Komal','lastname' : 'Bhandari'},
   {'name':'Abdul','lastname' : 'Sattar'},
   {'name':'Priya','lastname' : 'Rajhans'},
])

conn.execute(addresses.insert(), [
   {'st_id':1, 'postal_add':'Shivajinagar Pune', 'email_add':'ravi@gmail.com'},
   {'st_id':1, 'postal_add':'ChurchGate Mumbai', 'email_add':'kapoor@gmail.com'},
   {'st_id':3, 'postal_add':'Jubilee Hills Hyderabad', 'email_add':'komal@gmail.com'},
   {'st_id':5, 'postal_add':'MG Road Bangaluru', 'email_add':'as@yahoo.com'},
   {'st_id':2, 'postal_add':'Cannought Place new Delhi', 'email_add':'admin@khanna.com'},
])

# Join
from sqlalchemy import join
from sqlalchemy.sql import select
j = students.join(addresses, students.c.id == addresses.c.st_id)
stmt = select([students]).select_from(j)
result = conn.execute(stmt)
result.fetchall()














