# may need to delete databse or table sometime
DROP DATABASE IF EXISTS hw2;
DROP TABLE Schedule;

# housekeeping code
SHOW FULL COLUMNS FROM Bus;

CREATE DATABASE hw2;
USE hw2;

# the schema of the tables
CREATE TABLE City(
	name VARCHAR(20),
    state CHAR(2),
    PRIMARY KEY (name)
);

CREATE TABLE Bus(
	company VARCHAR(20),
    number INT UNIQUE,
    price INT,
    fromCity VARCHAR(20),
    toCity VARCHAR(20),
    PRIMARY KEY (company, number),
    FOREIGN KEY (fromCity) REFERENCES City(name),
    FOREIGN KEY (toCity) REFERENCES City(name)
    );
    
CREATE TABLE Schedule(
	company VARCHAR(20),
    bnum INT,
    departureTime DATETIME,
    arrivalTime DATETIME,
    PRIMARY KEY (company, bnum, departureTime),
    FOREIGN KEY (company) REFERENCES Bus(company),
    FOREIGN KEY (bnum) REFERENCES Bus(number)
);

# insert value into the table/relation

INSERT INTO City (name,state) VALUES ("New York", "NY");
INSERT INTO City (name,state) VALUES ("Cincinnati", "OH");
INSERT INTO City (name,state) VALUES ("Columbus", "OH");
INSERT INTO City (name,state) VALUES ("Chicago", "IL");

INSERT INTO Bus (company, number,price, fromCity, toCity) VALUES ("Whitedog", 13, 210, "New York", "Cincinnati");
INSERT INTO Bus (company, number,price, fromCity, toCity) VALUES ("Whitedog", 102, 56, "Columbus", "Cincinnati");
INSERT INTO Bus (company, number,price, fromCity, toCity) VALUES ("Whitedog", 2, 115, "Chicago", "Cincinnati");
INSERT INTO Bus (company, number, price, fromCity, toCity) VALUES ("Kbus",3,560,"Chicago","Columbus");

INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Whitedog",13,"2020-01-12 08:13","2020-01-12 17:20");
INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Whitedog",102,"2020-01-13 12:15","2020-01-13 13:30");
INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Kbus",3,"2020-01-12 10:30","2020-01-12 15:20");
INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Kbus",3,"2020-01-13 11:30","2020-01-12 16:20");






