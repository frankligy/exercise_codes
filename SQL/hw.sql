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
INSERT INTO Bus (company, number, price, fromCity, toCity) VALUES ("Kbus",3,560,"Cincinnati","Chicago");

INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Whitedog",13,"2020-01-12 08:13","2020-01-12 17:20");
INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Whitedog",102,"2020-01-13 12:15","2020-01-13 13:30");
INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Kbus",3,"2020-01-12 10:30","2020-01-12 15:20");
INSERT INTO Schedule (company,bnum, departureTime,arrivalTime) VALUES ("Kbus",3,"2020-01-13 11:30","2020-01-12 16:20");


# let's start to finish the homework

# q1
SELECT DISTINCT b1.fromCity
FROM Bus b1 JOIN Bus b2
ON b1.fromCity = b2.toCity;

# q2
SELECT Bus.number, Bus.toCity,Schedule.arrivalTime
FROM Bus,Schedule
WHERE (Bus.company = "Whitedog") AND (Schedule.departureTime BETWEEN "2020-01-01" AND "2020-02-01");

# q3
SELECT MIN(Bus.price)
FROM Bus
WHERE (Bus.fromCity="New York" AND Bus.toCity="Cincinnati") OR (Bus.fromCity="Cincinnati" AND Bus.toCity="New York");

# q4
SELECT Count.company
FROM (
SELECT COUNT(Bus.number) AS c, Bus.company
FROM Bus
GROUP BY Bus.company) Count
WHERE Count.c > 5;

# q5
SELECT City.name
FROM City
WHERE NOT City.name in (
SELECT Bus.fromCity AS c FROM Bus
UNION
SELECT Bus.toCity AS c FROM Bus);

# q6
SELECT MAX(final.price)
FROM (
SELECT MIN(Bus.price) as price
FROM Bus
WHERE Bus.fromCity="Cincinnati" AND Bus.toCity="Houston"
UNION
SELECT (temp.fp + temp.sp) AS price
FROM (
SELECT b1.price AS fp, b1.fromCity AS departure, b2.price AS sp, b2.toCity AS destination
FROM Bus b1 JOIN Bus b2
ON b1.toCity = b2.fromCity) AS temp
WHERE temp.departure = "Cincinnati" AND temp.departure = 'Houston') AS final;









