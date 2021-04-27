DELETE from student
WHERE depart_name = "biology";

SELECT ALL student
FROM award;

# Join, you can either say FROM s,t (cross product) or say FROM s NATURAL JOIN t 
NATURAL JOIN
INNER Join
LEFT OUTER Join
RIGHT OUTER Join
FULL OUTER Join

# ALTER
ALTER table award ADD teacher,salary;
ALTER table award DROP student;

# UNION,INTERSECT,EXCEPT   (set semantics)
# UNION ALL, INTERSECT ALL, EXCEPT ALL  (bag semantics)

# ALL, ANY, SOME, EXISTS, NOT EXISTS, IN, NOT IN
# ALL will filter out items that do not say > all in (subquery)
# same explanation goes to ANY and SOME
# EXISTs will filter out items that do not exist in (subquery)

# HAVING, it is used after WHERE











