#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:33:57 2020

@author: ligk2e
"""
#############basic requests, retrieve and query###########################

# seems liek requests module is more popular than urllib,urllib2,urllib3, it supports RESTful API (post,get,put,delete,etc)
import requests
import json

# retrieve(endpoint, this could be obtained by modifying URL)
response = requests.get("http://api.open-notify.org/astros.json")
status_code = response.status_code   # 200 is a int
response.content    #directly peek the format of data in raw bytes
response.text
response.headers
content = response.json()   # content will be a dictionary
content_json = json.dumps(content)  # serilize dictionary (python object) to JSON string(python string), all the keys will be coerced to stirng
content_back = json.loads(content_json)

# query(need guidance provided by the server)
parameters = {
        'lat': 39.10, 
        'lon': -84.51}  # coordinates of cincinnati, general notation, not degree-minute-second notation. W longitute means minus
response = requests.get("http://api.open-notify.org/iss-pass.json", params=parameters)
content = response.json()


# if return xml format
import xmltodict

my_xml = """
    <audience>
      <id what="attribute">123</id>
      <name>Shubham</name>
    </audience>
"""                    # XML is a string in python
# if xml in a file, using with open() as fd: xmltodict.parse(fd.read())

my_dict = xmltodict.parse(my_xml)  #{'audience':{'id':{'#text':'123','@what':'attribute'},'name':'Shubham'}}
# it is a OrderedDict type, All the keys and values are string.     
# tag > key; attribute > @attr; text > #text

my_json = json.dumps(my_dict)

xml_back = xmltodict.unparse(my_dict,pretty=True)  # use dict to convert back and must has a root key/tag


#################### Some real examples and endeavors#######################################

# An real example for accessing uniprot(retrieve example)

response = requests.get('https://www.uniprot.org/uniprot/P12345.xml')  #endpoint, if no .xml, it will be a web page
r_dict = xmltodict.parse(response.content)
r_json = json.dumps(r_dict)


# Another example for accessing uniprot, do ID mapping(query example)
parameters = {
        'from':'EMBL_ID',
        'to':'ACC+ID',
        'format':'tab',
        'query':'ENSG00000105699'}

response = requests.get('https://www.uniprot.org/uploadlists/',data=parameters) # other like auth, header 
response.text
r_dict = xmltodict.parse(response.content)
response.headers['content-type']


# uniprot official instruction, to my best knowledge, requests doesn't work, urllib work.
import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'ENSEMBL_ID',
'format': 'tab',
'query': 'P40925 P40926 O43175 Q9UM73 P97793'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8'))

a = response.decode('utf-8')
















