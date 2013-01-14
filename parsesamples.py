#!/cygdrive/f/python27/python.exe
import json
from jsons import *

file = open('samples.jsons')
jsonsFile = JsonsFile(file)

for obj in jsonsFile:
    print str(obj['rates']).replace("[","{").replace("]","}")
     
file.close()