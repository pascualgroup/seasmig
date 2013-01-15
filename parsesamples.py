#!/cygdrive/f/python27/python.exe
import json
from jsons import *

file = open('samples.jsons')
jsonsFile = JsonsFile(file)

first = jsonsFile.next()
result="{"
for key in first.keys():
    result+=key
    if key!=first.keys()[-1]:
        result+=","
result+="}\n"

for obj in jsonsFile:
    result+=str(obj.values()).replace("[","{").replace("]","}")+"\n"    

print result     
file.close()