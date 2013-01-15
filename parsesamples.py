#!/cygdrive/f/python27/python.exe
import json
import sys
from jsons import *

file = open(sys.argv[1])
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