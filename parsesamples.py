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

resultList=list();
for obj in jsonsFile:
    resultList.append(str(obj.values()).replace("[","{").replace("]","}")+"\n")

result=resultList.join()    

print result     
file.close()