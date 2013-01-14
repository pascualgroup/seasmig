#!/cygdrive/f/python27/python.exe
import json
from jsons import *

file = open('samples.jsons')
jsonsFile = JsonsFile(file)
# TODO: import header and data ....

result="";
for obj in jsonsFile:
    result+="{"    
    result+=str(obj['rates']).replace("[","{").replace("]","}")
    result+=","
    result+=str(obj['rates2']).replace("[","{").replace("]","}")
    result+="}\n"
print result     
file.close()