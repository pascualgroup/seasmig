import json
from collections import OrderedDict

class JsonsFile:
	def __init__(self, file, object_pairs_hook=OrderedDict):
		self.file = file
		self.object_pairs_hook = object_pairs_hook
	
	def readline(self):
		line = self.file.readline()
		if line == '':
			return None
		comment_index = line.find('#')
		if comment_index != -1:
			line = line[:comment_index]
		return line.strip()
	
	def __iter__(self):
		return self
	
	def next(self):
		lines = list()
		
		# Seek start
		while True:
			line = self.readline()
			if line == None:
				raise StopIteration()
			elif line == '---':
				pass
			elif line == '':
				pass
			else:
				lines.append(line)
				break
		
		# Read lines until end
		while True:
			line = self.readline()
			if line == None:
				break
			elif line == '---':
				break
			else:
				lines.append(line)
		
		return json.loads('\n'.join(lines), object_pairs_hook=self.object_pairs_hook)