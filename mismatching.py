 #!/usr/bin/python
import sys

(pattern,text,k) = (sys.argv[1],sys.argv[2],sys.argv[3])

class patternMatch(object):
		def __init__(self,k,pattern,text):
				self.pattern = pattern
				self.k = int(k)
				self.text = text
				self.subs = []
		def calMiss(self):
				for i in range(len(self.text)-len(self.pattern)+1):
						self.subt = self.text[i:i+len(self.pattern)]
						self.sum = 0
						for n,m in enumerate(self.subt):
								if self.subt[n] != self.pattern[n]:
										self.sum += 1
						if self.sum <= self.k:
								self.subs.append(i)
				print(self.subs)
pM1 = patternMatch(k,pattern,text)
pM1.calMiss()

