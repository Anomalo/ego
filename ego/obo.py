import urllib
import os


def getOBO(dir = 'ego/annotations/'):
		'''Retrieves an obo file from geneontology.org/ontology/go.obo
		Keyword arguments:
		dir -- directory where to save the obo file (default 'ego/annotations')
		'''
		url ='http://geneontology.org/ontology/go.obo'
		urllib.urlretrieve(url,dir+'go.obo')

class OBO:
		def __init__(self, fname='ego/annotations/go.obo'):
				''' initiates an OBO object
				Keyword arguments:
				fname -- location of an obo file (default is 'ego/annotations/go.obo')
				'''
				self.OBO=self.readOBO(fname)

		def findTopGOs(self):
				''' returns the root go annotations.
				in theory they should be 0003674 (molecular_function), 0005575 (cellular_component),
				and 0008150 (biological_process).
				'''
				out = {}
				for key, entry in self.OBO.iteritems():
						if 'is_obsolete' in entry:continue
						if not 'is_a' in entry:
								out[key]=entry
				return out

		def stepList(self, GO):
				''' returns the go ids between the GO provided and the root GO
				Keyword arguments:
				GO -- GO id 
				'''
				steps = []
				if not GO in self.OBO: return []
				while True:
						steps.append((GO,self.OBO[GO]['name']))
						if not 'is_a' in self.OBO[GO]:break
						GO = self.OBO[GO]['is_a']
						if type(GO)==list:GO=GO[-0]

				return steps
		def readOBO(self, fname):
				''' reads an obo file returning a dictionary of go ids with all the information 
				of the obo file.
				Keyword arguments:
				fname -- location of an obo file (default is 'ego/annotations/go.obo')
				'''
				f = open(fname,'r')
				obo = f.read()
				f.close()
				obo=obo.split('\n\n')
				OBO={}
				for block in obo:
						lines = block.splitlines()
						if len(lines)<2: continue
						term = lines.pop(0)
						if term != '[Term]':
								continue
						entry={}
						for line in lines:
								key, value = line.split(': ')[:2]
								if key == 'is_a':
										value = value.split()[0]
								if key in entry:
										try: entry[key].append(value)
										except: entry[key]=[entry[key],value]
								else:entry[key]=value
						OBO[entry['id']]=entry
				return OBO
		def levelUp(self, GO,steps=1):
				''' returns a tuple of the parent go id and the go label 
				Keyword arguments:
				GO -- GO id
				steps -- number of steps to go up the tree (default 1)
				'''
				ladder = self.stepList(GO)
				if len(ladder)<steps:return steps[0]
				else: return ladder[steps]

		def levelMax(self, GO,steps=1):
				''' returns a tuple of the top parent go id and the go label, 
				or however many steps before the top parent annotation
				Keyword arguments:
				GO -- GO id
				steps -- number of steps from the top of the tree (default 1)
				'''
				ladder= self.stepList(GO)
				if len(ladder)<steps:return steps[0]
				else: return ladder[-steps]

		def filterGMT(self,gmt,level=3,save=True,v=False):
				'''
				filters a gmt ensuring all go annotations are at the same level
				saves a gmtf file and returns the path to the new gmt
				reccomended to use:
						gmt = obo.filterGMT(gmt[args])
				Keyword arguments:
				gmt -- gmt file to use
				level -- Level to which modify GMTS (default 3)
				save -- wether to save results (default True)
				v -- verbose (default False)
				'''
				f = open(gmt)
				out = []
				outName = gmt.replace('.gmt','_'+str(level)+'.gmtf')
				while True:
						line = f.readline()
						if line =='':break
						id = line.split('=')[-1].split()[0]
						if len(self.stepList(id))==level:
								if v:
										print line.split()[0]
								out.append(line)


				f.close()
				if save:
						f = open(outName,'w')
						f.write('\n'.join(out))
						f.close()
						return outName

#getOBO()
#obo = 'annotations/go.obo'

o = OBO()
def stepList(GO):
	''' returns the go ids between the GO provided and the root GO
	Keyword arguments:
	GO -- GO id 
	'''
	return o.stepList(GO)

'''
for i in o.OBO:
	print i
	print o.stepList(i)
	break
'''
	
#gmt = glob.glob('annotations/*.gmt')[0]
#print o.filterGMT(gmt,level=2,v=True,save=False)

