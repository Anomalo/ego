import glob
import urllib
import os
import zipfile,os.path
import shutil


def unzip(source_filename):
	'''unzips file
	Keyword arguments:
	source_filename -- filename to unzip
	'''
    with zipfile.ZipFile(source_filename) as zf:
        for member in zf.infolist():
            # Path traversal defense copied from
            # http://hg.python.org/cpython/file/tip/Lib/http/server.py#l789
            words = member.filename.split('/')
            path = source_filename.split('.zip')[0]
            for word in words[:-1]:
                drive, word = os.path.splitdrive(word)
                head, word = os.path.split(word)
                if word in (os.curdir, os.pardir, ''): continue
                path = os.path.join(path, word)
            zf.extract(member, path)
	return path

			
def getOnlineGMTs():
	'''
	makes a dictionary from http://www.go2msig.org/ with available .gmt (go files) 
	and respective URLS
	'''
	origin = 'http://www.go2msig.org/cgi-bin/prebuilt.cgi?'
	html = urllib.urlopen(origin).read().split('<br>')[3:-2]
	taxids = {}
	for line in html:
		name = line.split('>')[1].split('<')[0]
		url = line.split('"')[1]#.split('"')[1]
		url = 'http://www.go2msig.org/cgi-bin'+url[1:]
		lastURL = urllib.urlopen(url).read().split('a href="')[-1]
		gmtURL = lastURL.split('"')[0]
		gmtName= gmtURL.split('/')[-1].replace('.zip','')
		taxids.update({gmtName:gmtURL})
	return taxids


def getGMT(taxon='Mus musculus',purge = False):
	'''
	Checks available .gmt files online and downloads more if available
	Keyword arguments:
	taxon -- taxon of interest
	purge -- set True to remove all existing .gmt files and redownload them
	'''
	prefix = 'ego/annotations/'

	gmtsAvail = getOnlineGMTs()
	if purge: os.system('rm -rf '+prefix+'*.gmt')
	
	gmtsDownloaded = glob.glob(prefix+'*.gtm')
	for gmt in gmtsDownloaded:
		name = gmt.split('/')[-1]
		if name in gmtsAvail:
			del gmtsAvail[name]
	for name, url in gmtsAvail.items():
		taxon = taxon.upper().replace(' ','_')
		if not taxon in name.upper(): continue
		name = prefix + name+'.zip'
		zipName = name
		urllib.urlretrieve(url,name)
		name = unzip(name)
		ultimate_name = name
		moved_name=name+'.temp'
		name = glob.glob(name+'/*.gmt')[0]
		print name

		print moved_name
		os.rename(name,moved_name)
		shutil.rmtree(ultimate_name)
		os.remove(zipName)
		os.rename(moved_name,ultimate_name)

		'''
		os.system(' '.join(['unzip',
				    name,
				    '-d',
				    prefix]))
		os.system('rm '+name+'.zip')
		'''
#updateGMTs(purge = True)


class GO:
	def __init__(self,gmt='default'):
		'''
		GO object to analyse go data
		Keyword arguments:
		gmt -- gmt file to use (default will find a gmt in ego/annotations)
		'''
		if gmt == 'default':
			gmt = glob.glob('ego/annotations/*.gmt')[0]
		self.gmt=gmt
		self.GOfile=gmt
		self.readGO()
		self.mkGOgenes()
		self.mkGOfreq()
		self.mkGOdef()
		self.mkGenesGO()
		#GOgenes = {} #{gene:[goIDs]...}
		#GOfreq = {} #{GOID: frequency...}
		#genesGO = {'GOlabel':[genes]}
		#labels_ID = {'GOlabel':'ID'}
	def mkGenesGO(self):
		'''makes a self.GenesGO dictionary of {GOname:[genes]}'''
		genesGO={}
		labels_ID={}
		for line in self.GOTXT:
			items = line.split()
			if items==[]:continue
			GOlabel = items.pop(0)
			url =  items.pop(0)
			ID = url.split('=')[-1]
			labels_ID[GOlabel]=ID
			genesGO[GOlabel]=[]
			for gene in items:
				gene = gene.upper()			
				genesGO[GOlabel].append(gene)
		self.genesGO=genesGO
		self.labels_ID=labels_ID
	def GOlabel_ID(self,label):
		'''given a go label it returns the GO ID
		Keyword arguments:
		label -- go label 
		'''
		return self.labels_ID[label]
	def mkGOdef(self):
		'''makes a self.GOdef dictionary of {GOID:GO long name}'''
		GOdef = {}
		for line in self.GOTXT:
			items = line.split()
			if items == []:continue
			GOID = items[1].split('=')[-1]
			GOdefinition = items[0]
			GOdef.update({GOID:GOdefinition})
		self.GOdef = GOdef

	def mkGOfreq(self):
		'''makes a self.GOfreq dictionary of {GOID : 0.45 (frequency of that annotation)}'''
		GOfreq = {}
		geneNumber = len(self.GOgenes)
		for line in self.GOTXT:
			items = line.split()
			if items==[]:continue
			GOID = items[1].split('=')[-1]
			amount = len(items)-2
			frequency = amount/float(geneNumber)
			GOfreq.update({GOID:frequency})
		self.GOfreq = GOfreq

		
	
	def mkGOgenes(self):
		''' makes a self.GOgenes dictionary of {gene : [GOID, GOID ...]}'''
		GOgenes = {}
		for line in self.GOTXT:
			items = line.split()
			if items==[]:continue
			GOlabel = items.pop(0)
			GOID = items.pop(0).split('=')[-1]
			for gene in items:
				gene = gene.upper()
				if not gene in GOgenes:
					GOgenes.update({gene:[]})
				GOgenes[gene].append(GOID)
		self.GOgenes = GOgenes
	def GOgeneNames(self, gene):
		'''returns a list of go names for that gene'''
		goids = self.GOgenes[gene.upper()]
		gonames = []
		for id in goids:
			gonames.append(self.GOdef[id])
		return gonames
	def GenesWithGO(self, GO):
		genes = self.genesGO[GO]
		return genes

	def readGO(self):
		''' makes a self.GOTXT wich is a list of all the lines of the .gmt file'''
		f = open(self.GOfile,'r')
		self.GOTXT = f.read().splitlines()
		f.close()
	
	def enrich(self, genes, label=True, P_cut_off = 0.05,v=True):
		'''
		Returns a dictionary of go annotaions present with their P value 
		
		Keyword arguments:
		genes -- a list of genes to enrich.
		label --  To annotate results with go label, leave false to annotate with GO id. (default True)
		P_cut_off -- the P value to filter with. (default 0.05)
		v -- verbose (default True)
		'''
		#raise NameError('not functional')
		from scipy import stats

		
		GOfound = []
		n = 0
		for gene in genes:
			if gene.upper() in self.GOgenes:
				GOfound += self.GOgenes[gene.upper()]
				n+=1
		if v:print 'Mapped',n,'from',len(genes),'given.'
		GOenriched = {}#dict.fromkeys(set(GOfound))
		for GOID in set(GOfound):
			x = GOfound.count(GOID)
			#n = len(genes)
			frequency = self.GOfreq[GOID]
			allgenes = len(self.GOgenes)
			background = frequency*allgenes			
			#print GOID
			P = stats.fisher_exact([[n,x],[allgenes,background]])	
			if label: GOID  = self.GOdef[GOID]
			GOenriched[GOID]={}
			GOenriched[GOID]['p'] = P[-1]
			GOenriched[GOID]['mapped'] = n
			GOenriched[GOID]['n'] = x

			#print GOID,self.GOdef[GOID],P
		keys = GOenriched.keys()
		for GOID in keys:
			if GOenriched[GOID]['p'] > P_cut_off:
				del GOenriched[GOID]
			
		return GOenriched


'''
go = GO()
def enrich(genes,key='GOID', P_cut_off = 0.05):
	return go.enrich(genes, key, P_cut_off)
'''

