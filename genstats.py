class genstats:
	#init
	def __init__(self):
		pass

	#sort a dictionary, return tuple sorted by values by ascent
	def sort_dict(self, d):
		return sorted(d.items(), key=lambda item: item[1])
	
	#pop an item from a list
	def pop(self, item, mylist):
		l = []
		for i in mylist:
			if i != item:
				l.append(i)
		return l	
	
	#return square value
	def sq(self, n):
		return n*n
		
	#return sum of list
	def sum(self, mylist):
		total = 0
		for item in mylist:
			total += item
		return total	
		
	#return average of a list)
	def mean(self, mylist):
		return float(sum(mylist))/len(mylist)
		
	#return median of a list
	def median(self, mylist):
		tmp = mylist
		tmp.sort()
		return tmp[len(tmp)/2]
		
	#return mod of a list
	def mod(self, mylist):
		d = {}
		for item in mylist:
			try:
				d[item] += 1
			except KeyError:
				d[item] = 1
		key, value = self.sort_dict(d)[-1]
		return key
	
	#calculate allele frequencies given genotype count, input is a dictionary with at least 2 genotypes, return a tuple (minor allele, minor allele frequencies)
	def maf(self, d):
		k = d.keys()
		if len(k) < 2:
			return 1
		het = ''
		total = 0
		for geno in k:
			geno = geno.replace('/','')
			total += d[geno]
			if geno[0] != geno[1]:
				het = geno
				k = self.pop(geno, k)
		if total > 0:
			if het != '':
				a1 = (2.0*d[k[0]] + d[het])/(2*total)
				a2 = 1 - a1
				if a1 <= a2:
					return k[0][0], a1
				else:
					return self.pop(k[0][0], list(het))[0], a2
			else:
				a1 = 1.0*d[k[0]] / total
				a2 = 1 - a1
				if a1 <= a2:
					return k[0][0], a1
				else:
					return k[1][0], a2
		else:
			return 'Err:total = 0'
		
	#return chi-square value of hwe calculation of a SNP, input is a dictionary with 3 genotypes as keys and genotype counts as values, eg {'AA': 10, 'AB': 5, 'BB': 1}, key value must have length of 2 characters, such as AA, 00, 01, AB, et cetera.
	def hwe(self, d):
		k = d.keys()
		if len(k) != 3:
			return 'Err:dictionary should have 3 keys'
		het = ''
		total = 0
		for geno in k:
			geno = geno.replace('/','')
			total += d[geno]
			if geno[0] != geno[1]:
				het = geno
				k = self.pop(geno, k)
		hom1 = k[0]
		hom2 = k[1]
		if total > 0:
			pHom1 = (2.0*d[hom1] + d[het])/(2*total)
			pHom2 = 1 - pHom1
			oAA, oAB, oBB = d[hom1], d[het], d[hom2]
			eAA, eAB, eBB = self.sq(pHom1)*total, 2*pHom1*pHom2*total, self.sq(pHom2)*total
			if eAA > 0 and eAB > 0 and eBB > 0:
				chisq = self.sum([1.0*self.sq(oAA - eAA)/eAA, 1.0*self.sq(oAB - eAB)/eAB, 1.0*self.sq(oBB - eBB)/eBB])
				return chisq
			else:
				return 'Err:no variation'
		else:
			return 'Err:total = 0'
		
	#grep sequence from fasta file that match the search id
	def grep_fasta(self, fid, fasta_file):
		f = open(fasta_file, 'r')
		match = 0
		sid = ''
		for line in f:
			line = line.strip()
			if match == 1:
				string = '%s\n%s' %(sid, line)
				match = 0
				return string
			if '>' in line:
				if fid in line:
					match = 1
					sid = line
		return ''
		
	
