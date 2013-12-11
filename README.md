GeneModule
==========

=====
Files
=====
setup.py
genstats.py
README.md

====================
Installing genstats
====================
python setup.py install

=================================
Using the methods in genstats.py
=================================
After installation, to use the package in the python code, type:
from genestats import genstats
p = genstats()   #initialise class variable
p.sq(4)         #call class method using class variable  

=======
Methods
=======
sort_dict(dictionary) 
#sort a dictionary, return tuple sorted by values by ascent

pop(item, list)
#pop an item from a list

sq(n)
#return n*n

sum(list)
#return sum of all items in the list, item has to be numeric

mean(list)
#return average values of all items in the list, item has to be numeric

median(list)
#return median values of all items in the list, item has to be numeric

mod(list)
#return mod values of all items in the list, item has to be numeric

maf(dictionary)
#calculate allele frequencies given genotype count, input is a dictionary with at least 2 genotypes, eg {'AA':10, 'AB':5}, key value must have length of 2 characters, such as AA, 00, 01, AB, et cetera. Return a tuple (minor allele, minor allele frequencies)

hwe(dictionary)
#return chi-square value of hwe calculation of a SNP, input is a dictionary with 3 genotypes as keys and genotype counts as values, eg {'AA': 10, 'AB': 5, 'BB': 1}, key value must have length of 2 characters, such as AA, 00, 01, AB, et cetera.


