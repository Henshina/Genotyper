import itertools, argparse
import string
from collections import Counter
from sys import argv
import pdb

fileref = open("sorted_table_2.txt","r"); fh = open("GENOS.txt", "w")
##Table and Nucleotide Lists
line_num = 0; header = []; table_part = []; qindexes = []; lister=[]; table= {}; ref_table=[]; old_table=[]; RB= {}; Mol={}; CG={}
##Outputer lists
new_table = []; Qpositions = []; bases = []; one = []; data = []; data2 = []; data3 = []
##GT_Num Lists
match = []; gt = []; patternum = 1

for line in fileref.readlines():
    if line_num == 0:        
        header = line.split()
        for i,k in enumerate(header):
            if "qbase:" in k:
                qindexes.append(i)
##counter               
        pattern = 'qbase:'; count = 0; flag = True; start = 0; y = str(header)
        while flag:
            a = y.find(pattern,start)
            if a == -1:
                flag = False
            else:
                count += 1
                start = a+1
    else:
        line_parts = line.split(); ref_bases = [line_parts[3]]

        old_line_parts = line.split()
        old_table.append(old_line_parts)
##Nucleotides
        def get_pos(nuc):
            i = 0
            if nuc in ['--', 'indel', 'NA']:
                return 'ERR'
            for i,k in enumerate(ref_bases):
                if k == nuc:
                    return i
            ref_bases.append(nuc)
            return i+1
        
        for i in qindexes:
            line_parts[i] = str(get_pos(line_parts[i]))
        line_parts[3] = '0'; s = count + 4
        sorted(line_parts[3:s])
        table_part.append(line_parts)#converted table
##Table count
        obj = "".join(line_parts[3:s])
        if table.has_key(obj):
            table[obj] += 1
        else:
            table[obj] = 1
##Nucleotide Key
        nkey = ''.join(line_parts[3:s])
##Reference position Dictionary
        if RB.has_key(nkey):
            RB[nkey].append(line_parts[1])
        else:
            RB[nkey] = [line_parts[1]]
##Molecule Dictionary
        if Mol.has_key(nkey):
            Mol[nkey].append(line_parts[0])
        else:
            Mol[nkey] = [line_parts[0]]
##Numerical GT
        groups = {}
        for i,k in enumerate(table.keys()):
            groups[k] = i
    line_num += 1
##Output
Qpositions = header[4:s]
def get_no(p, tp):
    for q in range(0, len(tp)):
        if p == tp[q]:
            data.append(str(Qpositions[q]))
    return(data)

def get_no2(p, tp):
    for q in range(0, len(tp)):
        if p == tp[q]:
            data2.append(str(Qpositions[q]))
    return(data2)

def get_no3(p, tp):
    for q in range(0, len(tp)):
        if p == tp[q]:
            data3.append(str(Qpositions[q]))
    return(data3)

for i,line in enumerate(table_part):

    for q in qindexes:
        bases.append(table_part[i][q])
        
    one = get_no('1', bases)
    two = get_no2('2', bases)
    three = get_no3('3', bases)
    
    ph = []
    nucleotide = ''.join(line[3:s])
    cnt = table[nucleotide]
    g = str(groups[nucleotide])
    r = RB[nucleotide]
    M = Mol[nucleotide]
    ph = ["Group: %s" % '_'.join(g), "C." + str(cnt), "ref:%s" % '_'.join(r), "Molecule:%s" % '_'.join(M),"1: %s" % '_'.join(one), "2: %s" % '_'.join(two), "3: %s" % '_'.join(three)]
    nl = line[0:s] + ph + line[s+1:]
    new_table.append(nl)
    bases=[]
    data=[]
    data2=[]

#Writer
bloom = ["Count", "Reference Positions", "Molecule Number", "Headers With 1", "Headers With 2", "Headers With 3"]
hd=header
new_header = hd[:s] + bloom + hd[s:]
nh = str(new_header)

fh.write(nh + '\n')#prints out the header and breaks to the next line
for i,l in enumerate(new_table):
    fh.write('\t'.join(l)+'\n')#prints out the wrapper for line_parts
fh.close()
fileref.close()
