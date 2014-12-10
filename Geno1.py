import itertools
import string
from collections import Counter
from sys import argv
import pdb
#script, filename = argv[1], filename2 = argv[2]
#fileref = open(filename, 'r')
#fh = open(filename2, 'w')

fileref = open("Phages.txt","r"); fh = open("GENOS.txt", "w")
#ids = [] #for the letter increment
line_num = 0; header = []; table_part = []; qindexes = []; lister=[]; table= {}; ref_table=[]; old_table=[]; RB= {}; Mol={}; CG={}

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
#        line_parts.sort(key=lambda x: x[2])
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
    line_num += 1
##Outputer
new_table = []; Qpositions = []; bases = []; one = []; data = []
Qpositions = header[4:s]

def get_no(p, tp):
    for q in range(0, len(tp)):
        if p == tp[q]:
            data.append(str(Qpositions[q]))
    return(data)

for i,line in enumerate(table_part):

    for q in qindexes:
        bases.append(table_part[i][q])
        
    one = get_no('1', bases)
    two = get_no('2', bases)
    three = get_no('3', bases)
    #print(one)
    
    ph = []
    nucleotide = ''.join(line[3:s])
    cnt = table[nucleotide]
    r = RB[nucleotide]
    M = Mol[nucleotide]
    
    ph = ["C." + str(cnt), "ref:%s" % '_'.join(r), "Molecule:%s" % '_'.join(M),"Headers w 1: %s" % '_'.join(one), "Headers w 2: %s" % '_'.join(two), "Headers w 3: %s" % '_'.join(three)]
    nl = line[0:s] + ph + line[s+1:]
    new_table.append(nl)
    bases=[]
#for line in new_table:
#    print('\t'.join(line))

##Writer
#print len(header)
fh.write(y + '\n')#prints out the header array and makes a break to the next line
for i,l in enumerate(new_table):
    fh.write('\t'.join(l)+'\n')#prints out the wrapper for line_parts
fh.close()
fileref.close()

#pdb.set_trace()
########Character Group not working
#        ordered_list = reduce(lambda x,y:x+y, map(lambda N:[''.join(x) for x in itertools.product(string.lowercase, repeat=N)], range(1,4)))#adding unicode
#        if CG.has_key(nkey):
#            CG[nkey].append(ordered_list[0])
#        else:
#            CG[nkey] = [ordered_list[0]]
