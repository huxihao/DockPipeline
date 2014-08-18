
def case1():
    infile = open('../data/500_HBHQ_Docking_Candidates.txt', 'r')
    otfile = open('../work/500_HBHQ_Docking_Candidates.txt', 'w')
    for line in infile:
        ele = line.split('\t')
        pdbch1 = ele[2]
        pdbch2 = ele[3]
        if pdbch1.find('_') < 0:
            pdb1 = pdbch1
            ch1 = ' '
        else:
            pdb1, ch1 = pdbch1.split('_')
        if pdbch2.find('_') < 0:
            pdb2 = pdbch2
            ch2 = ' '
        else:
            pdb2, ch2 = pdbch2.split('_')
        otfile.write('\t'.join(ele[:2]+[pdb1,ch1,pdb2,ch2]+ele[4:]))
    infile.close()
    otfile.close()

case1()

