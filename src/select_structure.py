from tools import *

def get_all_cry(filename='../data/pdbresiduemapping.txt.gz'):
    pdb_map = {}
    pro_map = {}
    import gzip
    infile = gzip.open(filename, 'r')
    infile.readline() ## header
    for line in infile:
        ele = line.strip().split()
        pdb, ch, pro = ele[:3] 
        ## chains in pdb
        pdb_chs = pdb_map.get(pdb, set())
        pdb_chs.add(ch)
        pdb_map[pdb] = pdb_chs
        ## chains for protein
        pro_chs = pro_map.get(pro, set())
        pro_chs.add((pdb,ch))
        pro_map[pro] = pro_chs
    infile.close()
    all_cry = {}
    for pro in pro_map:
        pdb_chs = set()
        for pdb,ch in pro_map[pro]:
            if len(pdb_map[pdb]) == 1: ## single protein chain
                pdb_chs.add((pdb,ch))
        if len(pdb_chs) > 0:
            all_cry[pro] = pdb_chs
    print 'Read in', len(all_cry), 'proteins with single crystal structures.'
    return all_cry

def get_cocry_id(filename='../data/excluded_pdbs.txt'):
    data = {}
    infile = open(filename, 'r')
    for line in infile:
        p1, p2, pdbs = line.split()
        data[(p1,p2)] = pdbs.strip().split(';')
        data[(p2,p1)] = data[(p1,p2)]
    infile.close()
    return data

def get_all_cocry(filename='../data/Intres_052914.txt'):
    pairs = [] ## protein pairs
    ints = [] ## interface residue
    pdbs = [] ## pdb chain pairs
    infile = open(filename, 'r')
    infile.readline() ## header
    for line in infile:
        ele = line.split()
        g1 = ele[0]
        g2 = ele[1]
        p1 = ele[3]
        p2 = ele[4]
        pairs.append((p1,p2))
        i1 = ele[9]
        i2 = ele[10]
        ints.append(([int(i) for i in i1.split(',')],
                     [int(i) for i in i2.split(',')]))
        pd = ele[11]
        pdbs.append([i.replace('-',':').split(':') for i in pd.split(';')])
    infile.close()
    print 'Read in', len(pairs), 'protein pairs in co-crystal structures.'
    return pairs, ints, pdbs

def get_all_mod(filename='../data/modbase_units.txt', sifts='../data/SIFTS'):
    data = {} ## all chains
    infile = open(filename, 'r')
    infile.readline() ## header
    for line in infile:
        ele = line.strip().split()
        protein = ele[0]
        pdb = ele[3].upper()
        begin = int(ele[4])
        end = int(ele[5])
        score = float(ele[8])
        chain = ele[11]
        one = data.get(protein, [])
        one.append((end-begin+1, score, pdb, chain))
        one.sort(reverse=True)
        data[protein] = one
        if True: ## save the residue map
            mapfile = sifts + '/'+chain+'_map.txt'
            if not os.path.exists(mapfile):
                print 'Save mapping file', mapfile
                with open(mapfile, 'w') as tempfile:
                    for i in xrange(begin, end+1):
                        tempfile.write('%s\t\t%s\t-\t%s\t%s\t-\n'%
                                       (chain, i, protein, i))
    infile.close()
    print 'Read in', len(data), 'proteins with model structures.'
    return data

def pdb_chain_info(filename='../data/pdb_chain_uniprot.tsv.gz'):
    import gzip
    infile = gzip.open(filename, 'r')
    infile.readline() ## time mark
    infile.readline() ## header
    data = {}
    for line in infile:
        pdb, ch, sp, res1, res2, pdb1, pdb2, sp1, sp2 = line.split()
        reg = data.get((pdb,ch), [])
        reg.append((sp, int(sp1), int(sp2)))
        data[(pdb,ch)] = reg
    infile.close()
    print 'Read in data for', len(data), 'pdb chains'
    return data

def both_from_other(old_pairs, old_pdbs, new_pairs, ban_ids=None):
    for p1p2 in new_pairs:
        if ban_ids == None or p1p2 not in ban_ids:
            in_test = []
        else:
            in_test = ban_ids[p1p2]
        get_p1 = set()
        get_p2 = set()
        for P1P2, PDBCH in zip(old_pairs, old_pdbs):
            if p1p2 == P1P2:
                continue ## itself
            p1,p2 = p1p2
            P1,P2 = P1P2
            if p1 == P1 or p1 == P2:
                for PDB, CH1, CH2 in PDBCH:
                    if PDB not in in_test:
                        if p1 == P1:
                            get_p1.add((PDB, CH1))
                        if p1 == P2:
                            get_p1.add((PDB, CH2))
            if p2 == P1 or p2 == P2:
                for PDB, CH1, CH2 in PDBCH:
                    if PDB not in in_test:
                        if p2 == P1:
                            get_p2.add((PDB, CH1))
                        if p2 == P2:
                            get_p2.add((PDB, CH2))
        yield (p1, p2, in_test, get_p1, get_p2)

def max_cover_pdb(protein, pdbch, top=-1):
    ''' Input a list of PDB chain pair for a protein
        Return a sorted list with the same coverage but less members
    '''
    pdbch = sorted(list(pdbch))
    from map_pdb_res import pdblist_to_uniprot
    pdbmap = pdblist_to_uniprot([pdb for pdb,ch in pdbch])
    newlist = []
    covres = set()
    while len(pdbch) > 0: ## all checked
        best_cov = covres
        best_pdbch = None
        for pdb, ch in pdbch:
            newres = set()
            for pdbres in pdbmap:
                if pdbres.startswith(pdb+':'+ch+':'):
                    prores = pdbmap[pdbres]
                    if prores.startswith(protein): ## same protein
                        newres.add(prores)
            if len(covres | newres) > len(best_cov):
                best_cov = covres | newres
                best_pdbch = pdb, ch
        if best_pdbch == None: ## maximum coverage
            break
        pdbch.remove(best_pdbch)
        if len(best_cov) > len(covres): ## improved
            newlist.append(best_pdbch)
            if len(newlist) == top:
                break
            covres = best_cov
    return newlist

def save_set1(para, max_case=1):
    if os.path.exists(para['Set1File']): return
    listfile = open(para['Set1File'], 'w')
    pp, ints, pdbs = get_all_cocry(para['DataPath']+'/Intres_052914.txt')
    pdb_ids = get_cocry_id(para['DataPath']+'/excluded_pdbs.txt')
    cc = 0
    for p1, p2, pdbs, get_p1, get_p2 in both_from_other(pp, pdbs, pp, pdb_ids):
        if len(get_p1) == 0 or len(get_p2) == 0:
            continue
        cc += 1
        info = ','.join(pdbs)
        for pdb1 in max_cover_pdb(p1, get_p1, max_case):
            for pdb2 in max_cover_pdb(p2, get_p2, max_case):
                listfile.write('%s\t%s\t'%(p1,p2))
                listfile.write('%s\t%s\t'%tuple(pdb1))
                listfile.write('%s\t%s\t'%tuple(pdb2))
                listfile.write(info+'\n')
    listfile.close()
    show(['Number of pairs from other co-crystal is', cc], True)

def save_set2(para, max_case=1):
    if os.path.exists(para['Set2File']): return
    listfile = open(para['Set2File'], 'w')
    pdb_cry = get_all_cry(para['DataPath']+'/pdbresiduemapping.txt.gz')
    pdb_ids = get_cocry_id(para['DataPath']+'/excluded_pdbs.txt')
    pp, ints, pdbs = get_all_cocry(para['DataPath']+'/Intres_052914.txt')
    cc = 0
    for p1,p2 in pp:
        if p1 not in pdb_cry or p2 not in pdb_cry:
            continue
        cc += 1
        info = ','.join(pdb_ids[(p1,p2)])
        for pdb1 in max_cover_pdb(p1, pdb_cry[p1], max_case):
            for pdb2 in max_cover_pdb(p2, pdb_cry[p2], max_case):
                listfile.write('%s\t%s\t'%(p1,p2))
                listfile.write('%s\t%s\t'%tuple(pdb1))
                listfile.write('%s\t%s\t'%tuple(pdb2))
                listfile.write(info+'\n')
    listfile.close()
    show(['Number of pairs from single crystal is', cc], True)

def save_set3(para, max_case=1):
    if os.path.exists(para['Set3File']): return
    listfile = open(para['Set3File'], 'w')
    pdb_mod = get_all_mod(para['DataPath']+'/modbase_units.txt')
    pp, ints, pdbs = get_all_cocry(para['DataPath']+'/Intres_052914.txt')
    pdb_ids = get_cocry_id(para['DataPath']+'/excluded_pdbs.txt')
    cc = 0
    for p1,p2 in pp:
        if p1 not in pdb_mod or p2 not in pdb_mod:
            continue
        cc += 1
        info = ','.join(pdb_ids[(p1,p2)])
        p1_mods = pdb_mod[p1]
        p2_mods = pdb_mod[p2]
        if max_case > 0 and len(p1_mods) > max_case:
            p1_mods = p1_mods[:max_case]
        if max_case > 0 and len(p2_mods) > max_case:
            p2_mods = p2_mods[:max_case]
        for pdb1 in p1_mods:
            for pdb2 in p2_mods:
                listfile.write('%s\t%s\t'%(p1,p2))
                listfile.write('%s\t%s\t'%tuple(pdb1[-2:]))
                listfile.write('%s\t%s\t'%tuple(pdb2[-2:]))
                listfile.write(info+'\n')
    listfile.close()
    show(['Number of pairs from model database is', cc], True)

def show_ratio(pp):
    i=0; j=0
    for p1,p2 in pp:
        if p1 == p2:
            i += 1
        else:
            j += 1
    show(['Ratio of Homo v.s. Hetero is', i/float(j)], True)

def show_overlap(para):
    s1 = set()
    s2 = set()
    s3 = set()
    with open(para['Set1File'],'r') as infile:
        for line in infile:
            s1.add(tuple(line.split('\t')[:2]))
    with open(para['Set2File'],'r') as infile:
        for line in infile:
            s2.add(tuple(line.split('\t')[:2]))
    with open(para['Set3File'],'r') as infile:
        for line in infile:
            s3.add(tuple(line.split('\t')[:2]))
    show('Set1'); show(len(s1),True)
    show('Set2'); show(len(s2),True)
    show('Set3'); show(len(s3),True)
    show('Set1&Set2'); show(len(s1 & s2),True)
    show('Set1&Set3'); show(len(s1 & s3),True)
    show('Set2&Set3'); show(len(s2 & s3),True)
    show('Set1&Set2&Set3'); show(len(s1 & s2 & s3),True)
    show('Set1|Set2'); show(len(s1 | s2),True)
    show('Set1|Set3'); show(len(s1 | s3),True)
    show('Set2|Set3'); show(len(s2 | s3),True)
    show('Set1|Set2|Set3'); show(len(s1 | s2 | s3),True)
    show('Set1'); show_ratio(s1)
    show('Set2'); show_ratio(s2)
    show('Set3'); show_ratio(s3)

def cover_by_others(fullset, lapset):
    otherset = set()
    with open(lapset, 'r') as tempfile:
        for line in tempfile:
            p1,p2 = line.strip().split()
            otherset.add((p1,p2))
            otherset.add((p2,p1))
    infile = open(fullset, 'r')
    newset = fullset.replace('.txt','')+'_new.txt'
    outfile = open(newset, 'w')
    savedset = set()
    for line in infile:
        p1,p2 = line.split('\t')[:2]
        if (p1,p2) in otherset:
            outfile.write(line)
            savedset.add((p1,p2))
            savedset.add((p2,p1))
    infile.close()
    outfile.close()
    print 'Number of pairs that are not saved', len(otherset-savedset)
    return newset

def main(para):
    if 'Set1File' not in para:
        para['Set1File'] = 'set1_all.txt'
    if 'Set2File' not in para:
        para['Set2File'] = 'set2_all.txt'
    if 'Set3File' not in para:
        para['Set3File'] = 'set3_top3.txt'
    save_set1(para, -1)
    save_set2(para, -1)
    save_set3(para, 3)
    show_overlap(para)
    para['Set1File'] = cover_by_others(para['Set1File'], para['DataPath']+'/set1.txt')
    para['Set2File'] = cover_by_others(para['Set2File'], para['DataPath']+'/set2.txt')
    para['Set3File'] = cover_by_others(para['Set3File'], para['DataPath']+'/set3.txt')
    show_overlap(para)
    para['Set1File'] = cover_by_others(para['Set1File'], para['DataPath']+'/set4.txt')
    para['Set2File'] = cover_by_others(para['Set2File'], para['DataPath']+'/set4.txt')
    para['Set3File'] = cover_by_others(para['Set3File'], para['DataPath']+'/set4.txt')
    show_overlap(para)

if __name__ == '__main__': main_fun(main)
