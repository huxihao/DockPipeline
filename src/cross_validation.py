from tools import *
from multiprocessing import Pool ## required by each file

def train_model(datafile, model='LR', mfile=None):
    if mfile == None:
        mfile = datafile+'.'+model
    r_file = 'temp_script_train_'+model+'.r'
    f = open(r_file, 'w')
    f.write(r'''
data = read.table("'''+datafile+r'''", header=FALSE)
na_idx = apply(is.na(data), 1, any)
vdata = data[!na_idx,] ## remove rows containing NA
    ''')
    if model == 'LR':
        f.write(r'''
m = glm(y ~ ., data=data.frame(y=vdata[,2], x=vdata[,c(-1,-2)]))
write.table(t(m$coefficients), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE)
        ''')
    elif model == 'RF-bin':
        f.write(r'''
require('randomForest')
m = randomForest(x=vdata[,c(-1,-2)], y=factor(vdata[,2]>=1), sampsize=min(5000,nrow(vdata)))
write.table(t(m$importance), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE)
        ''')
    elif model == 'RF':
        f.write(r'''
require('randomForest')
m = randomForest(x=vdata[,c(-1,-2)], y=vdata[,2], sampsize=min(5000,nrow(vdata)))
write.table(t(m$importance), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE)
        ''')
    elif model == 'SVM':
        f.write(r'''
require('e1071')
m = svm(x=vdata[,c(-1,-2)], y=vdata[,2], kernel="linear")
write.table(m$SV, file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE)
        ''')
    elif model == 'Meta':
        f.write(r'''
require('randomForest')
require('e1071')
## 1 name, 2 label, 3 surface, 4 RCF, 5:14 ZDOCK, 15:20 seq

m1 = randomForest(x=vdata[,15:20], y=vdata[,2], sampsize=min(5000,nrow(vdata)))
p1 = predict(m1, vdata[,15:20])
write.table(t(m1$importance), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)

m2 = naiveBayes(x=data.frame(rcf=vdata[,4], zdock=vdata[,5], seq=p1), y=factor(vdata[,2]>=1))
write.table(t(m2$tables), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)

m = list(m1, m2)
        ''')
    elif model == 'Meta-bin':
        f.write(r'''
require('randomForest')
require('e1071')
## 1 name, 2 label, 3 surface, 4 RCF, 5:14 ZDOCK, 15:20 seq

m1 = glm(y ~ ., data=data.frame(y=factor(vdata[,2]>=1), x=vdata[,5:14]), family = "binomial")
p1 = predict(m1, data.frame(x=vdata[,5:14]), type='response')
write.table(t(m1$coefficients), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE)

m2 = randomForest(x=vdata[,15:20], y=factor(vdata[,2]>=1), sampsize=min(5000,nrow(vdata)))
p2 = predict(m2, vdata[,15:20], type='prob')[,2]
write.table(t(m2$importance), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)

m3 = naiveBayes(x=data.frame(sur=vdata[,3], rcf=vdata[,4], zdock=p1, seq=p2), y=factor(vdata[,2]>=1))
write.table(t(m3$tables), file="'''+mfile+r'''.log", sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)

m = list(m1, m2, m3)
        ''')
    else:
        raise ValueError('Unknown model name %s'%model)
    f.write(r'''
save(m, file="'''+mfile+r'''") 
    ''')
    f.close()
    os.system('Rscript %s'%r_file)
    if os.path.exists(mfile):
        os.remove(r_file)
    return mfile

def model_predict(datafile, model='LR', mfile='model_file', outfile=None):
    if outfile == None:
        outfile = datafile + '.out'
    r_file = 'temp_script_predict_'+model+'.r'
    f = open(r_file, 'w')
    f.write(r'''
data = read.table("'''+datafile+r'''", header=FALSE)
all = rep(NA, nrow(data))
na_idx = apply(is.na(data), 1, any)
vdata = data[!na_idx,]
load(file="'''+mfile+r'''")
    ''')
    if model == 'LR':
        f.write(r'''
pred = predict(m, data.frame(x=vdata[,c(-1,-2)]))
        ''')
    elif model == 'RF-bin':
        f.write(r'''
require('randomForest')
pred = predict(m, vdata[,c(-1,-2)], type='prob')[,2]
        ''')
    elif model == 'RF':
        f.write(r'''
require('randomForest')
pred = predict(m, vdata[,c(-1,-2)])
        ''')
    elif model == 'SVM':
        f.write(r'''
require('e1071')
pred = predict(m, vdata[,c(-1,-2)])
        ''')
    elif model == 'Meta':
        f.write(r'''
require('randomForest')
require('e1071')
m1 = m[[1]]
m2 = m[[2]]
## 1 name, 2 label, 3 surface, 4 RCF, 5:14 ZDOCK, 15:20 seq
p1 = predict(m1, vdata[,15:20])
pred = predict(m2, data.frame(rcf=vdata[,4], zdock=vdata[,5], seq=p1), type='raw')[,2]
        ''')
    elif model == 'Meta-bin':
        f.write(r'''
require('randomForest')
require('e1071')
m1 = m[[1]]
m2 = m[[2]]
m3 = m[[3]]
## 1 name, 2 label, 3 surface, 4 RCF, 5:14 ZDOCK, 15:20 seq
p1 = predict(m1, data.frame(x=vdata[,5:14]), type='response')
p2 = predict(m2, vdata[,15:20], type='prob')[,2]
pred = predict(m3, data.frame(sur=vdata[,3], rcf=vdata[,4], zdock=p1, seq=p2), type='raw')[,2]
        ''')
    else:
        raise ValueError('Unknown model name %s'%model)
    f.write(r'''
all[!na_idx] = pred
write.table(data.frame(data[,1], all), file="'''+outfile+r'''", sep='\t', 
            quote=FALSE, append=FALSE, row.names=FALSE, col.names=FALSE) 
    ''')
    f.close()
    os.system('Rscript %s'%r_file)
    if os.path.exists(outfile):
        os.remove(r_file)
    return outfile

def split_data(datafile, fold, seed=2014):
    '''
    datafile from select_structure.py
    Each line is:
        p1, p2, pdb1, ch1, pdb2, ch2, pdbch
    in tab-seperated format 
    '''
    import random
    random.seed(seed)
    upp = [] ## a unique list of protein pairs
    ipp = [] ## index of protein pairs in pdb pairs
    fdx = [] ## fold index for protein pairs
    findex = 0
    infile = open(datafile, 'r')
    data = []
    for line in infile:
        ele = line.split('\t')
        data.append(ele)
        p1, p2 = ele[:2]
        if (p1,p2) not in upp: ## update set
            upp.append((p1,p2))
            fdx.append(findex) ## within [0, fold)
            findex += 1
            if findex == fold: ## roll back
                findex = 0
        ipp.append(upp.index((p1,p2))) ## append index
    infile.close()
    ipp_back = {} ## map protein pair back to pdb pair
    for i in xrange(len(ipp)):
        j = ipp[i]
        if j in ipp_back:
            ipp_back[j].append(i)
        else:
            ipp_back[j] = [i]
    random.shuffle(fdx) ## random permutation
    print len(upp), 'protein pairs from',
    print len(set(zip(*upp)[0]+zip(*upp)[1])), 'proteins'
    print 'Fold indices:', fdx[:20], '...'
    flist = []
    for f in set(fdx):
        ## names
        train_name = 'cv%s_fold%s_train.txt'%(fold, f+1)
        test_name = 'cv%s_fold%s_test.txt'%(fold, f+1)
        ## get PDB list from test
        ban_pdb = set()
        for i in xrange(len(fdx)):
            if fdx[i] == f: ## test set
                for j in ipp_back[i]:
                    for pdbid in data[j][6].split(','):
                        ban_pdb.add(pdbid)
        ## save them
        trainfile = open(train_name, 'w')
        testfile = open(test_name, 'w')
        bancc1 = 0; bancc2 = 0
        for i in xrange(len(fdx)):
            if fdx[i] != f: ## train set
                save1 = 0
                for j in ipp_back[i]:
                    pdb1,ch1,pdb2,ch2 = data[j][2:6]
                    if pdb1 not in ban_pdb and pdb2 not in ban_pdb:
                        trainfile.write('\t'.join(data[j]))
                        save1 += 1
                bancc1 += int(save1==0)
            else: ## test set
                save2 = 0
                for j in ipp_back[i]:
                    pdb1,ch1,pdb2,ch2 = data[j][2:6]
                    if pdb1 not in ban_pdb and pdb2 not in ban_pdb:
                        testfile.write('\t'.join(data[j]))
                        save2 += 1
                bancc2 += int(save2==0)
        show(['Remove', bancc1, 'protein pairs in train and', bancc2, 'protein pairs in test'])
        trainfile.close()
        testfile.close()
        yield train_name, test_name

def map_pdb_residue(filename, listname, useidx=1):
    ''' Input file format:
            for each line:
                p1=p2,index,pdb:chain:pos   value others
        Output format:
            a list:
                p1=p2, residue, max(value)
    '''
    data = []
    pdblist = set()
    infile = open(filename, 'r')
    for line in infile:
        ele = line.split('\t')
        g1g2, pdbidx, pdbres = ele[0].split(',')
        has_na = False
        for v in ele:
            if v.lower().startswith('na'):
                has_na = True
        if has_na:
            continue
        value = ele[useidx].strip()
        data.append((g1g2, pdbres, float(value)))
        pdblist.add(pdbres.split(':')[0])
    infile.close()
    from map_pdb_res import pdblist_to_uniprot
    res_map = pdblist_to_uniprot(pdblist)
    sup_map = {} ## supplimentary residue map from the input list
    with open(listname, 'r') as tempfile:
        for line in tempfile:
            p,s,c = line.split('\t')
            sup_map[s] = p
    comb = {}
    for pp, res, val in data:
        if res in res_map:
            res = res_map[res]
        elif res.count(':') == 2:
            pdb, ch, pos = res.split(':')
            if pdb in sup_map:
                res = sup_map[pdb]+':'+pos
        if (pp, res) in comb and comb[(pp, res)] > val:
            continue ## no need to update if having a larger value
        comb[(pp, res)] = val
    output = []
    for pp, res in comb:
        output.append((pp, res, comb[(pp, res)]))
    output.sort()
    return output

def combine_pdb_residue(filename, listname=None, outname=None):
    ''' Combine the features vectors from prepare_feature functions
        by taking the maximum values for PDB residues mapped to the same
        resiude in a protein pair.
        Output is a file with the same format
    '''
    if listname == None:
        listname = filename.replace('.fea','.map')
    if outname == None:
        outname = filename + '.max'
    data = []
    pdblist = set()
    with open(filename, 'r') as tempfile:
        for line in tempfile:
            ele = line.split('\t')
            info = ele[0]
            pp, idx, res = info.split(',')
            pdblist.add(res.split(':')[0])
            vals = [float(val) for val in ele[1:]]
            data.append([pp, res, vals])
    from map_pdb_res import pdblist_to_uniprot
    res_map = pdblist_to_uniprot(pdblist)
    sup_map = {} ## supplimentary residue map from the input list
    with open(listname, 'r') as tempfile:
        for line in tempfile:
            p,s,c = line.split('\t')
            sup_map[s] = p
    comb = {}
    for pp, res, vals in data:
        if res in res_map:
            res = res_map[res]
        if res.count(':') == 2:
            pdb, ch, pos = res.split(':')
            if pdb in sup_map:
                res = sup_map[pdb]+':'+pos
        if (pp,res) in comb:
            maxv = []
            for i,j in zip(comb[(pp,res)], vals):
                if i == float('nan'):
                    maxv.append(j)
                elif j == float('nan'):
                    maxv.append(i)
                else:
                    maxv.append(max(i,j))
            comb[(pp,res)] = maxv
        else:
            comb[(pp,res)] = vals
    with open(outname, 'w') as tempfile:
        for pp,res in sorted(comb.keys()):
            tempfile.write('%s,0,%s'%(pp,res))
            for val in comb[(pp,res)]:
                tempfile.write('\t%s'%val)
            tempfile.write('\n')
    return outname

def add_residue_label(resfile, pp_int=None):
    ## update file
    data = []
    with open(resfile, 'r') as infile:
        for line in infile:
            ele = line.split('\t')
            data.append(ele)
    cc = 0
    with open(resfile, 'w') as oufile:
        for ele in data:
            p1p2, idx, res = ele[0].split(',')
            both = tuple(p1p2.split('='))
            if pp_int == None:
                new = [ele[0],'-1']+ele[1:]
            elif both in pp_int:
                intres = pp_int[both]
                if res in intres:
                    new = [ele[0],'1']+ele[1:]
                else:
                    new = [ele[0],'0']+ele[1:]
            oufile.write('\t'.join(new))
            cc += 1
    return cc

def add_more_info(filename, ref_ddi, outname=None):
    ''' Add additional columns to show protein information, such as domains
    '''
    if outname == None:
        outname = filename + '.add'
    from domain_map import load_maps
    domains = load_maps(['protein2domain'])['protein2domain']
    infile = open(filename, 'r')
    oufile = open(outname, 'w')
    lastp = ''; lastd = []
    lastP = ''; lastD = []
    for line in infile:
        ele = line.strip().split('\t')
        p1,p2 = ele[0].split(',')[0].split('=')
        res = ele[0].split(',')[-1].split(':')
        if len(res) != 2:
            continue
        pro, pos = res
        pos = int(pos)
        ## to reduce search cost in a sorted list
        if pro != lastp:
            lastp = pro
            lastd = domains.get(pro,[])
            ## update the partner
            if pro == p1:
                lastP = p2
                lastD = domains.get(p2,[])
            elif pro == p2:
                lastP = p1
                lastD = domains.get(p1,[])
            else: ## unmatched protein names
                continue
        info_lv = 0 ## 0: not domain, 
                    ## 1: within domain,
                    ## 2: interacting domain.
        for name, start, end in lastd:
            if start <= pos and pos <= end:
                info_lv = max(info_lv, 1)
                for namE, starT, enD in lastD:
                    for P1,S1,E1,P2,S2,E2 in ref_ddi.get((name, namE), []):
                        if lastp != P1 and lastP != P2: ## not self
                            info_lv = max(info_lv, 2)
        ele.append(str(info_lv))
        oufile.write('\t'.join(ele) + '\n')
    infile.close()
    oufile.close()
    return outname

def intersect_poslist(real, pred):
    ''' Return the intersection of two position lists '''
    m1 = {}; m2 = {}
    for idx, res, val in real:
        m1[(idx, res)] = val
    for idx, res, val in pred:
        m2[(idx, res)] = val
    r = []
    for idx, res in m1:
        if (idx, res) in m2:
            r.append((idx, res, m1[(idx,res)], m2[(idx,res)]))
        else:
            r.append((idx, res, m1[(idx,res)], -1))
    r.sort()
    print len(m1), '&', len(m2), '=', len(r)
    return r

def get_res_labels(para):
    real_value = []
    from select_structure import get_all_cocry
    pp, ints, pdbs = get_all_cocry(para['DataPath']+'/Intres_052914.txt')
    pairs = set()
    with open(para['ListFile'], 'r') as infile:
        for line in infile:
            p1,p2 = line.split('\t')[:2]
            pairs.add((p1,p2))
    show(['Target to predict interface residues among', len(pairs), 'protein pairs with'])
    pro_len = {}
    with open(para['DataPath']+'/Uniprot_longest_052914.txt', 'r') as tempfile:
        for line in tempfile:
            protein, length = line.split()
            pro_len[protein] = int(length)
    pp_int = {}
    for p1p2, int1int2 in zip(pp, ints):
        if p1p2 not in pairs:
            continue
        p1, p2 = p1p2
        int1, int2 = int1int2
        intres = set()
        for i in int1: intres.add(p1+':'+str(i))
        for i in int2: intres.add(p2+':'+str(i))
        pp_int[p1p2] = intres
        if p1 == p2:
            for i in xrange(1,pro_len[p1]+1):
                res = p1+':'+str(i)
                real_value.append((p1+'='+p2, res, int(res in intres)))
        else:
            for i in xrange(1,pro_len[p1]+1):
                res = p1+':'+str(i)
                real_value.append((p1+'='+p2, res, int(res in intres)))
            for i in xrange(1,pro_len[p2]+1):
                res = p2+':'+str(i)
                real_value.append((p1+'='+p2, res, int(res in intres)))
    show([len(real_value), 'labelled residues'])
    return real_value, pp_int

def main(para):
    if 'ListFile' not in para:
        para['ListFile'] = '../data/set1_cmp_all.txt' 
    if 'SplitFold' not in para:
        para['SplitFold'] = '10'
    if 'SolutionNum' not in para:
        para['SolutionNum'] = '10'
    if 'RandomSeed' not in para:
        para['RandomSeed'] = '2014'
    if 'ThreadNum' not in para:
        para['ThreadNum'] = '1'
    if 'FeatureType' not in para:
        para['FeatureType'] = 'SaveZDOCK'
    if 'ModelName' not in para:
        para['ModelName'] = 'RF-bin'
    train_list = [] ## unique pdb list for training
    real_value = [] ## from cocrystal
    pred_value = [] ## from learning model
    other_vals = {}
    if para['FeatureType'] == 'SaveResidue':
        other_vals = {2:[], 3:[], 4:[], 5:[]}
#    elif para['FeatureType'] == 'SaveZDOCK':
#        other_vals = {2:[], 3:[], 4:[]}
    elif para['FeatureType'] == 'SavePatchDock':
        other_vals = {2:[], 3:[], 4:[]}
    elif para['FeatureType'] == 'SaveSequence':
        other_vals = {2:[], 3:[], 4:[]}

    all_res, pp_int = get_res_labels(para)
    import prepare_feature
    for train, test in split_data(para['ListFile'], int(para['SplitFold']), int(para['RandomSeed'])):
        if para['SplitFold'] == '1':
            train = para['ListFile']
            test = para['ListFile']
        if False: ## compare DDI network
            from generate_hSIN import generate_ddi2, get_pdb_subset
            generate_ddi2(output_file = train + '.ddi',
                          pdb_subset = get_pdb_subset(train))
            generate_ddi2(output_file = test + '.ddi',
                          pdb_subset = get_pdb_subset(test))
            from domain_map import reduced_ddi
            train_ddi = reduced_ddi(train + '.ddi')
            test_ddi = reduced_ddi(test + '.ddi')
        ###############################################################
        ## Train
        show(train, False)
        para2 = para.copy() ## copy parameters
        para2['ListFile'] = train
        para2['ListSize'] = '-1'
        para2['ListFormat'] = 'p1/p2/pdb1/ch1/pdb2/ch2'
        prepare_feature.main(para2)
        resfile = combine_pdb_residue(para2['OutFile'])
        #resfile = add_more_info(resfile, train_ddi)
        show(add_residue_label(resfile, pp_int), False)
        mfile = train_model(resfile, model=para['ModelName'])

        ###############################################################
        ## Test
        show(test, False)
        para2['ListFile'] = test
        prepare_feature.main(para2)
        resfile = combine_pdb_residue(para2['OutFile'])
        #resfile = add_more_info(resfile, train_ddi)
        show(add_residue_label(resfile, None), False)
        outfile = model_predict(resfile, model=para['ModelName'], mfile=mfile)

        ## Save values
        if not os.path.exists(outfile):
            continue ## skip this fold
        pred_value += map_pdb_residue(outfile, para2['MapFile'])
        for idx in other_vals:
            values = other_vals[idx]
            values += map_pdb_residue(resfile, para2['MapFile'], useidx=idx)
            other_vals[idx] = values
        show()
        if para['SplitFold'] != '1':
            ## clean files
            os.system('rm %s*'%train)
            os.system('rm %s*'%test)
    ## Get the real labels of residues in predicted protein pairs
    pred_pp = set([pp for pp, res, val in pred_value])
    show('Performance based on %s protein pairs'%len(pred_pp))
    for pp, res, val in all_res:
        if pp in pred_pp:
            real_value.append((pp, res, val))
        
    ## Comparison between know_value and pred_value agaist real_value
    save_list = [('cocry', real_value), ('pred'+'-'+para['FeatureType']+'-'+para['ModelName'], pred_value)]
    for idx in other_vals:
        save_list.append((para['FeatureType']+'-idx%s'%idx, other_vals[idx]))
    for name, value in save_list:
        head = []; real = []; pred = [];
        for idx, res, val1, val2 in intersect_poslist(real_value, value):
            p1, p2 = idx.split('=')
            head.append([p1, p2, res])
            real.append(val1)
            pred.append(val2)
        with open('cv_cocry_'+name+'_'+para['RandomSeed']+'.txt', 'w') as tempfile:
            for _a, _b, _c in zip(head, real, pred):
                tempfile.write('%s\t%s\t%s\n'%('\t'.join(_a), _b, _c))
        ## evaluate
        show(name, False)
        area, px, py, pc = performance(real, pred, x='FPR', y='TPR')
        show(area, False)
        #show(); show('FPR'); show(px, True); 
        #show('TPR'); show(py, True);
        area, px, py, pc = performance(real, pred, x='TPR', y='PPV')
        show(area, False)
        #show(); show('Recall'); show(px, True); 
        #show('Precision'); show(py, True);
        #show('Cutoff'); show(pc, True);
        show()

if __name__ == '__main__': main_fun(main)

