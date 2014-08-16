from tools import *

def read_residue_data(filename):
    idx = []
    val1 = []
    val2 = []
    infile = open(filename, 'r')
    for line in infile:
        p1, p2, res, real, pred = line.split('\t')
        idx.append((p1,p2,res))
        val1.append(float(real) >= 1)
        val2.append(float(pred))
    infile.close()
    return idx, val1, val2

def group_residue(idx, val):
    pp = {}
    for (p1, p2, res), val in zip(idx, val):
        both = (p1,p2)
        m = {} ## new entry
        if both in pp:
            m = pp[both]
        m[res] = max(m.get(res,0), val)
        pp[both] = m ## put back
    print 'Find', len(pp), 'pairs of proteins.'
    return pp


def unique_add(map1, key, item):
    list1 = map1.get(key, [])
    if item not in list1:
        map1[key] = list1 + [item]

def protein_domain(map_file):
    infile = open(map_file, 'r')
    data = {}
    cc = 0
    infile.readline() ## header
    for line in infile:
        cc += 1
        protein, domain, start, stop = line.split()
        start = int(start)
        stop = int(stop)
        unique_add(data, protein, (domain, start, stop))
    print cc, 'entries from', map_file
    return data

def read_ddi_network(sin_pair_file):
    infile = open(sin_pair_file, 'r')
    infile.readline() ## header
    gene = {}
    protein = {}
    domain = {}
    cc = 0; cc_diff = 0
    for line in infile:
        cc += 1
        ele = line.split('\t')
        g1 = ele[0]
        g2 = ele[1]
        p1 = ele[2]
        p2 = ele[3]
        d1 = ele[4]
        d1_st = int(ele[5])
        d1_ed = int(ele[6])
        d2 = ele[7]
        d2_st = int(ele[8])
        d2_ed = int(ele[9])
        if p1 != p2:
            cc_diff += 1
        ## map gene to protein
        gene[g1] = p1
        gene[g2] = p2
        ## also protein to gene
        gene[p1] = g1
        gene[p2] = g2
        ## map protein to interacting ranges
        unique_add(protein, p1, (d1_st, d1_ed, p2, d2_st, d2_ed))
        unique_add(protein, p2, (d2_st, d2_ed, p1, d1_st, d1_ed))
        ## map range to domain
        domain[(p1, d1_st, d1_ed)] = d1
        domain[(p2, d2_st, d2_ed)] = d2
    print cc, 'entries from', sin_pair_file
    print '\t', cc_diff, 'entries for two different proteins'
    print '\t', cc-cc_diff, 'entries for the same protein'
    return protein

def residue_domain(p1p2, res, domains, network, combine='top5', byInteract=True):
    p1, p2 = p1p2

    if p1 not in domains or p2 not in domains:
        return [], [], []

    ## 4.a: get all domains in the two proteins
    p1d = [(p1, d_n, d_st, d_ed) for d_n, d_st, d_ed in domains[p1]]
    p2d = [(p2, d_n, d_st, d_ed) for d_n, d_st, d_ed in domains[p2]]

    ## 4.b: find contact domains based on predicted contact residues
    contact = []
    for protein, dname, p_start, p_end in p1d + p2d:
        stack = []
        for position in xrange(p_start, p_end+1):
            pos = '%s:%s'%(protein, position)
            if pos in res:
                stack.append(res[pos])
        stack.sort(reverse=True)
        if len(stack) < 0: #(p_end - p_start + 1): ## check coverage
            contact.append(None)
        else:
            if combine == 'sum':
                contact.append(sum(stack))
            elif combine == 'prod':
                contact.append(1-reduce(lambda x,y:(1-x)*(1-y), stack))
            elif combine == 'top5':
                if len(stack) < 5:
                    contact.append(None) ## coverage issue
                else:
                    contact.append(stack[5-1]) ## the 5-th value
            else:
                raise ValueError('Unknown combining method %s'%combine)
    ## all pairwise
    p1d_p2d = []
    p1d_p2d_dd = []
    for ii in xrange(len(p1d)):
        for jj in xrange(len(p2d)):
            if contact[ii] == None or contact[len(p1d)+jj] == None:
                p1d_p2d.append(None)
            else:
                p1d_p2d.append(min(contact[ii], contact[len(p1d)+jj]))
            p1d_p2d_dd.append((p1d[ii], p2d[jj]))

    ## 4.c: known contact domains in the protein network
    real_contact = [False for i in p1d+p2d] ## N1 + N2
    real_p1d_p2d = [False for i in p1d_p2d] ## N1 * N2
    for p1_st, p1_ed, P2, p2_st, p2_ed in network.get(p1,[]):
        if P2 != p2: continue ## must be the same protein pair
        for ii in xrange(len(p1d)):
            _p1, _d1, _p1_st, _p1_ed = p1d[ii]
            if not (p1_st == _p1_st and _p1_ed == p1_ed):
                continue
            for jj in xrange(len(p2d)):
                _p2, _d2, _p2_st, _p2_ed = p2d[jj]
                if not (p2_st == _p2_st and _p2_ed == p2_ed):
                    continue
                real_contact[ii] = True
                real_contact[len(p1d)+jj] = True
                real_p1d_p2d[ii*len(p2d) + jj] = True
                #show([p1d[ii], p2d[jj], contact[ii], contact[len(p1d)+jj]], True)
    #show([p1, p2, len(real_p1d_p2d), '\n'])
    real = []; pred = []; info = []
    if byInteract:
        compare = zip(real_p1d_p2d, p1d_p2d, p1d_p2d_dd)
    else:
        compare = zip(real_contact, contact, p1d+p2d)
    for r,p,dd in compare:
        if p != None:
            real.append(r)
            pred.append(p)
            info.append(dd)
            #if (r and p <5) or (not r and p >=5): show(dd); show([r, p], True)
    return real, pred, info

def comparison(para, outname='evaluate_raw.txt'):
    domains = protein_domain(para['DataPath']+'/9606_SPROT_domains_061714.txt')
    network = read_ddi_network(para['DataPath']+'/hSIN_cc_052914.txt')
    outfile = open(outname, 'w')
    map_name = {#'cocry':'Random', 
                'SaveResidue-idx2':'ZDOCK-RCF',
                'SaveResidue-idx3':'ZDOCK-Top1',
                #'SaveResidue-idx4':'ZDOCK-Top2',
                #'SaveResidue-idx5':'ZDOCK-Top3',
                #'SaveZDOCK-idx2':'ZDOCK-Top1',
                #'SaveZDOCK-idx3':'ZDOCK-Top2',
                #'SaveZDOCK-idx4':'ZDOCK-Top3',
                'pred-SaveZDOCK-RF-bin':'RF-ZDOCK',
                #'pred-SavePatchDock-RF-bin':'PatchDock RF',
                #'pred-SaveSequence-RF-bin':'Sequence RF',
                'pred-SaveResidue-RF-bin':'RF-All',
                #'pred-SaveResidue-Meta-bin':'All Meta',
                }
    show('Method name')
    show('All residues')
    show('Positive residues')
    show('Negative residues')
    show('Area under ROC')
    show('Area under PR')
    show('All domain pairs')
    show('Postive domain pairs')
    show('Negative domain pairs')
    show('Area under ROC')
    show('Area under PR')
    show('All domains')
    show('Postive domains')
    show('Negative domains')
    show('Area under ROC')
    show('Area under PR')
    show()
    for infile in sorted(os.listdir('.')):
        if 'FeatureType' in para:
            if not infile.startswith('cv_cocry_pred-%s'%para['FeatureType']):
                continue
        else:
            if not infile.startswith('cv_cocry_'):
                continue
        ele = infile.split('_')
        if ele[2] not in map_name:
            continue
        dname = map_name[ele[2]]
        show(dname)

        #######################################################################
        ## Evalute interface residues
        n, x, y = read_residue_data(infile)
        if dname == 'Random':
            import random
            random.shuffle(y)
        show(len(x))
        show(x.count(True))
        show(x.count(False))

        outfile.write('Interface Residues\t'+dname)
        for i in x: outfile.write('\t'+str(i))
        outfile.write('\nInterface Residues\t'+dname)
        for i in y: outfile.write('\t'+str(i))
        outfile.write('\n')

        auroc, FPR, TPR, pc = performance(x, y, x='FPR', y='TPR')
        show(auroc)
        auprc, TPR, PPV, pc = performance(x, y, x='TPR', y='PPV')
        show(auprc)
        with open(dname+'_InterfaceResidue.txt', 'w') as tempfile:
            for a in zip(PPV, TPR, FPR):
                tempfile.write('%s\t%s\t%s\n'%a)

        #######################################################################
        ## Evalute domain-domain interaction
        pp = group_residue(n, y)
        real = []; pred = []; info = []
        for p1p2 in pp:
            r,p,f = residue_domain(p1p2, pp[p1p2], domains, network)
            real += r
            pred += p
            info += f
        show(len(real))
        show(real.count(True))
        show(real.count(False))

        outfile.write('Interacting Domains\t'+dname)
        for i in real: outfile.write('\t'+str(i))
        outfile.write('\nInteracting Domains\t'+dname)
        for i in pred: outfile.write('\t'+str(i))
        outfile.write('\n')

        auroc, FPR, TPR, pc = performance(real, pred, x='FPR', y='TPR')
        show(auroc)
        auprc, TPR, PPV, pc = performance(real, pred, x='TPR', y='PPV')
        show(auprc)
        with open(dname+'_DomainPair.txt', 'w') as tempfile:
            for a in zip(PPV, TPR, FPR):
                tempfile.write('%s\t%s\t%s\n'%a)
        ####### Repeat using interface domains
        real = []; pred = []; info = []
        for p1p2 in pp:
            r,p,f = residue_domain(p1p2, pp[p1p2], domains, network, byInteract=False)
            real += r
            pred += p
            info += f
        show(len(real))
        show(real.count(True))
        show(real.count(False))
        auroc, FPR, TPR, pc = performance(real, pred, x='FPR', y='TPR')
        show(auroc)
        auprc, TPR, PPV, pc = performance(real, pred, x='TPR', y='PPV')
        show(auprc)
        with open(dname+'_InterfaceDomain.txt', 'w') as tempfile:
            for a in zip(PPV, TPR, FPR):
                tempfile.write('%s\t%s\t%s\n'%a)
        show()
    outfile.close()

def main(para):
    if 'RepeatNum' not in para:
        para['RepeatNum'] = '1'
    for i in xrange(2014, 2014 + int(para['RepeatNum'])):
        if 'FeatureType' in para:
            import cross_validation
            para1 = para.copy()
            para1['RandomSeed'] = str(i)
            cross_validation.main(para1)
        else:
            import cross_validation
            para1 = para.copy()
            para1['RandomSeed'] = str(i)
            para1['FeatureType'] = 'SaveZDOCK'
            cross_validation.main(para1)
            para1['FeatureType'] = 'SaveResidue'
            cross_validation.main(para1)
            #print 'Please add parameter "FeatureType" to run cross validation.'
    comparison(para)
    os.system('Rscript %s/evaluate_perform.R . evaluate_raw.txt "ROC Curves for Predicting Interface Residues"'%para['SrcPath'])
    os.system('mv evaluate_raw.txt.pdf evaluate_residue_roc.pdf')
    os.system('Rscript %s/evaluate_perform.R . evaluate_raw.txt "Precision Curves for Predicting Interface Residues"'%para['SrcPath'])
    os.system('mv evaluate_raw.txt.pdf evaluate_residue_pre.pdf')
    os.system('Rscript %s/evaluate_perform.R . evaluate_raw.txt "ROC Curves for Predicting Interacting Domains"'%para['SrcPath'])
    os.system('mv evaluate_raw.txt.pdf evaluate_domain_roc.pdf')
    os.system('Rscript %s/evaluate_perform.R . evaluate_raw.txt "Precision Curves for Predicting Interacting Domains"'%para['SrcPath'])
    os.system('mv evaluate_raw.txt.pdf evaluate_domain_pre.pdf')

if __name__ == '__main__': main_fun(main)

