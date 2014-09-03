from tools import *
from multiprocessing import Pool ## required by each file
from functools import partial
import random

def format_chains(id, pdb1='../data/modbase/mlh1.pdb', ch1=' ', pdb2='../data/modbase/pms1.pdb', ch2=' '):
    if os.path.exists('%s.pdb'%id):
        return str(id), 'A', 'B'
    warnings.filterwarnings("ignore")
    from Bio.PDB import PDBParser
    par = PDBParser()
    pdb1file = open(pdb1, 'r')
    pdb2file = open(pdb2, 'r')
    import gzip
    if pdb1.endswith('.gz'):
        pdb1file = gzip.open(pdb1, 'rb')
    if pdb2.endswith('.gz'):
        pdb2file = gzip.open(pdb2, 'rb')
    s1 = par.get_structure('pdb1', pdb1file)
    s2 = par.get_structure('pdb2', pdb2file)
    from use_dock import save_chain, combine_chains
    save_chain('chain1.pdb', s1, ch1, 'A')
    save_chain('chain2.pdb', s2, ch2, 'B')
    combine_chains(['chain1.pdb', 'chain2.pdb'], '%s.pdb'%id)
    os.remove('chain1.pdb')
    os.remove('chain2.pdb')
    return str(id), 'A', 'B'

def main(para):
    if 'ListSize' not in para:
        para['ListSize'] = '5'
    if 'ThreadNum' not in para:
        para['ThreadNum'] = '1'
    if 'RandomSeed' not in para:
        para['RandomSeed'] = '2014'
    if 'ListFile' not in para:
        para['ListFile'] = '../data/set1_cmp_all.txt'
    if 'ModelName' not in para:
        para['ModelName'] = 'RF-bin'
    if 'ModelFile' not in para:
        para['ModelFile'] = para['ListFile']+'.fea.max.'+para['ModelName']
    if 'FeatureType' not in para:
        para['FeatureType'] = 'SaveZDOCK'
    if 'SolutionNum' not in para:
        para['SolutionNum'] = '10'
    if 'PredictCutoff' not in para:
        para['PredictCutoff'] = '0.5'

    ## Step 1: Train a model
    if 'ModelName' in para and not os.path.exists(para['ModelFile']):
        import cross_validation
        para1 = para.copy()
        para1['SplitFold'] = '1'
        cross_validation.main(para1)
    else:
        print 'Load model file', para['ModelFile']

    ## Step 2: Dock and predict new
    if 'NewList' not in para:
        para['NewList'] = 'list_from_user.txt'
        para['ListFormat'] = 'p1/p2/pdb1/ch1/pdb2/ch2'
        with open(para['NewList'], 'w') as tempfile:
            tempfile.write('Hhp1\tTas3\t4HOK\tA\t3D1D\tA\n')
            #tempfile.write('Hhp1\tMoc3\t4HOK\tA\tMOC3_modbase\t \n')
            tempfile.write('Hhp1\tPpc89\t4HOK\tA\tPPC89_Modbase\t \n')

    ## docking them and generate features
    import prepare_feature
    para2 = para.copy()
    para2['ListFile'] = para['NewList']
    prepare_feature.main(para2)
    feature_file = para2['OutFile']

    ## predicted by machine learning
    from cross_validation import add_residue_label, model_predict
    add_residue_label(feature_file)
    predfile = model_predict(feature_file, model=para['ModelName'], mfile=para['ModelFile'])
    
    ## Step 3: Summary
    from cross_validation import map_pdb_residue
    residue_value = map_pdb_residue(predfile)
    with open(para['ExeFile']+'data.txt', 'w') as outfile:
        for g1g2, res, val in residue_value:
            outfile.write('%s\t%s\t-1\t%s\n'%('\t'.join(g1g2.split('=')), res, val))
    from evaluate_perform import read_residue_data, group_residue
    idx, val1, val2 = read_residue_data(para['ExeFile']+'data.txt')
    pp_val = group_residue(idx, val2) ## using predicted value

    listfile = open(para['NewList'], 'r')
    for line in listfile:
        ele = line.split()
        p1 = ele[0]; p2 = ele[1]
        s1 = ele[2]; s2 = ele[3]
        show([p1,p2,s1,s2])
        if (s1,s2) in pp_val:
            res = pp_val[(s1,s2)]
        elif (s2,s1) in pp_val:
            res = pp_val[(s2,s1)]
        else:
            show('Docking\tFailed\n')
            continue
        int1 = [r for r in res if r.split(':')[0]==p1 and res[r] >= float(para['PredictCutoff'])]
        int2 = [r for r in res if r.split(':')[0]==p2 and res[r] >= float(para['PredictCutoff'])]
        ord1 = sorted([int(r.split(':')[-1]) for r in int1])
        ord2 = sorted([int(r.split(':')[-1]) for r in int2])
        show(','.join([str(i) for i in ord1]))
        show(','.join([str(i) for i in ord2]))
        show()
    listfile.close()

if __name__ == "__main__": main_fun(main)
