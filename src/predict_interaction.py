from tools import *
from multiprocessing import Pool ## required by each file
from functools import partial
import random

def main(para):
    if 'ListSize' not in para:
        para['ListSize'] = '5'
    if 'ThreadNum' not in para:
        para['ThreadNum'] = '1'
    if 'RandomSeed' not in para:
        para['RandomSeed'] = '2014'
    if 'ListFile' not in para:
        para['ListFile'] = para['DataPath']+'/set1_all_new.txt'
    if 'ModelName' not in para:
        para['ModelName'] = 'RF-bin'
    if 'ModelFile' not in para:
        para['ModelFile'] = para['ListFile']+'.fea.max.'+para['ModelName']
    if 'FeatureType' not in para:
        para['FeatureType'] = 'SaveResidue'
    if 'SolutionNum' not in para:
        para['SolutionNum'] = '10'
    if 'PredictCutoff' not in para:
        para['PredictCutoff'] = '0.218'

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
        para['ListFormat'] = 'a1/a2/p1/p2/pdb1/ch1/pdb2/ch2'
        with open(para['NewList'], 'w') as tempfile:
            tempfile.write('Hhp1\tTas3\tP49674\tO94687\t4HOK\tA\t3D1D\tA\n')
            #tempfile.write('Hhp1\tMoc3\t4HOK\tA\tMOC3_modbase\t \n')
            #tempfile.write('Hhp1\tPpc89\t4HOK\tA\tPPC89_Modbase\t \n')

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
    residue_value = map_pdb_residue(predfile, para2['MapFile'])
    with open(para['ExeFile']+'data.txt', 'w') as outfile:
        for g1g2, res, val in residue_value:
            outfile.write('%s\t%s\t-1\t%s\n'%('\t'.join(g1g2.split('=')), res, val))
    from evaluate_perform import read_residue_data, group_residue
    idx, val1, val2 = read_residue_data(para['ExeFile']+'data.txt')
    pp_val = group_residue(idx, val2) ## using predicted value

    listfile = open(para['NewList'], 'r')
    outfile = open('PredictedInterfaceResidues_Cutoff%s.txt'%para['PredictCutoff'], 'w')
    cc1 = 0; cc2 = 0
    for line in listfile:
        ele = line.split()
        p1 = ele[0]; p2 = ele[1]
        s1 = ele[2]; s2 = ele[3]
        if (p1,p2) not in pp_val:
            cc2 += 1
            continue
        outfile.write('\t'.join([p1,p2,s1,s2,'%s--%s--ZDOCK-'%(s1,s2)]))
        cc1 += 1
        res = pp_val[(p1,p2)]
        int1 = [r for r in res if r.split(':')[0]==p1 and res[r] >= float(para['PredictCutoff'])]
        int2 = [r for r in res if r.split(':')[0]==p2 and res[r] >= float(para['PredictCutoff'])]
        ord1 = sorted([int(r.split(':')[-1]) for r in int1])
        ord2 = sorted([int(r.split(':')[-1]) for r in int2])
        outfile.write('\t'+','.join([str(i) for i in ord1]))
        outfile.write('\t'+','.join([str(i) for i in ord2])+'\n')
    listfile.close()
    outfile.close()
    show(['Saved',cc1,'Failed',cc2])

if __name__ == "__main__": main_fun(main)
