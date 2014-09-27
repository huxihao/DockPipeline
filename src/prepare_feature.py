from tools import *
from solvent_area import *
from multiprocessing import Pool, cpu_count
from functools import partial
import gzip
import warnings

def fun_map(p,f):
    ## Fake function for passing multiple arguments to Pool.map()
    return f(*p)

def zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num=5, pool='.'):
    from use_dock import UseZDOCK
    return dock_area(UseZDOCK(pool_path=pool), info, pdb1, ch1, pdb2, ch2, sol_num)

def zdock_clean(info, pdb1, ch1, pdb2, ch2, sol_num=5, pool='.'):
    from use_dock import UseZDOCK
    dock = UseZDOCK(pool_path=pool)
    infilename, linenum = dock_area(dock, info, pdb1, ch1, pdb2, ch2, sol_num, skip=True)[0]
    if linenum < sol_num:
        os.remove(infilename)
        dock.prepare_data(pdb1, ch1, pdb2, ch2)
        if os.path.exists(dock.solution) and dock.solution_number() < sol_num:
            os.remove(dock.solution)
            print 'Delete', dock.solution
        dock.clean_temp_path()
        show('The clean function has deleted %s'%infilename)
        return [[infilename, 'Deleted', linenum]]
    else:
        return []

def zdock_rcf(info, pdb1, ch1, pdb2, ch2, sol_num=None, pool='.'):
    from use_dock import UseZDOCK
    dock = UseZDOCK(pool_path=pool)
    outfile = '/rcf_%s_%s%s%s%s.txt'%(dock.TOOL.lower(),pdb1,ch1,pdb2,ch2)
    outfile = dock.POOL_PATH + outfile.replace(' ', '-')
    if os.path.exists(outfile) and len([l for l in open(outfile,'r')]) > 1:
        return [[outfile, 'Precomputed']]
    pdb1, ch1, pdb2, ch2 = dock.prepare_data(pdb1, ch1, pdb2, ch2)
    rcf_path = os.path.abspath(dock.TOOL_PATH + '/../proteins_journal')
    old_path = os.path.abspath('.')
    try:
        os.chdir(dock.TOOL_PATH)
        os.system('./mark_sur %s %s_m.pdb'%(dock._receptor, dock._receptor[:-4]))
        os.system('./mark_sur %s %s_m.pdb'%(dock._ligand, dock._ligand[:-4]))
        os.chdir(dock.TEMP_PATH)
        chain1 = os.path.basename(dock._receptor[:-4]+'_m.pdb')
        chain2 = os.path.basename(dock._ligand[:-4]+'_m.pdb')
        os.system("cp %s zout"%dock.solution)
        os.system('%s/source/reformat zout %s/source/euler.15 3.0'%(rcf_path, rcf_path))
        os.system('%s/source/contact_frequency -i zout.hout -o CG'%rcf_path)
        os.system('perl %s/source/Cal_RCF.pl %s.CG'%(rcf_path, chain1))
        os.system('perl %s/source/Cal_RCF.pl %s.CG'%(rcf_path, chain2))
        os.system('cat %s.CG.RCF %s.CG.RCF > %s_old'%(chain1, chain2, outfile))
    except Exception as e:
        print 'Failed to run RCF due to', e
        return [[outfile, 'Failed']]
    finally:
        dock.clean_temp_path()
        os.chdir(old_path)
    ## reformat file
    file1 = open(outfile + '_old', 'r')
    file2 = open(outfile, 'w')
    for line in file1:
        if line[0] in '-0123456789':
            ele = line.split('\t')
            pos, res, cha, val = ele[:4]
            val = float(val)
            if cha in dock.chain1_id:
                file2.write('%s\t%s\n'%(pdb1+':'+ch1.strip()+':'+pos, val))
            elif cha in dock.chain2_id:
                file2.write('%s\t%s\n'%(pdb2+':'+ch2.strip()+':'+pos, val))
            else:
                pass ## skip
    file1.close()
    file2.close()
    return [[outfile, 'New']]

def patchdock_area(info, pdb1, ch1, pdb2, ch2, sol_num=5, pool='.'):
    from use_dock import UsePatchDock
    return dock_area(UsePatchDock(pool_path=pool), info, pdb1, ch1, pdb2, ch2, sol_num)

def patchdock_clean(info, pdb1, ch1, pdb2, ch2, sol_num=5, pool='.'):
    from use_dock import UsePatchDock
    dock = UsePatchDock(pool_path=pool)
    infilename, linenum = dock_area(dock, info, pdb1, ch1, pdb2, ch2, sol_num, skip=True)[0]
    if linenum < sol_num:
        os.remove(infilename)
        dock.prepare_data(pdb1, ch1, pdb2, ch2)
        if os.path.exists(dock.solution) and dock.solution_number() < sol_num:
            os.remove(dock.solution)
            print 'Delete', dock.solution
        dock.clean_temp_path()
        show('The clean function has deleted %s'%infilename)
        return [[infilename, 'Deleted', linenum]]
    else:
        return []

def dock_area(dock, info, pdb1, ch1, pdb2, ch2, sol_num=5, skip=False):
    ''' Process DOCK output files and save solvent areas into a working file'''
    outfilename = '/solvent_%s_%s%s%s%s.txt.gz'%(dock.TOOL.lower(),pdb1,ch1,pdb2,ch2)
    outfilename = dock.POOL_PATH + outfilename.replace(' ', '-')
    if os.path.exists(outfilename): ## compressed by gzip
        line_num = -3 ## without headers
        tmpfile = gzip.open(outfilename, 'rb')
        for l in tmpfile:
            line_num += 1
        tmpfile.close()
        if line_num >= 0:
            print dock.TOOL, 'Precomputed\t', info, pdb1, ch1, pdb2, ch2, '\t', line_num
            if line_num < 2: show('May need to delete %s'%outfilename)
            return [[outfilename, line_num]]
    elif skip:
        return [[outfilename, sol_num]]
    ## Parse PDB and produce receptor and ligand files
    warnings.filterwarnings("ignore")
    pdb1, ch1, pdb2, ch2 = dock.prepare_data(pdb1, ch1, pdb2, ch2)

    ## solvent area change in real structure
    as1 = get_area(dock.chain1)
    as2 = get_area(dock.chain2)
    asc = get_area(dock.combine)
    miss = [-1,-1] ## define missing SASA

    po1 = sorted(dock.seq1.keys())
    po2 = sorted(dock.seq2.keys())

    ## save to file
    outfile = gzip.open(outfilename, 'wb')
    outfile.write('Residue')
    for c,p in po1: outfile.write('\t%s:%s:%s,%s'%(pdb1, dock.to_old_chain[c].strip(), p, dock.seq1[(c,p)]))
    for c,p in po2: outfile.write('\t%s:%s:%s,%s'%(pdb2, dock.to_old_chain[c].strip(), p, dock.seq2[(c,p)]))
    outfile.write('\nSeperate')
    for p in po1: outfile.write('\t%s,%s'%tuple(as1.get(p,miss)))
    for p in po2: outfile.write('\t%s,%s'%tuple(as2.get(p,miss)))
    outfile.write('\nCombined')
    for p in po1: outfile.write('\t%s,%s'%tuple(asc.get(p,miss)))
    for p in po2: outfile.write('\t%s,%s'%tuple(asc.get(p,miss)))
    outfile.close()

    cc = 0
    try:
        outline = ''
        dock.dock_them(sol_num)
        sfiles = dock.generate_complex(sol_num)
        score = dock.solution_scores()
        for i in xrange(min(sol_num, dock.solution_number())):
            cc += 1
            asp = get_area(sfiles[i])
            outline += '\n%.2f'%score[i]
            for p in po1: outline += '\t%s,%s'%tuple(asp.get(p,miss))
            for p in po2: outline += '\t%s,%s'%tuple(asp.get(p,miss))
        outfile = gzip.open(outfilename, 'ab') ## append
        outfile.write(outline) ## write all at once
        outfile.close()
    except Exception, v: ## any error will lead to fail
        show('Failed due to %s in:'%v)
        cc = 0 ## reset
    finally:
        dock.clean_temp_path()
    show([dock.TOOL, info, pdb1, ch1, pdb2, ch2, cc])
    return [[outfilename, cc]]

def feature_area(info, infilename, sol_num=100):
    ''' Save the change of accessibility area for each residue '''
    infile = gzip.open(infilename, 'rb')
    ## read header
    residue = infile.readline().strip().split('\t')[1:]
    seperate = infile.readline().strip().split('\t')[1:]
    combined = infile.readline().strip().split('\t')[1:]
    area = {}; area_c = {}
    for res, sep, com in zip(residue, seperate, combined):
        all_abs, all_rel = sep.split(',')
        area[res] = [float(all_abs), float(all_rel)]
        all_abs, all_rel = com.split(',')
        area_c[res] = [float(all_abs), float(all_rel)]
    from solvent_area import get_dSASA
    ## read docking results
    area_plist = []
    cc = 0
    scores = []
    for line in infile:
        if cc == sol_num: break
        cc += 1
        ele = line.split()
        scores.append(float(ele[0]))
        pred = ele[1:]
        area_p = {}
        for res, com in zip(residue, pred):
            all_abs, all_rel = com.split(',')
            area_p[res] = [float(all_abs), float(all_rel)]
        area_plist.append(area_p)
    infile.close()
    if cc < 2: return [] ## need to have at least two docking solutions
    ## encode
    norm_score = [score/max(scores) for score in scores]
    data_set = []
    for res in residue:
        pos, aa = res.split(',')
        entry = ['%s,%s'%(info, pos)]
        vals = [area[res][0] - area_p[res][0] for area_p in area_plist]
        entry += vals
        if len(vals) < sol_num:
            for i in xrange(sol_num - len(vals)):
                entry.append(0) ## fill in empty predictions
#        entry += scores
        ## save all
        data_set.append(entry)
    return data_set

def feature_residue(info, infilename):
    ''' Save the residue level features '''
    infile = gzip.open(infilename, 'rb')
    ## read header
    residue = infile.readline().strip().split('\t')[1:]
    seperate = infile.readline().strip().split('\t')[1:]
    combined = infile.readline().strip().split('\t')[1:]
    area = {}; area_c = {}
    for res, sep, com in zip(residue, seperate, combined):
        all_abs, all_rel = sep.split(',')
        area[res] = [float(all_abs), float(all_rel)]
        all_abs, all_rel = com.split(',')
        area_c[res] = [float(all_abs), float(all_rel)]
    infile.close()
    ## encode
    from residue_properties import get_vector
    data_set = []
    for res in residue:
        pos, aa = res.split(',')
        entry = ['%s,%s'%(info, pos), area[res][1]]
        ## features
        entry += get_vector(aa)
        ## save
        data_set.append(entry)
    return data_set

def feature_rcf(info, infilename, scorefile):
    ''' Save the residue level features '''
    infile = gzip.open(infilename, 'rb')
    ## read header
    residue = infile.readline().strip().split('\t')[1:]
    seperate = infile.readline().strip().split('\t')[1:]
    combined = infile.readline().strip().split('\t')[1:]
    area = {}; area_c = {}
    for res, sep, com in zip(residue, seperate, combined):
        all_abs, all_rel = sep.split(',')
        area[res] = [float(all_abs), float(all_rel)]
        all_abs, all_rel = com.split(',')
        area_c[res] = [float(all_abs), float(all_rel)]
    infile.close()
    from solvent_area import get_dSASA
    try:
        ## read data
        rcf_score = {}
        with open(scorefile, 'r') as tmpfile:
            for line in tmpfile:
                r, s = line.split('\t')
                rcf_score[r] = float(s)
        if len(rcf_score) < 1:
            return []
        ## save data
        data_set = []
        for res in residue:
            pos, aa = res.split(',')
            entry = ['%s,%s'%(info, pos), rcf_score[pos]]
            data_set.append(entry)
    except:
        show('Something wrong in %s'%scorefile)
        return []
    return data_set

def save_sequence(info, pdb1, ch1, pdb2, ch2, sol_num=5, dock_pool='.'):
    return feature_residue(info, zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0])

def save_zdock(info, pdb1, ch1, pdb2, ch2, sol_num=5, dock_pool='.'):
    'Binary ~ ZDOCK solution'
    return feature_area(info, zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0], sol_num)

def save_zdock_rcf(info, pdb1, ch1, pdb2, ch2, sol_num=5, dock_pool='.'):
    'Binary ~ ZDOCK solution'
    return feature_rcf(info, zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0],
                             zdock_rcf(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0])

def save_patchdock(info, pdb1, ch1, pdb2, ch2, sol_num=5, dock_pool='.'):
    'Binary ~ PatchDock solution'
    return feature_area(info, patchdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0], sol_num)

def save_all_features(info, pdb1, ch1, pdb2, ch2, sol_num=5, dock_pool='.'):
    f1 = feature_rcf(info, zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0],
                           zdock_rcf(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0])
    f2 = feature_area(info, zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0], sol_num)
    f3 = feature_residue(info, zdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0])
#    f4 = feature_area(info, patchdock_area(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0], sol_num)
    if f1 == [] or f2 == [] or f3 == []: return []
    return [a1+a2[1:]+a3[1:] for a1,a2,a3 in zip(f1,f2,f3)]

def save_all_files(info, pdb1, ch1, pdb2, ch2, sol_num=5, dock_pool='.'):
    from use_dock import UseZDOCK
    dock = UseZDOCK(pool_path=dock_pool)
    area_file = dock_area(dock, info, pdb1, ch1, pdb2, ch2, sol_num, skip=False)[0][0]
    rcf_file = zdock_rcf(info, pdb1, ch1, pdb2, ch2, sol_num, dock_pool)[0][0]
    f1 = feature_rcf(info, area_file, rcf_file)
    f2 = feature_area(info, area_file, sol_num)
    f3 = feature_residue(info, area_file)
    if f1 == [] or f2 == [] or f3 == []:
        dock.clean_temp_path()
        return []
    dock.prepare_data(pdb1, ch1, pdb2, ch2)
    if not os.path.exists('zdock_pool'):
        os.mkdir('zdock_pool')
    if not os.path.exists('zdock_pdbs'):
        os.mkdir('zdock_pdbs')
    if not os.path.exists('zdock_pool/'+dock.solution.split('/')[-1]): ## check one
        sfiles = dock.generate_complex(sol_num)
        from shutil import copy
        copy(dock.solution, 'zdock_pool/'+dock.solution.split('/')[-1])
        copy(area_file, 'zdock_pool/'+area_file.split('/')[-1])
        copy(rcf_file, 'zdock_pool/'+rcf_file.split('/')[-1])
        for i in xrange(sol_num):
            copy(sfiles[i], 'zdock_pdbs/%s-%s-%s-%s-ZDOCK-%d.pdb'%(pdb1,ch1,pdb2,ch2,i+1))
#            break ## save only one for now
    dock.clean_temp_path()
    return [a1+a2[1:]+a3[1:] for a1,a2,a3 in zip(f1,f2,f3)]

def save_data_set(filename, data_set):
    cc = 0
    with open(filename, 'w') as tmpfile:
        for data in data_set:
            for e in data:
                cc += 1
                tmpfile.write('\t'.join([str(v) for v in e])+'\n')
    #show('Save %s data points to %s\n'%(cc, os.path.abspath(filename)))

def main(para):
    ## default parameter
    if 'PoolPath' not in para:
        para['PoolPath'] = para['DataPath'] + '/dock_pool'
    if 'ListFile' not in para:
        para['ListFile'] = para['DataPath'] + '/pdb_list_example.txt'
    if 'SkipSize' not in para:
        para['SkipSize'] = '0'
    if 'ListSize' not in para:
        para['ListSize'] = '2'
    if 'ListFormat' not in para:
        para['ListFormat'] = 'p1/p2/pdb1,pdb2/ch1/ch2/re1/re2'
        #para['ListFormat'] = 'p1/p2/pdb1/ch1/pdb2/ch2/re1/re2'
    if 'ListType' not in para:
        para['ListType'] = 'Both'
    if 'MinSupport' not in para:
        para['MinSupport'] = '-1'
    if 'FeatureType' not in para:
        para['FeatureType'] = 'SaveZDOCK'
    if 'SolutionNum' not in para:
        para['SolutionNum'] = '10'
    if 'ThreadNum' not in para:
        para['ThreadNum'] = '1'
    if 'OutFile' not in para:
        para['OutFile'] = para['ListFile'] + '.fea'
    if 'MapFile' not in para:
        para['MapFile'] = para['ListFile'] + '.map'
    if not os.path.exists(para['PoolPath']):
        os.mkdir(para['PoolPath'])

    ## Task options
    fun = {'ZDOCK': zdock_area,
           'ZDOCK_clean': zdock_clean,
           'PatchDock': patchdock_area,
           'PatchDock_clean': patchdock_clean,
           'SaveSequence': save_sequence,
           'SaveZDOCK': save_zdock,
           'SaveRCF': save_zdock_rcf,
           'SavePatchDock': save_patchdock,
           'SaveResidue': save_all_features,
           'SaveResidueFile': save_all_files,
           }
    if para['FeatureType'] not in fun:
        print 'Try to use', para['FeatureType']
        raise ValueError('FeatureType must be within %s'%fun.keys())

    ## Generate task list
    listfile = open(para['ListFile'], 'r')
    ppdbfile = open(para['MapFile'], 'w')
    cc = -1
    par = []
    for line in listfile:
        cc += 1
        if cc < int(para['SkipSize']): continue
        if cc == int(para['ListSize']): break
        ele = para['ListFormat'].split('/')
        val = line.split('\t')
        val[-1] = val[-1].strip() ## remove end of line
        ch1 = ''
        ch2 = ''
        re1 = None
        re2 = None
        for ee, v in zip(ele, val):
            for e in ee.split(','):
                if e not in ['p1','p2','pdb1','pdb2','ch1','ch2','re1','re2']:
                    continue
                exec('%s="%s"'%(e,v))
        if re1 != None and len(re1.split(',')) < int(para['MinSupport']):
            print 'Skip', p1, p2, pdb1, ch1, pdb2, ch2
            continue
        if re2 != None and len(re2.split(',')) < int(para['MinSupport']):
            print 'Skip', p1, p2, pdb1, ch1, pdb2, ch2
            continue
        if para['ListType'] == 'Hom' and p1 != p2:
            continue
        if para['ListType'] == 'Het' and p1 == p2:
            continue
        ppdbfile.write('%s\t%s\t%s\n'%(p1, pdb1, ch1))
        ppdbfile.write('%s\t%s\t%s\n'%(p2, pdb2, ch2))
        par.append(('%s=%s,%s'%(p1,p2,cc), pdb1, ch1, pdb2, ch2, 
                    int(para['SolutionNum']), para['PoolPath']))
    listfile.close()
    ppdbfile.close()
    ## Process the tasks robustly
    data_set = []
    try:
        print 'There are', cpu_count(), 'CUPs on the machine.'
        if int(para['ThreadNum']) == 1:
            raise ## using the main thread
        elif int(para['ThreadNum']) < 1:
            pl = Pool(max(1, cpu_count()+int(para['ThreadNum'])))
        else:
            pl = Pool(int(para['ThreadNum']))
        print 'We use', pl._processes, 'cores in this run.'
        data_set = pl.map(partial(fun_map, f=fun[para['FeatureType']]), par, chunksize=1)
        #data_set = pl.map(partial(fun_map, f=fun[para['FeatureType']]), par)
        pl.close()
        pl.join() ## important to avoid zombie processes
    except: ## abandon multiple thread to check error places
        data_set = [fun[para['FeatureType']](*p) for p in par]
    ## Save all results
    if para['FeatureType'].find('_clean') >= 0:
        save_data_set(para['OutFile']+'.log', data_set)
    else:
        save_data_set(para['OutFile'], data_set)

if __name__ == "__main__": main_fun(main)
