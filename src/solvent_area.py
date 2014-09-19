from tools import *
import warnings
import subprocess
from Bio.PDB import *
import Bio.PDB.NACCESS as nac

DEFINE_NACCESS_PATH = os.path.dirname(os.path.realpath(__file__))+'/../bin/naccess2.1.1/naccess'

def get_area(pdbfile, chainids=None, newids=None, naccess=None):
    if naccess == None:
        naccess = os.path.abspath(DEFINE_NACCESS_PATH)
    '''
    NACCESS BioPython Interface:
 84              naccess_rel_dict[(chain_id, res_id)] = { 
 85                  'res_name': res_name, 
 86                  'all_atoms_abs': float(line[16:22]), 
 87                  'all_atoms_rel': float(line[23:28]), 
 88                  'side_chain_abs': float(line[29:35]), 
 89                  'side_chain_rel': float(line[36:41]), 
 90                  'main_chain_abs': float(line[42:48]), 
 91                  'main_chain_rel': float(line[49:54]), 
 92                  'non_polar_abs': float(line[55:61]), 
 93                  'non_polar_rel': float(line[62:67]), 
 94                  'all_polar_abs': float(line[68:74]), 
 95                  'all_polar_rel': float(line[75:80]) } 
    '''
    warnings.filterwarnings("ignore")
    data = {}
    try:
        rsa, asa = nac.run_naccess(model=None, pdb_file=pdbfile, 
            naccess=naccess, temp_path=os.path.abspath('.'))
    except Exception as e:
        print 'Error:', e
        return data
    d = nac.process_rsa_data(rsa)
    for residue in d:
        chain, (a, pos, b) = residue
        all_abs = d[residue]['all_atoms_abs']
        all_rel = d[residue]['all_atoms_rel']
        if chainids != None and chain not in chainids: 
            continue
        if newids == None:
            data[(chain, pos)] = [all_abs, all_rel]
        elif chainids != None: ## map chain ids one by one
            data[(newids[chainids.find(chain)], pos)] = [all_abs, all_rel]
        else:
            raise ValueError('Please provide parameter: chainids')
    return data

def get_outer(seperate, cutoff=15):
    "return the indeces of residues that are on the outer surface"
    p = set()
    for r in seperate: ## according to rel-ASA
        if seperate[r][1] >= cutoff:
            p.add(r)
    return p

def get_interface(seperate, combine, cutoff=1):
    "return the indeces of residues that are on the interface"
    p = set()
    for r in seperate: ## according to abs-ASA
        if r not in combine:
            continue
        #if abs(seperate[r][0] - combine[r][0]) >= cutoff:
        if seperate[r][0] - combine[r][0] >= cutoff:
            p.add(r)
    return p

def get_dSASA(seperate, combine, cutoff=15):
    "return the indeces of residues that are on the interface"
    p = {}
    for r in seperate: ## according to abs-ASA
        if r not in combine:
            continue
        if seperate[r][1] < cutoff:
            p[r] = 0 #float('nan')
        else:
            p[r] = seperate[r][0] - combine[r][0]
    return p

def change_of_rel_asa(pdb, ch1='', ch2='', dock=None, file1=None, file2=None):
    if dock == None:
        from use_dock import DockTool
        dock = DockTool()
        ch1, ch2 = dock.prepare_data(pdb, ch1, pdb, ch2, file1=file1, file2=file2)
    area1 = get_area(dock.chain1, dock.chain1_id, ch1)
    area2 = get_area(dock.chain2, dock.chain2_id, ch2)
    areac = get_area(dock.combine, dock.chain1_id+dock.chain2_id, ch1+ch2)
    dock.clean_temp_path()
    data = []
    for p in area1:
        data.append([pdb, ch1, p, area1[p][1] - areac[p][1]])
    for p in area2:
        data.append([pdb, ch2, p, area2[p][1] - areac[p][1]])
    return data

def two_pdb_chains(pdb, ch1='', ch2='', dock=None, file1=None, file2=None,
                   surface=False, allres=False, cutoff1=15, cutoff2=2):
    if dock == None:
        from use_dock import DockTool
        dock = DockTool()
        ch1, ch2 = dock.prepare_data(pdb, ch1, pdb, ch2, file1=file1, file2=file2)
    area1 = get_area(dock.chain1, dock.chain1_id, ch1)
    area2 = get_area(dock.chain2, dock.chain2_id, ch2)
    areac = get_area(dock.combine, dock.chain1_id+dock.chain2_id, ch1+ch2)
    dock.clean_temp_path()
    res1 = area1.keys()
    res2 = area2.keys()
    outer1 = get_outer(area1, cutoff1)
    outer2 = get_outer(area2, cutoff1)
    inter1 = get_interface(area1, areac, cutoff2)
    inter2 = get_interface(area2, areac, cutoff2)
    if allres:
        return inter1 & outer1, inter2 & outer2, outer1, outer2, res1, res2
    if surface:
        return inter1 & outer1, inter2 & outer2, outer1, outer2
    return inter1 & outer1, inter2 & outer2

def ppi_intres(p1, p2, pdb, ch1, ch2):
    from map_pdb_res import pdb_to_uniprot
    pdb_map = pdb_to_uniprot(pdb)
    i1, i2 = two_pdb_chains(pdb, ch1, ch2)
    j1 = set()
    j2 = set()
    for ch,i in i1:
        pos = '%s:%s:%s'%(pdb, ch, i)
        if pos in pdb_map:
            pp1, j = pdb_map[pos].split(':')
            if pp1 == p1:
                j1.add((p1,j))
    for ch,i in i2:
        pos = '%s:%s:%s'%(pdb, ch, i)
        if pos in pdb_map:
            pp2, j = pdb_map[pos].split(':')
            if pp2 == p2:
                j2.add((pp2,j))
    return j1, j2

def format_set(s, t='comma'):
    a = sorted([int(i) for c,i in s])
    if t == 'list':
        return a
    elif t == 'string':
        return [str(i) for i in a]
    elif t == 'comma':
        return ','.join([str(i) for i in a])
    else:
        raise ValueError('Unknown format %s'%t)

def jaccard_index(set1, set2):
    u = len(set1 | set2)
    if u == 0: return 0
    else: return len(set1 & set2) / float(u)

def main(para):
    if 'P1' not in para:
        para['P1'] = 'P00519'
    if 'P2' not in para:
        para['P2'] = 'P46108'
    if 'PDB' not in para:
        para['PDB'] = '1JU5'
    if 'Ch1' not in para:
        para['Ch1'] = 'C'
    if 'Ch2' not in para:
        para['Ch2'] = 'A'
    i1, i2 = ppi_intres(para['P1'], para['P2'], 
                        para['PDB'], para['Ch1'], para['Ch2'])
    show(format_set(i1))
    show(format_set(i2))
    show()

if __name__ == "__main__": main_fun(main)
