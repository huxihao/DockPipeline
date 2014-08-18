from tools import *
from Bio.PDB import *
import re
import warnings

DEFINE_PDB_PATH = '../data/pdb'
DEFINE_MOD_PATH = '../data/modbase'
DEFINE_ZDOCK_PATH = '../bin/zdock3.0.2_linux_x64'
DEFINE_PATCHDOCK_PATH = '../bin/PatchDock'
DEFINE_GRAMM_PATH = '../bin/gramm'

def get_pdb_file(pdbid, pdbfile=None, pdbpath=None, savepath=None):
    locations = []
    if pdbpath == None:
        pdbpath = os.path.abspath(DEFINE_PDB_PATH)
    if pdbfile == None:
        ## try the current path
        pdbfile = '%s.pdb'%pdbid
    if os.path.exists(pdbfile):
        return pdbfile
    locations.append(pdbfile)
    ## try the file in the data path
    pdbfile = '%s/%s.pdb'%(pdbpath, pdbid)
    if os.path.exists(pdbfile):
        return pdbfile
    locations.append(pdbfile)
    ## try the default path format in PDB
    pdbfile = pdbpath + '/%s/pdb%s.ent'%(pdbid[1:3].lower(), pdbid.lower())
    if os.path.exists(pdbfile+'.gz'): ## compressed file
        import gzip
        return gzip.open(pdbfile+'.gz', 'rb')
    locations.append(pdbfile+'.gz')
    ## try the given path
    if os.path.exists(pdbfile):
        return pdbfile
    locations.append(pdbfile)
    ## try simpler path format
    pdbfile = pdbpath + '/pdb%s.ent'%pdbid.lower()
    if os.path.exists(pdbfile):
        return pdbfile
    locations.append(pdbfile)
    ## try to download
    if savepath != None:
        os.system('wget http://www.rcsb.org/pdb/files/%s.pdb.gz'%pdbid.upper())
        pdbsavepath = '%s/%s'%(savepath, pdbid[1:3].lower())
        if not os.path.exists(pdbsavepath):
            os.makedirs(pdbsavepath)
        os.system('mv %s.pdb.gz %s/pdb%s.ent.gz'%
                  (pdbid.upper(), pdbsavepath, pdbid.lower()))
    else:
        raise IOError('Fail to file %s at:\n%s\n'%(pdbid,'\n'.join(locations)))
    return get_pdb_file(pdbid, pdbfile, savepath, savepath=None)

class DockTool(object):
    def __init__(self, temp_path='', pool_path='.'):
        self.POOL_PATH = os.path.abspath(pool_path)
        self.TEMP_PATH = temp_path
        self.TOOL = 'DOCK'

    def create_temp_path(self):
        while self.TEMP_PATH == '': ## find a non-used path name
            self.TEMP_PATH = 'TEMP_'+randomword(10)
            if os.path.exists(self.TEMP_PATH):
                self.TEMP_PATH = ''
        self.TEMP_PATH = os.path.abspath(self.TEMP_PATH) ## full path
        if not os.path.exists(self.TEMP_PATH):
            os.mkdir(self.TEMP_PATH)

    def __str__(self):
        return self.TOOL

    def read_pdb(self, pdbid, thefile=None, download=True):
        pdb_path = os.path.abspath(DEFINE_PDB_PATH)
        if len(pdbid) != 4:
            pdb_path = os.path.abspath(DEFINE_MOD_PATH)
        warnings.filterwarnings("ignore")
        if download:
            pdbfile = get_pdb_file(pdbid, thefile, pdb_path, pdb_path)
        else:
            pdbfile = get_pdb_file(pdbid, thefile, pdb_path, None)
        par = PDBParser()
        try:
            return par.get_structure(pdbid, pdbfile)
        except AssertionError:
            ## Ref: https://gist.github.com/JoaoRodrigues/9689449
            pdbfile = get_pdb_file(pdbid, thefile, pdb_path, None)
            if isinstance(pdbfile, basestring):
                infile = open(pdbfile, 'r')
            else:
                infile = pdbfile
            atom_lines = [l.strip() for l in infile if l.startswith(('ATOM', 'HETATM'))]
            sorted_atoms = sorted(atom_lines, key = lambda l: (l[0:6], l[21], int(l[22:26]), l[26], l[16]))
            infile.close()
            with open(pdbid.upper()+'.pdb', 'w') as newfile:
                newfile.write('\n'.join(sorted_atoms))
            return par.get_structure(pdbid, pdbid.upper()+'.pdb')

    def read_seq(self, pdbid, pdbfile=None, download=True):
        structure = self.read_pdb(pdbid, pdbfile, download)
        seq = {}
        for mod in structure:
            for chain in mod:
                ch = chain.id
                for residue in chain:
                    hetatm_flag, resseq, icode=residue.get_id()
                    if hetatm_flag != ' ' or icode != ' ':
                        continue
                    seq[(ch, resseq)] = residue.get_resname()
        return seq

    def set_solution(self, input1, input2):
        if input1 > input2:
            input1, input2 = input2, input1
        f1 = self.POOL_PATH+'/%s_%s%s'%(self.TOOL, input1, input2)
        self.solution  = f1.replace(' ', '-')

    def prepare_data(self, pdb1, chain1, pdb2, chain2, rotate=False, file1=None, file2=None, swap=True):
        ## working files
        self.create_temp_path()
        self.set_solution(pdb1+chain1, pdb2+chain2)
        self.combine = self.TEMP_PATH+'/combine.pdb'

        ## process chain ids with start and end information
        c1_st = None; c1_ed = None
        c2_st = None; c2_ed = None
        if ',' in chain1: chain1, c1_st, c1_ed = chain1.split(',')
        if ',' in chain2: chain2, c2_st, c2_ed = chain2.split(',')

        ## long chain means to use modbase
        if len(chain1) > 2:
            pdb1 = chain1
            chain1 = ' '
        if len(chain2) > 2:
            pdb2 = chain2
            chain2 = ' '

        s1 = self.read_pdb(pdb1, file1)
        s2 = self.read_pdb(pdb2, file2)

        ## empty chain means to use all chains
        if chain1 == '':
            for ch in s1[0]:
                chain1 += ch.id
        if chain2 == '':
            for ch in s2[0]:
                chain2 += ch.id

        ## map the chain ids to make them unique
        self.chain1_id = 'ACEGIKMOQSUWY'[:len(chain1)]
        self.chain2_id = 'BDFHJLNPRTVXZ'[:len(chain2)]
        if file1 != None and file2 != None: ## with existing files
            self.chain1 = self.TEMP_PATH+'/'+os.path.basename(file1)
            self.chain2 = self.TEMP_PATH+'/'+os.path.basename(file2)
        else: ## generate new files
            self.chain1 = self.TEMP_PATH+'/%s_%s.pdb'%(pdb1, chain1.replace(' ','-'))
            self.chain2 = self.TEMP_PATH+'/%s_%s.pdb'%(pdb2, chain2.replace(' ','-'))
            if self.chain1 == self.chain2: ## avoid overlap
                self.chain2 = self.TEMP_PATH+'/%s_%s1.pdb'%(pdb2, chain2.replace(' ','-'))

        ## save pdb files with new chain names
        save_chain(self.chain1, s1, chain1, self.chain1_id, c1_st, c1_ed)
        save_chain(self.chain2, s2, chain2, self.chain2_id, c2_st, c2_ed)
        ## read chain sequences
        self.seq1 = self.read_seq('0000', self.chain1)
        self.seq2 = self.read_seq('0000', self.chain2)
        ## save combined chains
        combine_chains([self.chain1, self.chain2], self.combine)
        if rotate:
            ## rotate the two chains to make the inital status non-trivial
            ## force rotate for the two chains from the same PDB file
            save_chain(self.chain1, self.read_pdb(pdb1), chain1, self.chain1_id, c1_st, c1_ed,
                    rotateSeed=int(str2num(pdb1[-4:]+chain1)))
            save_chain(self.chain2, self.read_pdb(pdb2), chain2, self.chain2_id, c2_st, c2_ed,
                    rotateSeed=int(str2num(pdb2[-4:]+chain2)))
            print ' Apply random rotation on', pdb1, chain1, 'and', pdb2, chain2

        ## If using PDB ids, make sure the longer chain is the receptor
        if (not swap) or len(self.seq1) >= len(self.seq2):
            self._receptor = self.chain1
            self._receptor_id = self.chain1_id
            self._ligand = self.chain2
            self._ligand_id = self.chain2_id
        else:
            self._receptor = self.chain2
            self._receptor_id = self.chain2_id
            self._ligand = self.chain1
            self._ligand_id = self.chain1_id
        self.to_old_chain = {}
        for old_id, new_id in zip(list(chain1+chain2), list(self.chain1_id+self.chain2_id)):
            self.to_old_chain[new_id] = old_id
        return pdb1, chain1, pdb2, chain2 ## re-assigned values

    def dock_them(self, N=5):
        if os.path.exists(self.solution) and (self.solution_number() >= N):
            return True ## use pre-computed solutions
        print 'Docking', os.path.basename(self.chain1), 'and', 
        print os.path.basename(self.chain2), 'with', self.TOOL, '...'
        return False

    def _fix_path_in_solution(self):
        ## fix path problem in solution file
        temp_path = re.compile('[^/ \t]*/') 
        lines = []
        with open(self.solution, 'r') as outfile:
            for line in outfile:
                lines.append(temp_path.sub('', line))
        with open(self.solution, 'w') as outfile:
            for line in lines: 
                outfile.write(line)

    def generate_complex(self):
        print 'Generate solutions'

    def solution_scores(self):
        print 'Retrieve the scores for all solutions'
        return [0]

    def solution_number(self):
        return len(self.solution_scores())

    def clean_temp_path(self):
        import shutil
        shutil.rmtree(self.TEMP_PATH, ignore_errors=True)

class UseZDOCK(DockTool):
    def __init__(self, temp_path='', pool_path=None, tool_path=None):
        if pool_path == None: pool_path = '.'
        DockTool.__init__(self, temp_path, pool_path)
        self.TOOL = 'ZDOCK'
        self.OLD_PATH = os.path.abspath('.')
        if tool_path == None:
            self.TOOL_PATH = os.path.abspath(DEFINE_ZDOCK_PATH)
        else:
            self.TOOL_PATH = os.path.abspath(tool_path)
        if not os.path.exists(self.TOOL_PATH):
            raise IOError('Please make sure ZDOCK path is %s'%self.TOOL_PATH)
        assert os.path.exists(self.TOOL_PATH+'/zdock'), 'Please install zdock'

    def dock_them(self, N=5):
        if os.path.exists(self.solution) and (self.solution_number() >= N):
            return True ## use pre-computed solutions
        print 'Docking', os.path.basename(self.chain1), 'and', 
        print os.path.basename(self.chain2), 'with', self.TOOL, '...'
        try:
            N = max(N, 3000)
            os.chdir(self.TOOL_PATH)
            tempsolu = self.TEMP_PATH + '/dock_solution.txt'
            os.system('./mark_sur %s %s_m.pdb'%(self._receptor, self._receptor[:-4]))
            os.system('./mark_sur %s %s_m.pdb'%(self._ligand, self._ligand[:-4]))
            os.system('./zdock -o %s -R %s_m.pdb -L %s_m.pdb -N %s -F &> %s/log.txt'%
                     (tempsolu, self._receptor[:-4], self._ligand[:-4], N, self.TEMP_PATH))
            #os.system('head -n %s %s > %s'%(N+5, tempsolu, self.solution))
            os.system('cp %s %s'%(tempsolu, self.solution))
            self._fix_path_in_solution()
            print 'Save result to', self.solution
            return True
        except:
            print 'Error in docking:'
            os.system('cat %s/log.txt'%self.TEMP_PATH)
            return False
        finally:
            os.chdir(self.OLD_PATH)

    def generate_complex(self, n=10):
        if n < 0: n = self.solution_number()
        os.chdir(self.TOOL_PATH)
        os.system('./mark_sur %s %s_m.pdb'%(self._receptor, self._receptor[:-4]))
        os.system('./mark_sur %s %s_m.pdb'%(self._ligand, self._ligand[:-4]))
        os.chdir(self.TEMP_PATH)
        try:
            os.system('cp %s/create_lig ./'%self.TOOL_PATH)
            os.system('%s/create.pl %s %s &> %s/create.out'%
                     (self.TOOL_PATH, self.solution, n, self.TEMP_PATH))
        except:
            print 'Error in generating solutions'
            if os.path.exists(self.TEMP_PATH+'/create.out'):
                with open(self.TEMP_PATH+'/create.out','r') as outfile:
                    for line in outfile:
                        print 'ZDOCK:', line.strip()
        os.chdir(self.OLD_PATH)
        return ['%s/complex.%s.pdb'%(self.TEMP_PATH, i+1) for i in range(n)]

    def solution_scores(self, skip=4):
        scores = []
        with open(self.solution, 'r') as tempfile:
            for i in xrange(skip):
                tempfile.readline()
            for l in tempfile:
                scores.append(float(l.split()[-1]))
        return scores

class UsePatchDock(DockTool):
    def __init__(self, temp_path='', tool_path=None, pool_path=None):
        if pool_path == None: pool_path = '.'
        DockTool.__init__(self, temp_path, pool_path)
        self.TOOL = 'PatchDock'
        self.OLD_PATH = os.path.abspath('.')
        if tool_path == None:
            self.TOOL_PATH = os.path.abspath(DEFINE_PATCHDOCK_PATH)
        else:
            self.TOOL_PATH = os.path.abspath(tool_path)
        if not os.path.exists(self.TOOL_PATH):
            raise IOError('Please make sure PatchDock path is %s'%self.TOOL_PATH)
        assert os.path.exists(self.TOOL_PATH+'/patch_dock.Linux'), 'Please install patch_dock'

    def dock_them(self, N=5):
        if os.path.exists(self.solution) and (self.solution_number() >= N):
            return True ## use pre-computed solutions
        print 'Docking', os.path.basename(self.chain1), 'and', 
        print os.path.basename(self.chain2), 'with', self.TOOL, '...'
        try:
            N = max(N, 5000)
            os.chdir(self.TEMP_PATH)
            tempsolu = self.TEMP_PATH + '/patchdock_solution.txt'
            os.system('%s/buildParams.pl %s %s &> log.txt'%(self.TOOL_PATH, self._receptor, self._ligand))
            os.system('%s/buildMS.pl %s %s &> log.txt'%(self.TOOL_PATH, self._receptor, self._ligand))
            os.system('%s/patch_dock.Linux params.txt %s &> log.txt'%(self.TOOL_PATH, tempsolu))
            os.system('head -n %s %s > %s'%(N+35, tempsolu, self.solution))
            self._fix_path_in_solution()
            print 'Save result to', self.solution
            return True
        except:
            print 'Error in docking:'
            os.system('cat log.txt patch_dock.log')
            return False
        finally:
            os.chdir(self.OLD_PATH)

    def generate_complex(self, n=10):
        if n < 0: n = self.solution_number()
        os.chdir(self.TEMP_PATH)
        try:
            os.system('cp %s/chem.lib ./'%self.TOOL_PATH)
            os.system('cp %s output.txt'%self.solution)
            os.system('%s/transOutput.pl output.txt 1 %s &> log.txt'%(self.TOOL_PATH, n))
        except:
            print 'Error in generating solutions'
            os.system('cat log.txt')
        finally:
            os.chdir(self.OLD_PATH)
        return ['%s/output.txt.%s.pdb'%(self.TEMP_PATH, i+1) for i in range(n)]

    def solution_scores(self, skip=35):
        scores = []
        with open(self.solution, 'r') as tempfile:
            for i in xrange(skip):
                tempfile.readline()
            for l in tempfile:
                scores.append(float(l.split('|')[1]))
        return scores

class UseGRAMM(DockTool):
    def __init__(self, temp_path='', tool_path=None, pool_path=None):
        if pool_path == None: pool_path = '.'
        DockTool.__init__(self, temp_path, pool_path)
        self.chain1_id = 'A' ## must be A and B
        self.chain2_id = 'B'
        self.TOOL = 'GRAMM'
        self.OLD_PATH = os.path.abspath('.')
        if tool_path == None:
            self.TOOL_PATH = os.path.abspath(DEFINE_GRAMM_PATH)
        else:
            self.TOOL_PATH = os.path.abspath(tool_path)
        if not os.path.exists(self.TOOL_PATH):
            raise IOError('Please make sure GRAMM path is %s'%self.TOOL_PATH)
        assert os.path.exists(self.TOOL_PATH+'/gramm'), 'Please install gramm'

    def dock_them(self, N=5):
        if os.path.exists(self.solution) and (self.solution_number() >= N):
            return True ## use pre-computed solutions
        print 'Docking', os.path.basename(self.chain1), 'and', 
        print os.path.basename(self.chain2), 'with', self.TOOL, '...'
        N = max(N, 1000)
        os.chdir(self.TEMP_PATH)
        os.putenv('GRAMMDAT', self.TOOL_PATH)
        tmpfile = open('rpar.gr', 'w')
        tmpfile.write(
'''Matching mode (generic/helix) ....................... mmode= generic
Grid step ............................................. eta= 6.8
Repulsion (attraction is always -1) .................... ro= 6.5
Attraction double range (fraction of single range) ..... fr= 0.
Potential range type (atom_radius, grid_step) ....... crang= grid_step
Projection (blackwhite, gray) ................ ....... ccti= gray
Representation (all, hydrophobic) .................... crep= all
Number of matches to output .......................... maxm= %s
Angle for rotations, deg (10,12,15,18,20,30, 0-no rot.)  ai= 20
            '''%N)
        tmpfile.close()
        tmpfile = open('rmol.gr', 'w')
        tmpfile.write(
'''# Filename  Fragment  ID      Filename  Fragment  ID     [paral/anti  max.ang]
# ----------------------------------------------------------------------------
%s %s ch1 %s %s ch2
            '''%(self.chain1.split('/')[-1], self.chain1_id,
                 self.chain2.split('/')[-1], self.chain2_id))
        tmpfile.close()
        os.system('%s/gramm scan'%self.TOOL_PATH)
        if os.path.exists('ch1-ch2.res'):
            os.system('cp ch1-ch2.res %s'%self.solution)
            self._fix_path_in_solution()
            print 'Save result to', self.solution
            return True
        else: ## show error
            print 'Error in docking:'
            os.system('cat gramm.log')
            return False
        os.chdir(self.OLD_PATH)

    def generate_complex(self, n=10):
        if n < 0: n = self.solution_number()
        os.chdir(self.TEMP_PATH)
        tmpfile = open('wlist.gr', 'w')
        tmpfile.write(
'''# File_of_predictions   First_match   Last_match   separate/joint  +init_lig
# ----------------------------------------------------------------------------
ch1-ch2.res 1 %s separ no
            '''%(n))
        tmpfile.close()
        os.system('cp %s ch1-ch2.res'%self.solution)
        os.system('%s/gramm coord'%self.TOOL_PATH)
        os.system('cat gramm.log')
        os.chdir(self.OLD_PATH)
        return ['%s/ch1-ch2_%s.pdb'%(self.TEMP_PATH, i+1) for i in range(n)]

    def solution_scores(self, skip=31):
        scores = []
        with open(self.solution, 'r') as tempfile:
            for i in xrange(skip):
                tempfile.readline()
            for l in tempfile:
                scores.append(float(l.split()[1]))
        return scores

## Ref: http://biopython.org/DIST/docs/api/Bio.PDB.Dice-pysrc.html
import re
_hydrogen=re.compile("[123 ]*H.*")
class SelectChain(Dice.ChainSelector):
    """ Only accepts residues within the right chainids.
        Rename the chainids to be [new_id, '1'+new_id, '2'+new_id, ...].
        Slice all chains to be between start and end.
        Rotate the structure by some degree on x axel
    """ 

    def __init__(self, chain_id, new_id=None, start=None, end=None, rotate_seed=None):
        self.chain_id = chain_id
        if new_id == None:
            self.new_id = chain_id
        else:
            self.new_id = new_id
        self.start = int(start) if start != None else start
        self.end = int(end) if end != None else end
        self.model_id = 0
        if rotate_seed == None:
            self.rotate_x = 0
            self.rotate_y = 0
            self.rotate_z = 0
        else:
            import random
            random.seed(rotate_seed)
            self.rotate_x = random.randint(1, 360)
            self.rotate_y = random.randint(1, 360)
            self.rotate_z = random.randint(1, 360)

    def accept_atom(self, atom):
        # do rigid body transformation
        atom.transform(rotaxis(self.rotate_x, Vector(1,0,0)), [0,0,0])
        atom.transform(rotaxis(self.rotate_y, Vector(0,1,0)), [0,0,0])
        atom.transform(rotaxis(self.rotate_z, Vector(0,0,1)), [0,0,0])
        return 1
        # remove hydrogen atom
        #name=atom.get_id()
        #if _hydrogen.match(name): 
        #    return 0 
        #else: 
        #    return 1 

    def accept_chain(self, chain):
        idx = self.chain_id.find(chain.get_id())
        if idx >= 0:
            chain.id = self.new_id[idx] ## one to one map
            return 1
        return 0

    def accept_residue(self, residue):
        hetatm_flag, resseq, icode=residue.get_id()
        if hetatm_flag != " ":
            return 0
        if icode != " ":
            return 0
        if self.start != None and self.start > resseq:
            return 0
        if self.end != None and resseq > self.end:
            return 0
        return 1

def save_chain(filename, structure, chainID, newID=None, st=None, ed=None, rotateSeed=None):
    sel=SelectChain(chainID, newID, st, ed, rotateSeed)
    io=PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)

def save_pdb_chain(pdbfile, chain, outfile, st=None, ed=None):
    warnings.filterwarnings("ignore")
    par = PDBParser()
    s = par.get_structure('PDBID', pdbfile)
    save_chain(outfile, s, chain, None, st, ed, None)

def combine_chains(inlist, outfile):
    with open(outfile, 'w') as out_obj:
        for infile in inlist:
            with open(infile, 'r') as in_obj:
                for line in in_obj:
                    if line.startswith('ATOM'):
                        out_obj.write(line)

def main(para):
    if 'PDB1' not in para:
        para['PDB1'] = '1JU5'
    if 'Ch1' not in para:
        para['Ch1'] = 'A'
    if 'PDB2' not in para:
        para['PDB2'] = '1JU5'
    if 'Ch2' not in para:
        para['Ch2'] = 'C'
    dock1 = UseZDOCK(temp_path='dock_temp')
    dock2 = UsePatchDock(temp_path='dock_temp')
    #dock3 = UseGRAMM(temp_path='dock_temp')
    for dock in [dock1, dock2]:
        show(dock)
        dock.prepare_data(para['PDB1'], para['Ch1'], para['PDB2'], para['Ch2'])
        show(dock.to_old_chain)
        dock.dock_them()
        slist = dock.generate_complex(3)
        assert os.path.exists(slist[0])

        dock.prepare_data(para['PDB2'], para['Ch2'], para['PDB1'], para['Ch1'])
        show(dock.to_old_chain)
        dock.dock_them()
        slist = dock.generate_complex(3)
        assert os.path.exists(slist[0])

        show(dock.solution)
        show(dock.solution_number())
        show(dock.solution_scores()[:3])
        show()

if __name__ == "__main__": main_fun(main)
