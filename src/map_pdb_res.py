from tools import *
import os, gzip, urllib

DEFINE_SIFTS_PATH = os.path.dirname(os.path.realpath(__file__)) + '/../data/SIFTS'

def download_map(pdb):
    print 'Downloading the residue map for', pdb
    header = '{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}'
    local_file = pdb.lower()+'.xml.gz'
    try:
        urlstr = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'+local_file
        os.system('wget ' + urlstr)
    except:
        print 'Fail to download', urlstr
        return
    try:
        siftsfile = gzip.open(local_file,'rb')
        siftsxml = siftsfile.read()
        siftsfile.close()
        os.remove(local_file)
    except:
        print 'Fail to read', local_file
        return
    ## show contents
    from lxml import etree
    siftsdata = etree.fromstring(siftsxml)
    for i in siftsdata.findall(".//"+header+"crossRefDb[@dbSource='UniProt']"):
        if i.getparent().find(header+"residueDetail[@property='Annotation']") is not None:
            continue
        chainid = i.getparent().find(header+"crossRefDb[@dbSource='PDB']").get('dbChainId')
        pdbresnum = i.getparent().find(header+"crossRefDb[@dbSource='PDB']").get('dbResNum')
        pdbresname = i.getparent().find(header+"crossRefDb[@dbSource='PDB']").get('dbResName')
        chainuniprot = i.get('dbAccessionId')
        unpresnum = i.get('dbResNum')
        unpresname = i.get('dbResName')
        yield [pdb,chainid,pdbresnum,pdbresname,chainuniprot,unpresnum,unpresname]

PDB_RES_MAP = {}
def parse_maps(pdblist, mapfile=DEFINE_SIFTS_PATH+'/../pdbresiduemapping.txt.gz'):
    global PDB_RES_MAP
    if len(PDB_RES_MAP) == 0:
        if not os.path.exists(mapfile):
            print 'Can not find the processed file:', mapfile
            return
        if mapfile.endswith('.gz'):
            import gzip
            tfile = gzip.open(mapfile, 'rb')
            for line in tfile:
                ele = line.strip().split()
                maps = PDB_RES_MAP.get(ele[0].upper(), [])
                maps.append(ele)
                PDB_RES_MAP[ele[0].upper()] = maps
            tfile.close()
        else: ## assuming .txt
            with open(mapfile, 'r') as tfile:
                for line in tfile:
                    ele = line.strip().split()
                    maps = PDB_RES_MAP.get(ele[0].upper(), [])
                    maps.append(ele)
                    PDB_RES_MAP[ele[0].upper()] = maps
    data = []
    for pdb in pdblist:
        if pdb.upper() in PDB_RES_MAP:
            data += PDB_RES_MAP[pdb.upper()]
    for pdb, ch, p, map1, map2, map3 in data:
        r1 = map2[1:-1].split(',') ## PDB
        r2 = map1[1:-1].split(',') ## Uniprot
        l1 = []; l2 = []
        for r in r1:
            if r.find('-') > 0:
                st, ed = r.split('-')
                l1.extend(range(int(st), int(ed)+1))
            else:
                l1.append(r)
        for r in r2:
            if r.find('-') > 0:
                st, ed = r.split('-')
                l2.extend(range(int(st), int(ed)+1))
            else:
                l2.append(r)
        for j1, j2 in zip(l1, l2):
            yield [pdb, ch, str(j1), '', p, str(j2), '']

def process_map(pdb, data_path=DEFINE_SIFTS_PATH, download=True):
    if len(pdb) == 4:
        pdb = pdb.upper()
    num = 0
    if os.path.exists(data_path+'/'+pdb+'_map.txt'): ## processed map
        with open(data_path+'/'+pdb+'_map.txt', 'r') as tfile:
            for l in tfile:
                num += 1
                yield l.strip().split('\t')
    if num == 0 and download and len(pdb) == 4: ## get a new map from web server
        data = [d for d in download_map(pdb)] 
        with open(data_path+'/'+pdb+'_map.txt', 'w') as tfile:
            for d in data:
                tfile.write('\t'.join(d)+'\n')
                yield d

def pdblist_to_uniprot(pdb_list, download=True):
    ''' Search parsered pdb maps for batch mode '''
    pdblist = set()
    for pdb in pdb_list:
        if len(pdb) == 4:
            pdblist.add(pdb.upper())
        else:
            pdblist.add(pdb)
    m = {}
    known = set()
    for pdb in pdblist:
        if pdb not in known:
            add_chain = False
            for pdb,chid,pdbnum,tA,uni,unum,tB in process_map(pdb, download=False):
                m['%s:%s:%s'%(pdb, chid, pdbnum)] = '%s:%s'%(uni, unum)
                add_chain = True
            if add_chain:
                m['%s:%s'%(pdb, chid)] = uni ## chain map
    for pdb,chid,pdbnum,tA,uni,unum,tB in parse_maps([i for i in pdblist if i not in known]):
        m['%s:%s:%s'%(pdb, chid, pdbnum)] = '%s:%s'%(uni, unum)
        m['%s:%s'%(pdb, chid)] = uni ## chain map
        known.add(pdb)
    for pdb in pdblist:
        if pdb not in known:
            add_chain = False
            for pdb,chid,pdbnum,tA,uni,unum,tB in process_map(pdb, download=download):
                m['%s:%s:%s'%(pdb, chid, pdbnum)] = '%s:%s'%(uni, unum)
                add_chain = True
            if add_chain:
                m['%s:%s'%(pdb, chid)] = uni ## chain map
    return m

def pdb_to_uniprot(pdb, download=True):
    ''' Search individual pdb map first '''
    if len(pdb) == 4:
        pdb = pdb.upper()
    m = {}
    for pdb,chid,pdbnum,tA,uni,unum,tB in process_map(pdb, download=False):
        m['%s:%s:%s'%(pdb, chid, pdbnum)] = '%s:%s'%(uni, unum)
    if len(m) == 0:
        for pdb,chid,pdbnum,tA,uni,unum,tB in parse_maps([pdb]):
            m['%s:%s:%s'%(pdb, chid, pdbnum)] = '%s:%s'%(uni, unum)
    if len(m) == 0:
        for pdb,chid,pdbnum,tA,uni,unum,tB in process_map(pdb, download=download):
            m['%s:%s:%s'%(pdb, chid, pdbnum)] = '%s:%s'%(uni, unum)
    if len(m) != 0:
        m['%s:%s'%(pdb, chid)] = uni ## chain map
    return m

def main(para):
    if 'PDB' not in para:
        para['PDB'] = '2YL2'
    sifts_path = para['DataPath'] + '/SIFTS'
    if not os.path.exists(sifts_path):
        os.mkdir(sifts_path)
    # example list of pdb ids
    listofpdbs = ['1gc1','3k71','1JU5','1AK7','2DVS','1L2I','1EGL']
    for pdb in listofpdbs:
        maps = pdb_to_uniprot(pdb)
        print pdb, len(maps)
        if len(maps) > 0:
            print maps.keys()[0], maps[maps.keys()[0]]
    maps = pdblist_to_uniprot(listofpdbs)
    print pdb_to_uniprot('2xai') ## obsoleted PDB ID
    print 'Test Over!'
    for a in process_map(para['PDB'], data_path=sifts_path, download=True):
        print a

if __name__ == '__main__': main_fun(main)

