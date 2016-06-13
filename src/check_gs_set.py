from tools import *

def get_dock_list(filename):
    out = []
    infile = open(filename, 'r')
    for line in infile:
        ele = line.split('\t')
        out.append(ele)
    infile.close()
    return out

def main(para):
    known = []
    known += get_dock_list('set1_all.txt')
    known += get_dock_list('set2_all.txt')
    known += get_dock_list('set3_top3.txt')
    gs_set = []
    with open('../data/GSset_docking.txt', 'r') as infile:
        for line in infile:
            p1, p2 = line.strip().split()
            show(p1)
            show(p2)
            for P1, P2, PDB1, CH1, PDB2, CH2, Mark in known:
                if p1 == P1 and p2 == P2:
                    gs_set.append((P1,P2,PDB1,CH1,PDB2,CH2))
                    show([PDB1, CH1, PDB2, CH2])
                    break
                elif p1 == P2 and p2 == P1:
                    gs_set.append((P2,P1,PDB2,CH2,PDB1,CH1))
                    show([PDB2, CH2, PDB1, CH1])
                    break
            show()
    show('We have\t%s\tknown docking pairs.\n'%len(gs_set))
    if not os.path.exists('GSset'):
        os.mkdir('GSset')
    from use_dock import UseZDOCK
    for p1, p2, pdb1, ch1, pdb2, ch2 in gs_set:
        dock1 = UseZDOCK(temp_path='GSset/%s_%s'%(p1,p2),
                         pool_path='../data/dock_pool')
        dock1.prepare_data(pdb1, ch1, pdb2, ch2)
        dock1.dock_them(10) 
        dock1.generate_complex(10)

if __name__ == '__main__': main_fun(main)
