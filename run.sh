mkdir work/zdock_pdbs
mkdir work/zdock_pool
cat data/500_HBHQ_Docking_Candidates.txt data/rest_HBHQ_Docking_Candidates.txt > work/full_list.txt
python src/predict_interaction.py WorkPath=work NewList=full_list.txt ListFormat=p1/p2/pdb1/pdb2/a1/a2 ThreadNum=1 ListSize=2

