mkdir work/zdock_pdbs
mkdir work/zdock_pool
cat data/set1_all_new.txt > work/set_list.txt
cat data/500_HBHQ_Docking_Candidates.txt data/rest_HBHQ_Docking_Candidates.txt > work/full_list.txt
python src/predict_interaction.py WorkPath=work ListFile=set_list.txt NewList=full_list.txt ListFormat=p1/p2/pdb1/pdb2/a1/a2 FeatureType=SaveResidue ThreadNum=10 ListSize=5
python src/predict_interaction.py WorkPath=work ListFile=set_list.txt NewList=full_list.txt ListFormat=p1/p2/pdb1/pdb2/a1/a2 FeatureType=SaveFinal ThreadNum=10 ListSize=-1 PredictCutoff=0.218

