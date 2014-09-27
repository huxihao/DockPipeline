cat data/500_HBHQ_Docking_Candidates.txt data/rest_HBHQ_Docking_Candidates.txt > work/test_list.txt
predict_interaction.py PoolPath=~/pub/dock_pool WorkPath=work NewList=test_list.txt ListFormat=p1/p2/pdb1/pdb2/a1/a2 FeatureType=SaveResidue ThreadNum=1 ListSize=5
predict_interaction.py PoolPath=~/pub/dock_pool WorkPath=work NewList=test_list.txt ListFormat=p1/p2/pdb1/pdb2/a1/a2 FeatureType=SaveResidueFile ThreadNum=5 ListSize=-1 OutMethod=RandomForest PredictCutoff=0.218


