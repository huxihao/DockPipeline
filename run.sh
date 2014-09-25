cat data/set1_all_new.txt > work/train_list.txt
cat data/500_HBHQ_Docking_Candidates.txt data/rest_HBHQ_Docking_Candidates.txt > work/test_list.txt
python src/predict_interaction.py WorkPath=work ListFile=train_list.txt NewList=test_list.txt ListFormat=p1/p2/pdb1/pdb2/a1/a2 FeatureType=SaveResidueFile ThreadNum=5 ListSize=5 PredictCutoff=0.218

