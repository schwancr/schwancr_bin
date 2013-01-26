
DataDir="/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1/AnalyzeTrajs"
RMSD_Dir="/home/schwancr/CheY_smog/fahCheY/Calpha/SkipFrame1"

avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/Qnn_SkipFrame1.npy -o stateAvg_Qnn.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/Qcc_Raw_SkipFrame1.npy -o stateAvg_Qcc.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/QtotRaw_SkipFrame1.npy -o stateAvg_Qtot.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/b3b4_SkipFrame1.npy -o stateAvg_b3b4.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/b4b5_Raw_SkipFrame1.npy -o stateAvg_b4b5.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/a2a3_Raw_SkipFrame1.npy -o stateAvg_a2a3.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $DataDir/a3a4_Raw_SkipFrame1.npy -o stateAvg_a3a4.dat 
avg4stateMP.py -a Assignments.h5 -P 24 -d $RMSD_Dir/RawRMSDs.npy -o stateAvg_RMSD.dat 
