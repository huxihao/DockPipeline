DOCK_PATH=`pwd`
chmod 777 src/*
export PATH=$PATH:$DOCK_PATH/lib/bin:$DOCK_PATH/src

## Set Data Path
export DOCK_PDB_PATH=/home/local/CORNELL/resources/pdb/data
export DOCK_MOD_PATH=/home/local/CORNELL/resources/instruct_data/models/dockml/single_structures

## Set EXE Path
export DOCK_NACCESS_PATH=$DOCK_PATH/bin/naccess2.1.1/naccess
export DOCK_ZDOCK_PATH=$DOCK_PATH/bin/zdock3.0.2_linux_x64
export DOCK_PATCHDOCK_PATH=$DOCK_PATH/bin/PatchDock
export DOCK_GRAMM_PATH=$DOCK_PATH/bin/gramm
export DOCK_SIFTS_PATH=$DOCK_PATH/data/SIFTS


