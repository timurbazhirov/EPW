#!/bin/bash
#PBS -l mppwidth=480
#PBS -l walltime=02:00:00
#PBS -S /bin/bash
#PBS -m ae
#PBS -N test_epw
#PBS -q premium
#PBS -M PUTEMAIL
#PBS -V

cd $PBS_O_WORKDIR

# make folder on scratch where everything in the calculation is stored
# link that folder from scratch locally in link named "run"
# if link "run" already exists, then skip this step
if [ ! -L run ]; then
 TMPFOLDER=$SCRATCH/_tmp_epw_$RANDOM_$RANDOM
 mkdir $TMPFOLDER
 ln -s $TMPFOLDER run
fi

# make folder for stdout and stderr output
if [ ! -L out ]; then
 mkdir out
fi

# load modules needed for executables
module unload PrgEnv-pgi
module load PrgEnv-intel
module load mkl/12.1.3.293

# executables used in the PH part of the calculation
export   PH_SCF_EXEC="aprun -n  24 /global/homes/s/scoh/compile/epw/01/comp/bin/pw.x     -npool 4"
export  PH_NSCF_EXEC="aprun -n  24 /global/homes/s/scoh/compile/epw/01/comp/bin/pw.x     -npool 4"
export    PH_PH_EXEC="aprun -n  24 /global/homes/s/scoh/compile/epw/01/comp/bin/ph.x     -npool 4"
# executables used in the EPW part of the calculation
export EPW_NSCF_EXEC="aprun -n  24 /global/homes/s/scoh/compile/epw/01/comp/bin/pw.x      -npool 4"
export  EPW_EPW_EXEC="aprun -n 216 -N 12 /global/homes/s/scoh/compile/epw/01/comp/EPW/bin/epw.x -npool 216"


#############################
### do everything at once ###
#############################
#
#python do_epw.py >& out/_python_all


#############################
###  parallelize phonons  ###
#############################
#
# first, do scf
python do_epw.py scf >& out/_python__scf
#
# second, do parallel phonons
PHTMPA=`wc _phonon_list`; PHTMPB=($PHTMPA);NUMPH=${PHTMPB[0]}
for i in $(seq 1 $NUMPH)
do
  python do_epw.py ph_one $i >& out/_python__ph_one_$i &
done
wait
#
# third, do an epw calculation
python do_epw.py >& out/_python__epw
