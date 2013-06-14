# Parameters for the EPW run

# common parameters
ECUT=80.0 # energy cut in Ry
KMESH=[6,6,6] # automatic grid with no shift used for scf, nscf, and phonons

# stuff to be computed in the second EPW run
# and fine mesh to be used
# and smearings to be used
EPW_SECOND_RUN=\
"""
  phonselfen  = .true.
  a2f         = .false.
  elecselfen  = .false.

  nqf1 = 12
  nqf2 = 12
  nqf3 = 12
  nkf1 = 36
  nkf2 = 36
  nkf3 = 36

  fsthick     = 0.05
  eptemp      = 300.0
  degaussw    = 0.045
"""

# dft parameters for the QE runs: SCF, NSCF, PH
# only some parameters need to be specified here
# everything else will be added in the script itself
BLOCKS={
# scf and nscf data 
    "control":\
"""
pseudo_dir = '/global/homes/s/scoh/work/babio3/runs/psp',
""",
    "system":\
"""
ibrav = 0,
celldm(1) = 1.0,
nat=  5 ,
ntyp = 3 ,
ecutwfc = """+str(ECUT)+"""
occupations = 'smearing'
degauss = 0.005
""",
    "electrons":\
"""
mixing_beta = 0.7
conv_thr =  1.0d-11
""",
    "ATOMIC_SPECIES":\
"""
Xx  88.214 Ba_0.5__K_0.5.UPF
Bi  209.98  Bi.pbe-d-mt.UPF
O   15.999  08-O.GGA.fhi.UPF
""",
    "atomic_positions_crystal":\
"""
Xx   0.00  0.00  0.00
Bi   0.50  0.50  0.50
O    0.00  0.50  0.50
O    0.50  0.00  0.50
O    0.50  0.50  0.00
""",
    "cell_parameters_alat":\
"""
8.289090911   0.000000000   0.000000000
0.000000000   8.289090911   0.000000000
0.000000000   0.000000000   8.289090911
""",
# ph data 
    "inputph":\
"""
tr2_ph   =  1.0d-13
""",
# epw data 
    "inputepw":\
"""
nbndsub     =  9
nbndskip    =  13
num_iter    = 1000
dis_win_max = 20.0
dis_win_min = -20.0
dis_froz_max= 11.78
dis_froz_min= 2.5
proj(1)     = 'O:p'
""",
    }
