#!/bin/sh

# Automated checks for pw.x - PG 2007
# Some specific quantities are checked against a reference output
# Checks are implemented for the following calculations: 
#   'scf', 'relax', 'md', 'vc-relax', 'neb', 'nscf', 'meta'
# (see below for the three latter)
#
# Input data: *.in, reference results: *.res, output: *.out
# ./check-pw.x.j checks all *.in files
# ./check-pw.x.j "some file(s)" checks the specified files
# Example:
#    ./check-pw.x.j atom*.in lsda*
# If you want to save a copy in file "logfile":
#    ./check-pw.x.j atom*.in lsda* | tee logfile
#
# For 'meta' calculations, the quantity that is verified is
#    the averaged configurational energies
#
# For 'neb' calculations, the quantities that are verified are
#    the converged activation energy
#    the number of neb iterations
#
# For 'nscf' case, the data is in file $name.in2, where $name.in is the
# data for the scf calculation that must be executed before the nscf one.
# Output is written to $name.out2 and checked vs reference data $name.res2
# The quantities that are compared with reference ones are:
#    the Fermi energy, or
#    the HOMO and LUMO
#    the total polarization (for the Berry's phase calculation)
#
# For all other cases, the quantites that are verified are:
#    the converged total energy
#    the number of scf iterations
#    the module of the force ( sqrt(\sum_i f_i^2)) if calculated;
#    the pressure P if calculated

# taken from examples - not sure it is really needed
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

ESPRESSO_ROOT=`cd .. ; pwd`
PARA_PREFIX=
#PARA_PREFIX="mpirun -np 2"
PARA_POSTFIX=
ESPRESSO_TMPDIR=$ESPRESSO_ROOT/tmp/
ESPRESSO_PSEUDO=$ESPRESSO_ROOT/pseudo/

# no need to specify outdir and pseudo_dir in all *.in files
export ESPRESSO_TMPDIR ESPRESSO_PSEUDO

if test ! -d $ESPRESSO_TMPDIR
then
   mkdir $ESPRESSO_TMPDIR
fi

# this is the current directory, where the test is executed
TESTDIR=`pwd`

# With no arguments, checks all *.in files
# With an argument, checks files (ending with .in) matching the argument

if test $# = 0
then
    files=`/bin/ls *.in`
else
    files=`/bin/ls $*| grep "\.in$"`
fi

########################################################################
# function to test matadynamics - usage: check_meta "file prefix"
########################################################################
check_meta () {
  # get average configurational energy (truncated to 4 significant digits)
  e0=`grep 'Final energy' $1.ref | awk '{sum+=$4} END {printf "%8.4f\n", sum/NR}'`
  e1=`grep 'Final energy' $1.out | awk '{sum+=$4} END {printf "%8.4f\n", sum/NR}'`
  #
  if test "$e1" = "$e0"
  then
    $ECHO  "passed"
  fi
  if test "$e1" != "$e0"
  then
    $ECHO "discrepancy in average configurational energy detected"
    $ECHO "Reference: $e0, You got: $e1"
  fi
}
########################################################################
# function to test NEB calculations - usage: check_neb "file prefix"
########################################################################
check_neb () {
  # get reference number of neb iterations
  n0=`grep 'neb: convergence' $1.ref | awk '{print $1}'`
  # get reference activation energy (truncated to 4 significant digits)
  e0=`grep 'activation energy' $1.ref | tail -1 | awk '{printf "%8.4f\n", $5}'`
  #
  n1=`grep 'neb: convergence' $1.out | awk '{print $1}'`
  e1=`grep 'activation energy' $1.out | tail -1 | awk '{printf "%8.4f\n", $5}'`
  if test "$e1" = "$e0"
  then
    if test "$n1" = "$n0"
    then
      $ECHO  "passed"
    fi
  fi
  if test "$e1" != "$e0"
  then
    $ECHO "discrepancy in activation energy detected"
    $ECHO "Reference: $e0, You got: $e1"
  fi
  if test "$n1" != "$n0"
  then
    $ECHO "discrepancy in number of neb iterations detected"
    $ECHO "Reference: $n0, You got: $n1"
  fi
}
########################################################################
# function to test scf calculations - usage: check_scf "file prefix"
########################################################################
check_scf () {
  # get reference total energy (cut to 6 significant digits)
  e0=`grep ! $1.ref | tail -1 | awk '{printf "%12.6f\n", $5}'`
  # get reference number of scf iterations
  n0=`grep 'convergence has' $1.ref | tail -1 | awk '{print $6}'`
  # get reference initial force (cut to 4 significant digits)
  f0=`grep "Total force" $1.ref | head -1 | awk '{printf "%8.4f\n", $4}'`
  # get reference pressure
  p0=`grep "P= " $1.ref | tail -1 | awk '{print $6}'`
  #
  # note that only the final energy, pressure, number of iterations, 
  # and only the initial force are tested - hopefully this should 
  # cover the various MD and optimization cases as well as simple scf
  #
  e1=`grep ! $1.out | tail -1 | awk '{printf "%12.6f\n", $5}'`
  n1=`grep 'convergence has' $1.out | tail -1 | awk '{print $6}'`
  f1=`grep "Total force" $1.out | head -1 | awk '{printf "%8.4f\n", $4}'`
  p1=`grep "P= " $1.out | tail -1 | awk '{print $6}'`
  #
  if test "$e1" = "$e0"
  then
    if test "$n1" = "$n0"
    then
      if test "$f1" = "$f0"
      then
        if test "$p1" = "$p0"
        then
          $ECHO  "passed"
        fi
      fi
    fi
  fi
  if test "$e1" != "$e0"
  then
    $ECHO "discrepancy in total energy detected"
    $ECHO "Reference: $e0, You got: $e1"
  fi
  if test "$n1" != "$n0"
  then
    $ECHO "discrepancy in number of scf iterations detected"
    $ECHO "Reference: $n0, You got: $n1"
  fi
    if test "$f1" != "$f0"
    then
    $ECHO "discrepancy in force detected"
    $ECHO "Reference: $f0, You got: $f1"
  fi
  if test "$p1" != "$p0"
  then
    $ECHO "discrepancy in pressure detected"
    $ECHO "Reference: $p0, You got: $p1"
  fi
}
########################################################################
# function to test nscf calculations - usage: check_nscf "file prefix" "number"
########################################################################
check_nscf () {
  # get reference Fermi energy
  ef0=`grep Fermi $1.ref$2 | awk '{print $5}'`
  # get reference HOMO and LUMO
  eh0=`grep "highest occupied" $1.ref$2 | awk '{print $7}'`
  el0=`grep "highest occupied" $1.ref$2 | awk '{print $8}'`
  # get total polarization (for Berry's phase calculation)
  tf0=`grep " P = " $1.ref$2 | head -1 | awk '{printf "%7.5f", $3}'`
  #
  ef1=`grep Fermi $name.out$n | awk '{print $5}'`
  eh1=`grep "highest occupied" $1.out$2 | awk '{print $7}'`
  el1=`grep "highest occupied" $1.out$2 | awk '{print $8}'`
  tf1=`grep " P = " $1.out$2 | head -1 | awk '{printf "%7.5f", $3}'`
  #
  if test "$ef1" = "$ef0"
  then
    if test "$eh1" = "$eh0"
    then
      if test "$el1" = "$el0"
      then
        if test "$tf1" = "$tf0"
        then
          $ECHO  "passed"
        fi
      fi
    fi
  fi
  if test "$ef1" != "$ef0"
  then
    $ECHO "discrepancy in Fermi energy detected"
    $ECHO "Reference: $ef0, You got: $ef1"
  fi
  if test "$eh1" != "$eh0"
  then
    $ECHO "discrepancy in HOMO detected"
    $ECHO "Reference: $eh0, You got: $eh1"
  fi
  if test "$el1" != "$el0"
  then
    $ECHO "discrepancy in LUMO detected"
    $ECHO "Reference: $el0, You got: $el1"
  fi
  if test "$tf1" != "$tf0"
  then
    $ECHO "discrepancy in polarization detected"
    $ECHO "Reference: $tf0, You got: $tf1"
  fi
}
########################################################################
# function to get wall times - usage: get_times "file prefix"
########################################################################
get_times () {
  # convert from "1h23m45.6s" to seconds
  # the following line prevents cases such as "2m 7.5s"
  grep 'wall time' $1.ref | sed 's/m /m0/' > $1.tmp 
  # in order to get cpu instead of wall time, replace $3 to $6
  tref=`awk '/wall time/ \
                { str = $6; h = m = s = 0;
                  if (split(str, x, "h") == 2) { h = x[1]; str = x[2]; }
                  if (split(str, x, "m") == 2) { m = x[1]; str = x[2]; }
                  if (split(str, x, "s") == 2) { s = x[1]; str = x[2]; }
                  t += h * 3600 + m * 60 + s; }
                END { printf("%.2f\n", t); }' \
               $1.tmp`
  # as above for file *.out
  grep 'wall time' $1.out | sed 's/m /m0/' > $1.tmp 
  tout=`awk '/wall time/ \
                { str = $6; h = m = s = 0;
                  if (split(str, x, "h") == 2) { h = x[1]; str = x[2]; }
                  if (split(str, x, "m") == 2) { m = x[1]; str = x[2]; }
                  if (split(str, x, "s") == 2) { s = x[1]; str = x[2]; }
                  t += h * 3600 + m * 60 + s; }
                END { printf("%.2f\n", t); }' \
               $1.tmp`
  /bin/rm $1.tmp
  # accumulate data
  totref=`echo $totref $tref | awk '{print $1+$2}'`
  totout=`echo $totout $tout | awk '{print $1+$2}'`
}

for file in $files
do
  name=`basename $file .in`
  $ECHO "Checking $name...\c"
  ###
  # run the code in the scratch directory
  #
  cd $ESPRESSO_TMPDIR
  $PARA_PREFIX $ESPRESSO_ROOT/bin/pw.x $PARA_POSTFIX < $TESTDIR/$name.in \
                                                     > $TESTDIR/$name.out
  if test $? != 0; then
     $ECHO "FAILED with error condition!"
     $ECHO "Input: $name.in, Output: $name.out, Reference: $name.ref"
     $ECHO "Aborting"
     exit 1
  fi
  #
  cd $TESTDIR
  ###
  if test -f $name.ref ; then
     # reference file exists
     if grep 'neb: convergence achieved' $name.ref > /dev/null; then
        #
        # Specific test for NEB
        #
	check_neb $name
        #
     elif grep 'calculation of the mean force' $name.ref > /dev/null; then
        #
        # Specific test for metadynamics
        #
	check_meta $name
        #
     else
        #
        # Test for scf/relax/md/vc-relax
        #
	check_scf $name
        #
     fi
     #
     # extract wall time statistics
     #
     get_times $name
     #
  else
     $ECHO  "not checked, reference file not available "
  fi
  #
  # now check subsequent non-scf step if required
  # look for $name.in2
  n=2
  if test -f $name.in$n; then
     $ECHO "Checking $name, step $n ...\c"
     ###
     # run the code in the scratch directory
     #
     cd $ESPRESSO_TMPDIR
     $PARA_PREFIX $ESPRESSO_ROOT/bin/pw.x $PARA_POSTFIX < $TESTDIR/$name.in$n \
                                                        > $TESTDIR/$name.out$n
     if test $? != 0; then
        $ECHO "FAILED with error condition!"
        $ECHO "Input: $name.in$n, Output: $name.out$n, Reference: $name.ref$n"
        $ECHO "Aborting"
        exit 1
     fi
     #
     cd $TESTDIR
     ###
     if test -f $name.ref$n ; then
        # reference file exists
        check_nscf $name $n
        # extract wall time statistics
        get_times $name
     else
        $ECHO  "not checked, reference file not available "
     fi
  fi

done

$ECHO  "Total wall time (s) spent in this run: " $totout
$ECHO  "Reference                            : " $totref

