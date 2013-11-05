#!/bin/bash
# Test script for the SourceField program

# bc based calculator which accepts scientific notation
function float_eval()
{
	# the horrible expression here convers numbers from scientific notation to something readable by bc
	echo $* | sed 's/--/+/g' | sed 's/\([0-9]*\(\.[0-9]*\)\?\)[eE]+\?\(-\?[0-9]*\)/(\1*10^\3)/g;s/^/scale=15;/'| bc -l
}

if [ -z $1 ]
then
	echo Error, missing target
	exit 1
else
	target=$1
fi

##########################################################################################################
echo Testing against the van der Pauw method
DV=$(./$target -q Test/4probe_test)

# Corner voltages should be +/- ln(2)/(2*pi) = +/- 0.1103178000763260
DV=$(float_eval "($DV-0.2206356001526520)/0.2206356001526520")
DV=${DV#-} # dirty absolute value hack

# accuracy in 4probe_test set to 1e-8 for the error current
# This should be a measure for the error in voltage, I assume 
# we should have the same order of magnitude and test whether the
# error is smnaller than 10 x the error set in 4probe_test
if [ $(float_eval "$DV<1e-7") -eq 1 ]
then
	printf "Success, %.3e %% deviation\n" $(float_eval "$DV*100")
else
	printf "Failed!, %.3e %% deviation\n" $(float_eval "$DV*100")
fi
echo


##########################################################################################################
echo Testing two contacts on infinite lamina
RR=3.8853120030
R=$(./$target -q Test/InfiniteLamina_test)
DR=$(float_eval "($R-$RR)/$RR")
DR=${DR#-} # dirty absolute value hack
if [ $(float_eval "$DR<1e-3") -eq 1 ]
then
	printf "Success, %.3e %% deviation\n" $(float_eval "$DR*100")
else
	printf "Failed!, %.3e %% deviation\n" $(float_eval "$DR*100")
fi
echo


##########################################################################################################
echo 'Testing 90 degrees rotation (deviation compared to previous calculation)'
R9=$(./$target -q Test/InfiniteLamina_test_90degrees)
DR=$(float_eval "($R9-$R)/$R")
DR=${DR#-} # dirty absolute value hack
if [ $(float_eval "$DR<1e-12") -eq 1 ]
then
	printf "Success, %.3e %% deviation\n" $(float_eval "$DR*100")
else
	printf "Failed!, %.3e %% deviation\n" $(float_eval "$DR*100")
fi

