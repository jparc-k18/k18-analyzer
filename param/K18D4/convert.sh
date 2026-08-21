!/bin/bash

prefix=d4function_run

for i in *.txt;
do
	a=${i/${prefix}/};
#	echo $a;
	a=${a/.txt};
#	echo $a;
	b=${a/-*/};
	c=${a/*-/};
#	echo $i $b $c
	b=`printf %05d $b`;
	c=`printf %05d $c`;
#	echo $i $b $c
	echo "$i  -->  ${prefix}${b}-${c}.txt"	
	mv $i  ${prefix}${b}-${c}.txt
done
