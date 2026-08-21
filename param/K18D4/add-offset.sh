#!/bin/bash

tmp=./temp.txt
echo $2
echo "offset $2" > ${tmp}
cat $1 >> ${tmp}
mv $1 $1.bk
mv ${tmp} $1
