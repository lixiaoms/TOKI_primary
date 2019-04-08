#!/bin/bash
function usage(){
    echo "Usage: path/to/TOKI.sh <Hi-c matrix file> [options] 

   Options:
     -m <ctcf motif file> 
     -e <expression file> 
     -b <resolution of matrix> (default=40)
     -c <position of centromere> (default=0,0)
     -n <number of cores> (default=1)
     -o <output dir> (default=outs)"
}

if [ x$1 != x ]
then
    hic=$1
    shift
else
    usage
    exit
fi

b='40'
c='0,0'
n='1'
o='outs'
opt_num=0
path=$0
dir=${path%'/TOKI.sh'}

while getopts "h:m:e:r:c:n:o:" arg
do
        case $arg in
             m)
                m=$OPTARG
                opt_num=`expr $opt_num + 1`
                ;;
             e)
                e=$OPTARG
                opt_num=`expr $opt_num + 1`
                ;;
             b)
                b=$OPTARG
                opt_num=`expr $opt_num + 1`
                ;;
             c)
                c=$OPTARG
                opt_num=`expr $opt_num + 1`
                ;;
             n)
                n=$OPTARG
                opt_num=`expr $opt_num + 1`
                ;;
             o)
                o=$OPTARG
                opt_num=`expr $opt_num + 1`
                ;;
             ?)
                usage
                exit
        esac     
done


if [ $# != `expr $opt_num \* 2` ]
then
    usage
    exit
fi


if [ ! -d "$o" ]
then
  mkdir -p "$o"
fi


echo 'TOKI execution begins'
date +%r


state=`Rscript $dir/script/floyd.r $hic $b $c $n $o`
if [[ $state = '[1] 1' ]]
then
    echo 'distance feature matrix is created'
    date +%r
else
    echo 'distance feature matrix creation failed'
    date +%r
fi


if [ $m ]
then
    state=`Rscript $dir/script/fimotocsv.r $hic $m $b $c $o`
    if [[ $state = '[1] 1' ]]
    then
        echo 'motif feature matrix is created'
        date +%r
    else
        echo 'motif feature matrix creation failed'
        date +%r
    fi
fi

if [ $e ]
then
    state=`Rscript $dir/script/tabtocsv.r $hic $e $b $c $o`
    if [[ $state = '[1] 1' ]]
    then
        echo 'expression feature matrix is created'
        date +%r
    else
        echo 'expression feature matrix creation failed'
        date +%r
    fi
fi


state=`Rscript $dir/script/insulation.r $hic $b $c $o`
if [[ $state = '[1] 1' ]]
then
    echo 'TOKI execution successed'
        date +%r
else
    echo 'TOKI execution failed'
        date +%r
fi
