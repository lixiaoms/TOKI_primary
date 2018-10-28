#!/bin/bash
function usage(){
    echo "Usage: ./TOKI.sh <Hi-c matrix file> [options] 

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

mkdir $o
echo 'TOKI execution begins'
date +%r

Rscript script/floyd.r $hic $b $c $n $o
echo 'distance feature matrix created'
date +%r

if [ $m ]
then
    Rscript script/fimotocsv.r $hic $m $b $c $o
    echo 'motif feature matrix created'
    date +%r
fi

if [ $e ]
then
    Rscript script/tabtocsv.r $hic $e $b $c $o
    echo 'expression feature matrix created'
    date +%r
fi


Rscript script/insulation.r $hic $b $c $o
echo 'TOKI execution successed'
date +%r
