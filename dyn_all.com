#!/bin/sh
# Dynamics of RLF circuit.

fla=d4
flb=d1
flc=d2
fld=d5

./rlf.ex $fla
./rlf.ex $flb
./rlf.ex $flc
./rlf.ex $fld

echo "-1000.0 0.0" > start.xx
echo "-0.0 0.0" >> start.xx

cp start.xx rlf.col.$fla.xx.MR
awk '{print $1 / 1000.0, $2}' rlf.col.$fla >> rlf.col.$fla.xx.MR
cp start.xx rlf.col.$fla.xx.ML
awk '{print $1 / 1000.0, $3}' rlf.col.$fla >> rlf.col.$fla.xx.ML
cp start.xx rlf.col.$fla.xx.MF
awk '{print $1 / 1000.0, $4}' rlf.col.$fla >> rlf.col.$fla.xx.MF

cp start.xx rlf.col.$flb.xx.MR
awk '{print $1 / 1000.0, $2}' rlf.col.$flb >> rlf.col.$flb.xx.MR
cp start.xx rlf.col.$flb.xx.ML
awk '{print $1 / 1000.0, $3}' rlf.col.$flb >> rlf.col.$flb.xx.ML
cp start.xx rlf.col.$flb.xx.MF
awk '{print $1 / 1000.0, $4}' rlf.col.$flb >> rlf.col.$flb.xx.MF

cp start.xx rlf.col.$flc.xx.MR
awk '{print $1 / 1000.0, $2}' rlf.col.$flc >> rlf.col.$flc.xx.MR
cp start.xx rlf.col.$flc.xx.ML
awk '{print $1 / 1000.0, $3}' rlf.col.$flc >> rlf.col.$flc.xx.ML
cp start.xx rlf.col.$flc.xx.MF
awk '{print $1 / 1000.0, $4}' rlf.col.$flc >> rlf.col.$flc.xx.MF

cp start.xx rlf.col.$fld.xx.MR
awk '{print $1 / 1000.0, $2}' rlf.col.$fld >> rlf.col.$fld.xx.MR
cp start.xx rlf.col.$fld.xx.ML
awk '{print $1 / 1000.0, $3}' rlf.col.$fld >> rlf.col.$fld.xx.ML
cp start.xx rlf.col.$fld.xx.MF
awk '{print $1 / 1000.0, $4}' rlf.col.$fld >> rlf.col.$fld.xx.MF

xmgrace -graph  0 rlf.col.$fla.xx.MR \
        -graph  1 rlf.col.$fla.xx.MF \
        -graph  2 rlf.col.$fla.xx.ML \
        -graph  3 rlf.col.$flb.xx.MR \
        -graph  4 rlf.col.$flb.xx.MF \
        -graph  5 rlf.col.$flb.xx.ML \
        -graph  6 rlf.col.$flc.xx.MR \
        -graph  7 rlf.col.$flc.xx.MF \
        -graph  8 rlf.col.$flc.xx.ML \
        -graph  9 rlf.col.$fld.xx.MR \
        -graph 10 rlf.col.$fld.xx.MF \
        -graph 11 rlf.col.$fld.xx.ML \
        -hdevice EPS -p dyn_all.gr -printfile dyn_all.$fla.eps

/bin/rm start.xx
/bin/rm rlf.col.$fla.xx.MR rlf.col.$fla.xx.ML rlf.col.$fla.xx.MF
/bin/rm rlf.col.$flb.xx.MR rlf.col.$flb.xx.ML rlf.col.$flb.xx.MF
/bin/rm rlf.col.$flc.xx.MR rlf.col.$flc.xx.ML rlf.col.$flc.xx.MF
/bin/rm rlf.col.$fld.xx.MR rlf.col.$fld.xx.ML rlf.col.$fld.xx.MF
