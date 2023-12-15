#! /usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo " usage: $0 [runlist] "
    exit 1
fi

SDIR="/home/preghenella/EIC/sipm4eic-testbeam2023-analysis/recoana/root"
DDIR="/home/preghenella/EIC/sipm4eic-testbeam2023-analysis/data"

LIST=$1
for RUN in $(cat ${LIST}); do
    while [ $(ps -ef | grep hough.C | grep -v grep | wc -l) -ge 8 ]; do sleep 1; done
    DINP="${DDIR}/${RUN}/process-data-v1.0/recodata.root"
    DOUT="${DDIR}/${RUN}/process-data-v1.0/hough.root"
    echo " --- start processing run : ${RUN} "
    root -b -q -l "${SDIR}/hough.C(\"${DINP}\", \"${DOUT}\")" &> /dev/null &
    sleep 1
done

wait
