#!/bin/bash

rm *.o model Makefile

qmake simulation.pro
make 

iterations=10080
lactate=0
GOS=0
fl=0
lactose=211
ox=0
initox=0.1
drift=20
columns=225
rows=8
bmreplace=1
energysave=0
celldrift=0
bacplace=2
colonization=3.00005
death=0.15
growth=4
bacmove=1
initsize=0.01
expmode=0
flux=0.00025
dif=8
cellfrac=0.025

./model 1 ${iterations} ${lactate} ${GOS} ${fl} ${lactose} ${ox} ${initox} ${drift} ${columns} ${rows} ${bmreplace} ${energysave} ${celldrift} ${bacplace} ${colonization} ${death} ${growth} ${bacmove} ${initsize} ${expmode} ${flux} ${dif} ${cellfrac} samplerun &

