#!/bin/bash

od=/mnt/isilon/xing_lab/shens/ISO_ReRun/Simulation/jum_output/
mkdir $od/1st_SJ
mv $od/*SJ.out.tab $od/1st_SJ
mv $od/*Log.final.out $od/1st_SJ
mv $od/*Log.progress.out $od/1st_SJ
mv $od/*Log.out $od/1st_SJ
mv $od/*Aligned.out.sam $od/1st_SJ
