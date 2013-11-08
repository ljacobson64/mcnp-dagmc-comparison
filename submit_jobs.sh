#!/bin/bash

num_files=`ls -l | grep 'cmd' | wc -l | awk '{print $1}'`
echo $num_files

for (( i = 1 ; i <= $num_files ; i++)) ; do 
    job_names=`ls -l | grep 'cmd' | sed -n "$i"p | awk '{print $9}'`
    echo condor_submit $job_names
    condor_submit $job_names
done