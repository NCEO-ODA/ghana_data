#!/bin/bash
#BSUB -J "submit_lai[1-17]"
#BSUB -o /home/users/jlgomezdans//out.%I
#BSUB -e /home/users/jlgomezdans//err.%I
#BSUB -W 04:45
#BSUB –q short-serial
#BSUB -M 64000 

cd /gws/nopw/j04/odanceo/public/MCD15/
/home/users/jlgomezdans/miniconda3/bin/python ./to_tif.py   $((LSB_JOBINDEX + 2001)) 

