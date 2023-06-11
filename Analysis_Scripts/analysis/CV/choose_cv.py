#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os, sys
import glob
import shutil
import random

def choose_2cv(finp,fout,cv1,cv2,cvv1,cvv2):
    out = open(fout, 'w')
    with open(finp) as f: 
        lines = f.readlines()
    f.close()
    out.write(lines[0])
    line1 =lines[0].split()
    for i in range(len(line1)): 
        if line1[i] == cv1:
            idx1 = i
        if line1[i] == cv2:
            idx2 = i
    for l in lines[1:]:
        sp1 = float(l.split()[idx1-2])
        sp2 = float(l.split()[idx2-2]) 
        if sp1 < cvv1 and sp2 > cvv2:
            out.write(l)
    out.close()

path = '/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/013_water/006_256H2O_slab'
finp  = os.path.join(path,'COLVAR_tmp')
fout  = os.path.join(path,'COLVAR_new')
choose_2cv(finp,fout,'cp50','cm50',1.25,-1.25)
