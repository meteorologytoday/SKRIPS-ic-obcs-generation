#!/home/rus043/anaconda2/bin/python
import sys, os
import numpy as np
runFile = 'gen_ic_wave_megh.m'
# dayList = np.arange(1,15,1,dtype=int)
# dayList = np.arange(11,32,1,dtype=int)
dayList = np.arange(1,12,1,dtype=int)
print dayList

for dayI in dayList:

  # cmd = "sed -i \"18s@.*@eval([\'load \' hycompath \'0" + str(273+dayI) + '_' + str(20151000+dayI) + "00.mat\'])@\" " + runFile
  cmd = "sed -i \"18s@.*@eval([\'load \' hycompath \'0" + str(304+dayI) + '_' + str(20151100+dayI) + "00.mat\'])@\" " + runFile
  os.system(cmd)
  cmd = 'matlab -nodisplay <' + runFile
  os.system(cmd)
