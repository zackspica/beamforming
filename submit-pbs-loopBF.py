import commands
from popen2 import popen2
import time, glob, os, sys
from pyrocko import model
import subprocess


def run(days):
  k=0
  for id, d in enumerate(days):
      k+=1
#      if k < 8290: continue 
      print d, k,'p'
      #continue
      # Open a pipe to the qsub command.
      output, input = popen2('qsub')
      # Customize your options here
      job_name = str(d)+'.BF'
      nnodes = 1
      processors = "nodes=1:ppn=%s"%nnodes
      command = "/home/zspica/local/bin/python2.7 main_beamforming.py %s" %d
      node = 'default'
      if k % 100.==0.:node= 'Beroza'
    #  if ista % 100.==0.:node= 'Q2a'
    #  if ibin % 100.==0.:node= 'Q2a'

      job_string = """
      #!/bin/bash\n\
      #PBS -N %s\n\
      #PBS -q %s\n\
      #PBS -l %s\n\
      #PBS -o qout/%s.out\n\
      #PBS -e qout/%s.err\n\
      cd $PBS_O_WORKDIR\n\
      %s""" % (job_name, node, processors, job_name, job_name, command)
      
      # Send job_string to qsub
      input.write(job_string)
      input.close()
      
      # Print your job and the system response to the screen as it's submitted
      print job_string
      print output.read()
      time.sleep(0.1000)
      
      os.system('echo "%s %s %s" >> o.txt'%(bin, str(d), k ))
      njob = int(commands.getoutput('qstat -u zspica | wc -l'))
      while njob >= 306:
          njob = int(commands.getoutput('qstat -u zspica | wc -l'))
          time.sleep(30)


if __name__=='__main__':
    try: os.makedir('qout/')
    except: pass
    days = []
    for n in xrange(4,26,1):
        days.append(n)
    print days
    run(days)
    
