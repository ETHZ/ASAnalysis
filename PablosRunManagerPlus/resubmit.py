#!/usr/bin/python

#
# Resubmits specific jobs
#
# Takes a list of batch scripts as input

import sys,re
from os import path
import os,subprocess
from optparse import OptionParser

def getSeFile(file,id,opt):
   # open the script and get the remote directory
   
   sehome = ''
   sesubdir = ''
   hnname = ''
   sePat  = re.compile('^USER_SRM_HOME\s*=\s*"?(srm:[^"]*)')
   subPat = re.compile('^SEUSERSUBDIR\s*=\s*"?([^"]*)')
   hnPat  = re.compile('^HN_NAME\s*=\s*"?([^"]*)')
   f = open(file)
   try:
      for line in f:
         line = line.rstrip() # Remove trailing endline
         m = re.search(sePat,line)
         if m: sehome = m.group(1)
         m = re.search(subPat,line)
         if m: sesubdir = m.group(1)
         m = re.search(hnPat,line)
         if m: hnname = m.group(1)
         if sehome!='' and sesubdir != '' and hnname != '': break
   finally: f.close()
   sepath = path.join(sehome,hnname,sesubdir,'output_'+id+'.root')
   
   if opt.verb>2: 
      print '  Found remote home',sehome,', HN name',hnname,'and subdir',sesubdir
      print '  Returning',sepath
   return sepath


def process(file,opt):
   # Get job id
   id = 0
   filePat = re.compile('myScript_(\d+)\.sh')
   m = re.search(filePat,file)
   if m: id = m.group(1)
   else:
      print '*** FATAL: Couldn\'t get job id for',file
      sys.exit(-1)

   print "Getting remote file name",
   sefile = getSeFile(file,id,opt)             
   print sefile
   
   print "File to remove:",sefile
   if opt.dryrun: return

   # 1. Remove remote files
   # Check existence of file
   direct = sefile.replace('srm://t3se01.psi.ch:8443/srm/managerv2?SFN=','')
   cmd = ['ls',direct]
   err = subprocess.call(cmd)
   if err:
      print '*** File does not seem to exist:',err
   else:
      answ = 'y'
      cmd = ['lcg-del','-b','-D','srmv2','-l',sefile]
      if opt.verb>1: print "CMD =",' '.join(cmd)
      if not opt.force: 
         answ = raw_input('Remove file? [y/n] ')
         while answ != 'y' and answ != 'n': 
            answ = raw_input('Please answer y or n: ').rstrip()
      if answ == 'y':
         err = subprocess.call(cmd)
         if err:
            print '*** Error:',err
            sys.exit(-1)
         print "--> Removed remote file"
  
   # 2. Remove local file
   answ = 'y'
   log = file.replace('myScript','job').replace('sh','out')
   cmd = ['rm','-f','-v',log]
   if opt.verb>1: print "CMD =",' '.join(cmd)
   if not opt.force: 
      answ = raw_input('Remove local output file? [y/n] ')
      while answ != 'y' and answ != 'n': 
         answ = raw_input('Please answer y or n: ').rstrip()
   if answ == 'y':
      err = subprocess.call(cmd)
      if err:
         print '*** Error:',err
         sys.exit(-1)
      print "--> Removed local file"

   # 3. Resubmit
   cmd = ['qsub','-q','all.q','-N',"RMG"+str(id)+'resub','-o',log,'-j','y',file,str(id)]
   if opt.verb>1: print "CMD =",' '.join(cmd)
   if options.dryrun:
      return id
   err = subprocess.call(cmd)

if __name__ == '__main__' :
   usage = """%prog <myScript_1> [... <myScript_n>]
or     %prog <job_1.out> [... <job_n.out>]"""
   parser = OptionParser(usage=usage)
   parser.add_option("-v","--verbose",dest="verb",metavar="LEVEL",default=0,
                     help="Turn on verbosity")
   parser.add_option("-n","--dry-run",dest="dryrun",action="store_true",default=False,
                     help="Do not actually remove the files: just list them")
   parser.add_option("-f","--force",dest="force",action="store_true",default=False,
                     help="Force removal (don't ask questions)")
   (options,args) = parser.parse_args(args=None, values=None)

   if len(args)<1:
       parser.print_usage()
       sys.exit(-1)

   for file in args:
      if file.find('job_'):
         file = file.replace('job','myScript').replace('out','sh')
      print "+++ Resubmitting",file
      process(file,options)
