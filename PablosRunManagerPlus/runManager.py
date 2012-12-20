#!/usr/bin/python

### WATCH OUT: NEEDS TO BE LAUNCHED AS python runManager.py configuration.cfg and NOT ./runManager configuration.cfg

###############################################################
#                       runManager.py                         #
###############################################################
# + Python script which allows to manage multiple jobs in the #
# context of a lsfbatch system.                               #
#                                                             #
# + type runManager -h to see a description of the usage      #   
#                                                             #
# + Feedback to: pablom@cern.ch                               #
# + Feedback about this version: contact marco-andrea         #
#                                                             #
############################################################### 
   

domultiprocessing=True

# System libraries to interact with the shell
import sys
import os
from os import popen
import commands
import time
from sys import stdout
try:
        from multiprocessing import Pool
except ImportError:
	print "Importing multiprocessing failed. Please use the script like this: python runManager.py template.cfg and not ./runManager.py template.cfg"
	domultiprocessing=False

# Command-line options -----------------------------------
from optparse import OptionParser
usage = """%prog <configuration>
where configuration should look like this:
---------------------------------------
ReleasePath path
Executable name
NameOfConf path
NameOfSource path
SrmPath path
WorkPath path
HPName path
---------------------------------------"""
parser = OptionParser(usage=usage)
parser.add_option("-n","--dry-run",dest="dryrun",action="store_true",default=False,
                  help="Do not actually submit the jobs: just creates file to submit")
parser.add_option("-v","--verbose",dest="verbose",action="store_true",default=False,
                  help="Turn verbosity on")
parser.add_option("-r","--renew",dest="renew",action="store_true",default=False,
                  help="Renew cached files")
(options, args) = parser.parse_args(args=None, values=None)
#---------------------------------------------------------

jobnumbers=[]
totaljobnumber=0

#-------------------------------------------------------------#
# ensure_dir method: Make sure a path exists                  #
#-------------------------------------------------------------#

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

#-------------------------------------------------------------#
# checklist method: Check whether the jobs are still running  #
# to permit post-processing                                   #
#-------------------------------------------------------------#

def checklist(joblist,currlist):
  jobfinished=[]
  for i in joblist:
    if i in currlist: ## if the job is still active
      currlist.remove(int(i))## we don't need this anymore on the list of currently running jobs (we already know it's running)
    else:#job is finished
      jobfinished.append(int(i)) 
  for job in jobfinished:
    joblist.remove(int(job))
  stdout.write("\r Jobs left: "+str(len(joblist))+" / " +str(totaljobnumber))
  stdout.flush()

#-------------------------------------------------------------#
# getFiles method: Obtains a list with the name of the files  #
# that will be used as input for a given job "step"           #        
# returns a string with the list                              # 
#-------------------------------------------------------------#
def getFiles(step, FilesPerJob, NumberOfJobs, OverloadedJobs, listOfFiles):
  
  if(step <= NumberOfJobs-OverloadedJobs):
    first = step*FilesPerJob
  else:
    first = FilesPerJob*(NumberOfJobs-OverloadedJobs)+(FilesPerJob+1)*(step-(NumberOfJobs-OverloadedJobs))-1

  counter = 0
  files = []
  for line in listOfFiles:
    if(step <= NumberOfJobs-OverloadedJobs):
      if(counter >= first and first+FilesPerJob > counter):
        files.append(line)
    else:
      if(counter >= first and first+FilesPerJob+1 > counter):
        files.append(line) 
    counter = counter + 1    
  
  return files
  #Check the ending in case it has a ',' 
  #if(files[-1][-2] == ','):
  #  files[-1] = files[-1][0:-2]
  #stringCad = ""
  #for fileLine in files:
  #  stringCad = stringCad + fileLine
  #return stringCad 
      

#-------------------------------------------------------------#
# createDirectory method: Creates a directory                 #
# returns 0 if everything was all right                       #
#-------------------------------------------------------------# 
def createDirectory(nameOfDirectory):
  result = os.system("mkdir " + nameOfDirectory)
  result2 = os.system("cd " + nameOfDirectory)
  return result 

#-------------------------------------------------------------#
# fixPath method: Just makes sure the path is ending          #
# with slash. Otherwise appends it                            #
#-------------------------------------------------------------# 
def fixPath(path):
  if(path[-1] == '/'):
    return path
  else:
    return path + "/"


#-------------------------------------------------------------#
# createCMSConf method: creates CMS configuration file        #
#                                                             #
#-------------------------------------------------------------# 
def createCMSConf(step, nameOfDirectory, releasePath, nameOfConf, inputString, executablePath, nameOfSRM, hpname, nameOfFolder, arguments):
  CFGFile = open(nameOfConf)
  cfgText = CFGFile.read()

  cadFile = "\""
  for file in inputString:
    cadFile = cadFile + file + " " 
  cadFile = cadFile + "\"" 

  if(nameOfFolder[len(nameOfFolder)-1] == "/"):
    nameOfFolder = nameOfFolder[0:len(nameOfFolder)-1]
  value = nameOfFolder.rfind("/")
  taskName = nameOfFolder[value+1:len(nameOfFolder)] 


  newcfgText1 = cfgText.replace("sourcefile" , cadFile)
  newcfgText2 = newcfgText1.replace("jobdir" , taskName + "/" + nameOfSRM + "_" + str(step))
  newcfgText3 = newcfgText2.replace("exefile", executablePath)
  newcfgText4 = newcfgText3.replace("srmdir", nameOfSRM + "/" + taskName)
  newcfgText5 = newcfgText4.replace("hpname", hpname)
  newcfgText6 = newcfgText5.replace("argumentsname", "\"" + arguments + "\"")
  newcfgText7 = newcfgText6.replace("therelease", releasePath)

  nameOfConf2 = nameOfConf.replace(".", "_"+str(step)+ ".")
 
  ensure_dir(str(nameOfDirectory + taskName + "/"))
  outputCFGFile = open(nameOfDirectory + taskName + "/" + nameOfConf2, "w")
  outputCFGFile.write(newcfgText7)
  CFGFile.close()
  outputCFGFile.close()

  outputName = "output_" + str(step) + ".root"
  stderr = nameOfDirectory + taskName + '/job_' + str(step) + '.err'
  stdout = nameOfDirectory + taskName + '/job_' + str(step) + '.out'
  
  thisjobnumber=0

  cmd = " ".join(['qsub','-q all.q','-N',"RMG"+str(step)+taskName,'-o',stdout,'-e',stderr,nameOfDirectory+taskName+'/'+nameOfConf2+' '+str(step)])
#  print cmd
  if options.dryrun: return thisjobnumber

  pipe=popen(cmd)#"qsub -e /tmp/ -o /tmp/ -N " + "RMJ"+str(step)+taskName + " " + nameOfDirectory + taskName + "/" + nameOfConf2 + " " + str(step))
  for l in pipe.readlines():
	  if l.find("Your job")>-1:
		  thisjobnumber=int(l[l.index('job ')+4:l.index(' (')])
		  print str(taskName)+": Submitted job "+str(step)+" with job number: "+str(thisjobnumber)
  return thisjobnumber
###############################################################
# createJob method: Prepares everything for a single job      #
#                                                             #
############################################################### 
def createJob(step, FilesPerJob, NumberOfJobs, OverloadedJobs, stringInfo, listOfFiles, nameOfFolder, arguments):
 
  #Definition of the different substrings needed
  nameOfSourceFile = stringInfo[0]
  releasePath = fixPath(stringInfo[1])
  executablePath = stringInfo[2]
  nameOfConfigurationFile = stringInfo[3]
  nameOfSRMPath = stringInfo[4]
  nameOfCurrentDirectory = stringInfo[5]
  hpname = stringInfo[6]
  
  inputFilesForCMS = getFiles(step, FilesPerJob, NumberOfJobs, OverloadedJobs, listOfFiles)
  if(inputFilesForCMS == ""):
    return "Error: No input files available"
  result = createCMSConf(step, nameOfCurrentDirectory, releasePath, 
           nameOfConfigurationFile, inputFilesForCMS, executablePath, nameOfSRMPath, hpname, nameOfFolder, arguments)
  return result

###############################################################
# parseInputFile method: Parses the input file                #
# More options could be added if needed                       #
############################################################### 

def parseInputFile(name):   
  showMessage("Parsing " + name + " file ...")
  file = open(name)
  stringInfo = [[],[],[],[],[], [], []]
  for line in file.readlines():
    splitLine = line.split()
    if(splitLine[0] == "ReleasePath"):
      stringInfo[1] = splitLine[1]
    if(splitLine[0] == "Executable"):
      stringInfo[2] = splitLine[1]
    if(splitLine[0] == "NameOfConf"):
      stringInfo[3] = splitLine[1]
    if(splitLine[0] == "NameOfSource"):
      stringInfo[0] = splitLine[1]
    if(splitLine[0] == "SrmPath"):
      stringInfo[4] = splitLine[1]
    if(splitLine[0] == "WorkPath"):
      stringInfo[5] = splitLine[1]
    if(splitLine[0] == "HPName"):
      stringInfo[6] = splitLine[1]

  file.close()

  for i in range(0,len(stringInfo)):
    if(stringInfo[i] == []):
      return "Error"   
  showMessage("Configuration Parameters:")
  showMessage("ReleasePath = " + stringInfo[1] + sanitycheck("ReleasePath",stringInfo[1]))
  showMessage("Executable = " + stringInfo[2] + sanitycheck("Executable",stringInfo[2]))
  showMessage("NameOfConf = " + stringInfo[3] + sanitycheck("NameOfConf",stringInfo[3]))
  showMessage("NameOfSource = " + stringInfo[0] + sanitycheck("NameOfSource",stringInfo[0]))
  showMessage("SrmPath = " + stringInfo[4])
  showMessage("WorkPath = " + stringInfo[5] + sanitycheck("WorkPath",stringInfo[5]))
  showMessage("HPName = " + stringInfo[6])
  return stringInfo

  
###############################################################
# showMessage method: Adds label and print message            #
#                                                             #
############################################################### 
def showMessage(message):
  print "[runManager] " + message


###############################################################
# List of tasks: from the source file                         #
#                                                             #
###############################################################
def getListOfTasks(nameOfFile):

  file = open(nameOfFile)
  
  tasks = []

  for line in file.readlines():
    if(line[0] == "#"):
      continue
    splitedLine = line.split()
    argumentBegin = line.find("'")
    argumentEnd = line.rfind("'")
	    
    if(argumentBegin == -1 or argumentEnd == -1 or argumentBegin == argumentEnd):
      showMessage("Ignoring the line : " + line)
    else:  
      arguments = line[argumentBegin+1:argumentEnd]
      endLine = line[argumentEnd+1:len(line)]
      splitedEndLine = endLine.split()
      dataset = splitedLine[0] 
      numberOfJobs = splitedEndLine[0]
      numberOfFilesPerJob = splitedEndLine[1]
      task = [dataset, arguments, numberOfJobs, numberOfFilesPerJob]
      if(numberOfJobs == '-1' or numberOfFilesPerJob == '-1'):
        tasks.append(task)

  return tasks  


###############################################################
# process: process a single task                              #
#                                                             #
###############################################################
def process(task, conf):
  dcapPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat"
  dbsUrl = 'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'
  list = [] # Initialize file list

  dataset = task[0]

  # Check if cache file exists and if we should use it
  cacheFile = '.'+dataset.replace('/','_').replace('?','_')
  allmyfiles=[]
  if not options.renew and os.path.isfile(cacheFile):
      print "Reading from cached file"
      f = open(cacheFile)
      showMessage('Reading files pertaining to '+dataset+'from cache file '+cacheFile)
      allmyfiles = f.readlines()
      f.close()
  else:
      # If not: rebuild the list
      showMessage("Going to fetch list of files pertaining to "+str(dataset))
      command = 'dbs search --query="find file where dataset='+dataset+'" --noheader --url='+dbsUrl
      print command
      theList = os.popen(command)
      list = theList.readlines()
      for li in list:
              allmyfiles.append(li)
      f = open(cacheFile,'w')
      f.write(''.join(allmyfiles))
  
  numberOfFiles = len(allmyfiles)
  if(numberOfFiles == 0):
    showMessage("No files found in "+str(task[0]))
    return "Error"

  correctList = [];
  for fileName in allmyfiles:
    auxiliar = dcapPath + fileName
    auxiliar = auxiliar[0:len(auxiliar)-1]
    correctList.append(auxiliar)

  if(task[2] != "-1"):
    FilesPerJob = int(float(numberOfFiles)/int(task[2]))
    OverloadedJobs = numberOfFiles + 1 - FilesPerJob * int(task[2])
    NumberOfJobs = int(task[2])
  else:
    FilesPerJob = int(task[3])
    NumberOfJobs = int(numberOfFiles/FilesPerJob) + 1
    OverloadedJobs = numberOfFiles + 1 - FilesPerJob * NumberOfJobs  
    
  if(numberOfFiles<int(task[2])) : 
    FilesPerJob = 1;
    NumberOfJobs = numberOfFiles
    OverloadedJobs = 0

  #FR FIXME: this part is tricky and might be fragile
  # Task[0] is now a dataset: /DATASET/USER-TAG_DATASET-TYPE-AOD-HASH/USER
  # Extract useful information (DATASET-TYPE-AOD)
  primaryDataset = dataset.lstrip('/').split('/')[0]
  fullName = dataset.lstrip('/').split('/')[1]
  index1 = fullName.index(primaryDataset)
  index2 = fullName.rindex('-')
  nameOfFolder = fullName[index1:index2]
  taskName = nameOfFolder
  os.system("mkdir -p " + taskName)
#  os.system("mkdir -p ../" + taskName)

  showMessage(str(NumberOfJobs) + " jobs with " + str(FilesPerJob) + " files each will be created")
  jobnumbercollection=[]
  if(FilesPerJob > 0):
    for step in range(0, NumberOfJobs):
      result=createJob(step, FilesPerJob, NumberOfJobs, OverloadedJobs, conf, correctList, nameOfFolder, task[1])
      jobnumbercollection.append(result)
  return jobnumbercollection

##################################################################################
#                              POST PROCESSING                                   #
##################################################################################

def cb(r): #optional: callback function
    for listentry in r:
      jobnumbers.append(listentry)
    pass

def ensure_dir(f) :
	if not os.path.exists(f):
		os.makedirs(f)

def join_directory(path,filelist,username) :
	localpath="/scratch/"+username+"/ntuples/"+path
	ensure_dir(localpath)
	cleanpath=path;
	if (cleanpath[len(cleanpath)-1]=="/") : # remove trailing slash
		cleanpath=cleanpath[0:len(cleanpath)-2]
	fusecommand="hadd -f /scratch/"+username+"/ntuples/"+cleanpath+".root "
	for item in filelist:
		copycommand="dccp dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/"+username+"/ProcessedTrees/"+item+" /scratch/"+username+"/ntuples/"+item
		commands.getstatusoutput(copycommand)
		fusecommand=fusecommand+" /scratch/"+username+"/ntuples/"+item
	print commands.getoutput(fusecommand)
	deletecommand="rm -r /scratch/"+username+"/ntuples/"+path+"/" 
	print commands.getoutput(deletecommand)
		

def check_directory(path,username) :
	complete_path="/pnfs/psi.ch/cms/trivcat/store/user/"+username+"/ProcessedTrees/"+path
	print "\033[1;34m Going to checkout the subdirectory "+complete_path+" \033[0m "
	listoffiles=[]
	supposedtobejoined=False;
	commandline="lcg-ls -l srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+complete_path
	pipe=popen(commandline)
	for l in pipe.readlines():
		currentline=l.strip("\n")
		if(currentline[0]=="d") :
			check_directory(currentline[currentline.find(path):],username)
		else :
			if(currentline.count(path) > 0) :
				supposedtobejoined=True
				listoffiles.append(currentline[currentline.find(path):])
	if supposedtobejoined==True:
	 	join_directory(path,listoffiles,username)


##################################################################################
#                              POST PROCESSING                                   #
##################################################################################

def sanitycheck(name,path):
	if os.path.exists(path) == False:
		showMessage("Problem detected with "+name+" which is set to "+path+" -- please fix this. Aborting.")
		sys.exit(-1)
	else :
		return " ... ok !"

##################################################################################
#                                  MAIN                                          #
##################################################################################
if __name__ == '__main__' :

        if domultiprocessing==True:
                po = Pool()
        fusepath=""
        uname=""
        isdata=0
        if len(args) != 1:
                parser.print_usage()
                sys.exit(-1)
        
        timeleft=commands.getoutput("timeleft=`voms-proxy-info -valid -timeleft | grep timeleft | awk '{ print $3 }'` && pos=`expr index "+'"$timeleft" :`&& timelefth=${timeleft:0:$pos-1} && echo $timelefth');
        if(timeleft<2): 
		print "You need to refresh your proxy! (will run voms-proxy-init -voms cms for you)"
		os.system("voms-proxy-init -voms cms");
	else:
		print "Proxy lifetime is acceptable (more than "+str(timeleft)+" hours)"
        
        result = parseInputFile(args[0])
        if(result == "Error"):
                showMessage("Error parsing input file")
                sys.exit(-1)

        listOfTasks = getListOfTasks(result[0])
        fusepath=result[4]
        uname=result[6]
        for l in listOfTasks:
                if(l[0].find("/data/")>-1) :
                        isdata=1
                showMessage("Processing " + l[0])
                if domultiprocessing==True:
                        po.apply_async(process,(l, result),callback=cb)
                else :
                        print "At this point you could be saving a lot of time with multiprocessing ... "
                        jobnumber = process(l, result)
                        jobnumbers.extend(jobnumber)
      		
        if domultiprocessing==True :
                po.close()
                po.join()

        if options.dryrun: sys.exit(0) # We're done for the dry run

        totaljobnumber=len(jobnumbers)
        counter=0
        print 'Total number of jobs:',totaljobnumber
        while(len(jobnumbers)>0 and counter<300) :
                time.sleep(60)
                counter+=1
                currlist=[]
                pipe=popen("qstat | grep `whoami` | awk '{print $1}'")
                for line in pipe.readlines():
                        currlist.append(int(line))
                checklist(jobnumbers,currlist)
        print "All jobs have finished - need to merge everything and place it in your scratch directory!"
        check_directory(fusepath,uname)
        if isdata==1 and result[2].find("RunJZBAnalyzer")>-1:
                print "We're dealing with data - we still need to merge data files and check for duplicates!"
                pipe=popen("/shome/buchmann/material/flash_remove_duplicates.exec -d /scratch/"+uname+"/ntuples/"+str(fusepath)+"/ -o /scratch/"+uname+"/"+fusepath+".root")
                for l in pipe.readlines():
                        print l
