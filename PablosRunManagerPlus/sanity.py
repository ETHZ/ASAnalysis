#
# Sanity checks on standard output / standard error 
#
# For now just checks that all files are processed and number of events match


import sys,re
import os
from os import path
import commands
import time


from optparse import OptionParser
parser = OptionParser(usage="%prog <directory>")
parser.add_option("-v","--verbose",dest="verbosity",type="int",default=0,
                  help="Verbosity level")
parser.add_option("--doDBS",dest="doDBS",action="store_true",default=False,
                  help="Check DBS info as well")
parser.add_option("--noPrintErr",dest="noPrintErrors",action="store_true",default=False,
                  help="Print count of errors found in error output")
parser.add_option("--err-only",dest="doErrOnly",action="store_true",default=False,
                  help="Only check for suspicious errors")
(options, args) = parser.parse_args(args=None, values=None)

# Send a query to dbs
def dbs_query(query):
    dbsUrl = 'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'

    command = 'dbs search --query="'+query+'" --noheader --url='+dbsUrl
    if options.verbosity>2: print "DBS query:",command
    return commands.getoutput(command)

# Process error logs
def parse_errors(efiles):
    errors = dict()
    allOK = True
    badTransfers = []
    # skip lcg-cp
    lcgPat1 = re.compile('Source URL')
    lcgPat2 = re.compile('bytes.*avg.*inst')
    lcgPat3 = re.compile('File size:')
    lcgPat4 = re.compile('Destination [URL|SRM]')
    lcgPat5 = re.compile('Checksum type')
    lcgPat6 = re.compile('Using grid catalog')
    lcgPat7 = re.compile('# streams:')
    lcgPat8 = re.compile('VO name:')
    lcgPatOK = re.compile('Transfer took')
    minuitPat = re.compile('Info in <Minuit2>')
    for filename in efiles.values():
        with open(filename) as f:
            transferOK = False
            for line in f:
                line = line.rstrip()
                m = re.search(lcgPat1,line)
                if m: continue
                m = re.search(lcgPat2,line)
                if m: continue
                m = re.search(lcgPat3,line)
                if m: continue
                m = re.search(lcgPat4,line)
                if m: continue
                m = re.search(lcgPat5,line)
                if m: continue
                m = re.search(lcgPat6,line)
                if m: continue
                m = re.search(lcgPat7,line)
                if m: continue
                m = re.search(lcgPat8,line)
                if m: continue
                m = re.search(lcgPatOK,line)
                if m:
                    transferOK = True
                    continue
                m = re.search(minuitPat,line)
                if m: line = 'Info in <Minuit2>'
                if line in errors:
                    errors[line] += 1
                else:
                    errors[line] = 1
        if not transferOK: badTransfers.append(filename)
    # Parse errors
    if not options.noPrintErrors:
       print '   Errors found:'
    for (k,v) in errors.iteritems():
       if not options.noPrintErrors:
          print '    {0:10d} times {line}'.format(v,line=k)

    # Transfer errors
    if (not options.doErrOnly):
       if len(badTransfers)>0:
          print ' *** Bad transfers:'
          print badTransfers
          allOK = False
       else:
          print '   Transfers: OK'


    return not allOK
    
                

# Process task directory and check a few things
def process_dir(dir):
    errors = 0
    
    # Get list of output and scripts
    oFiles = dict()
    eFiles = dict()
    sFiles = dict()
    oPat = re.compile('job_(\d+)\.out')
    ePat = re.compile('job_(\d+)\.err')
    sPat = re.compile('myScript_(\d+)\.sh')
    for f in os.listdir(dir):
        m = re.search(oPat,f)
        if m: oFiles[int(m.group(1))] = path.join(dir,f)
        m = re.search(ePat,f)
        if m: eFiles[int(m.group(1))] = path.join(dir,f)
        m = re.search(sPat,f)
        if m: sFiles[int(m.group(1))] = path.join(dir,f)

    if options.verbosity>1: print 'Found:',len(oFiles),'output +',len(eFiles),'error +',len(sFiles),' script files'

    # Get list of actually processed data files
    procInputFiles = dict()
    dcapPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat"
    ofilePat = re.compile('Adding file: '+dcapPath+'(.*\.root)')
    prodEventsPat = re.compile('..\ pairs\s+(\d+)\s+')
    procEventsPat = re.compile('total events.*=\s+(\d+)')
    totProdEvents = 0
    for k,v in oFiles.iteritems():
        with open(v) as f:
            prodEvents = 0
            procEvents = 0
            inputFile = ''
            for line in f:
                m = re.search(ofilePat,line)
                if m:
                    inputFile = m.group(1)
                m = re.search(procEventsPat,line)
                if m:
                    procEvents = m.group(1)
                m = re.search(prodEventsPat,line) # This is last
                if m:
                    prodEvents = int(m.group(1))
                    totProdEvents = totProdEvents + prodEvents
                    tuple = [k,procEvents,prodEvents]
                    procInputFiles[inputFile] = tuple

    if options.verbosity>1:
        print 'Found',len(procInputFiles),'processed input files, total',totProdEvents,' events produced'
    if options.verbosity>2:
        for k,v in procInputFiles.iteritems():
            print 'job',k,': ',v
        
    ### Need to check they were actually processed (number of events?)

    # Get list of data files in the script files
    scriptInputFiles = dict()
    ofilePat = re.compile('^SOURCEFILES=\"(.*?)\s*\"')
    for k,v in sFiles.iteritems():
        with open(v) as f:
            for line in f:
                m = re.search(ofilePat,line)
                if m:
                    filestring = m.group(1).rstrip().strip()
                    filelist = m.group(1).split(' ')
                    for ifile in filelist:
                        ifile = ifile.replace('dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat','')
                        scriptInputFiles[ifile] = k
    if options.verbosity>1: print 'Found',len(scriptInputFiles),'input files in scripts'

    if not options.doErrOnly:
       if len(procInputFiles) != len(scriptInputFiles):
          print ' *** script and processsed input files do not match! (',len(procInputFiles),',',len(scriptInputFiles),')'
          errors+=1
       else:
          print '   Processed vs. script files count: OK'

    # Stop here if there's nothing to process
    if len(scriptInputFiles)==0: return
    
    # Get dataset name: use DBS
    if options.doDBS:
        ofile = scriptInputFiles.keys()[0]
        query = 'find dataset where file='+ofile
        dataset = dbs_query(query)
        if options.verbosity>0: print 'Dataset is',dataset

        # Get file info from DBS
        dbsFiles = dict()
        query = 'find file,file.numevents where dataset='+dataset
        dbs_result = dbs_query(query)
        for line in dbs_result.split('\n'):
            (file,events) = line.split()
            dbsFiles[file] = int(events)
        if options.verbosity>1: print 'Found',len(dbsFiles),'from DBS for dataset',dataset
        if len(dbsFiles) != len(scriptInputFiles):
            print ' *** script and DBS input files do not match!'
            errors += 1
        else:
            print '   DBS vs. script files count: OK'
    
    errors += parse_errors(eFiles)


    if errors>0:
        print '*** PROBLEMS WITH',directory
    else:
        print '--> ALL OK'

##################################################################################
#                                  MAIN                                          #
##################################################################################
if __name__ == '__main__' :
    if len(args)<1:
        parser.print_usage()
        sys.exit(-1)

    for directory in args:
        print '== Processing ==',directory
        process_dir(directory)
