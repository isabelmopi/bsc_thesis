#! /usr/bin/python

import os,sys
import subprocess

#only list here the routines that you have to compile yourself, all standard rivet routines do not have to be listed
ListRoutines = "RivetTEST_ROUTINE.so"  

OutputFlag = "H7_1k_1k5" # output flag cannot be too long, since the number of characters for the output on the grid is limited

# Change this to your username!
GridUser = "imolinos"

# Change this for testjob: nJobs = 2, nFilesPerJob = 5, otherwise larger (at least 10 files per job, 50 jobs)!
nJobs = 400
nFilesPerJob = 10

# here give the name of the containers you want to run over, and the cross-section of the sample (the ones below are just examples!!!)
ContList  = []
ContList.append(['user.aknue.703853.PhH7EG_H7UE_ttbar_hdamp258p75_bsf200_1k_1k5_SL__Bach3f_EXT0/', 2.73])


for Cont in ContList:
    Container = Cont[0]
    Xsec      = Cont[1]

    Sample     = Container.split(".")

    if "mc15" in Container or "mc16" in Container:
        NewCont    = Sample[1]+"."+Sample[2].replace("_EXT1/", "")+"_"+Sample[5].replace("/", "")
    else:
        NewCont    = Sample[2]+"."+Sample[3].replace("/", "")

    OutputFile = "Output_"+NewCont+".yoda"
     
     # now make JO script (python)
    submitFileNamePY = "JO_"+NewCont+".py"

    submitfile = open(submitFileNamePY, "w")
    codeLines = []
    codeLines.append("theApp.EvtMax = -1")
    codeLines.append("import AthenaPoolCnvSvc.ReadAthenaPool")
    codeLines.append("svcMgr.EventSelector.InputCollections =['bla']")
    codeLines.append("OutputYoda  = \""+OutputFile+"\"")
    codeLines.append("from AthenaCommon.AlgSequence import AlgSequence")
    codeLines.append("job = AlgSequence()")
    codeLines.append("from xAODEventInfoCnv.xAODEventInfoCnvConf import xAODMaker__EventInfoCnvAlg")
    codeLines.append("job += xAODMaker__EventInfoCnvAlg()")
    codeLines.append("from Rivet_i.Rivet_iConf import Rivet_i")
    codeLines.append("rivet = Rivet_i()")
    codeLines.append("rivet.AnalysisPath = os.environ['PWD']")    

    # here you have to give the list of rivet routines you want to use!!!
    codeLines.append("rivet.Analyses += ['TEST_ROUTINE'] ")
    codeLines.append("rivet.Analyses += ['ATLAS_2017_I1614149'] ")
    codeLines.append("rivet.Analyses += ['ATLAS_2022_I2037744'] ")

    codeLines.append("rivet.CrossSection = "+str(Xsec))
    codeLines.append("rivet.SkipWeights = True")
    codeLines.append("rivet.HistoFile = OutputYoda")
    codeLines.append("job += rivet")
    codeLines.append("from GaudiSvc.GaudiSvcConf import THistSvc")
    codeLines.append("svcMgr += THistSvc()") 

    for codeLine in codeLines:
        submitfile.write(codeLine+" \n")
        
    submitfile.close()

    os.system("pathena "+submitFileNamePY+" --extFile="+ListRoutines+" --nJobs="+str(nJobs)+" --nFilesPerJob="+str(nFilesPerJob)+" --extOutFile="+OutputFile+" --inDS="+Container+" --outDS=user."+GridUser+"."+NewCont+"_"+OutputFlag+"/ ")


