# RivetTopStudents


## Login via command line to CERN (lxplus cluster)

```
ssh -Y USERNAME@lxplus.cern.ch
```


## Download the example from git:

```
git clone https://gitlab.cern.ch/aknue/RivetTopStudents.git
cd RivetTopStudents
```

The test code lives in TEST_ROUTINE.cc .


## Run code locally before making any changes

First step: read the code to understand what it does.

Second step: compile the code:
```
source setup.sh
source compile.sh
```

Third step: run the code locally (I have a testfile that you should be able to read)
```
athena.py JO_TEST.py
```

Fourth step: plotting the output:

First check if a .yoda file is now in the folder:
```
ls 
```
If yes, then run:
```
rivet-mkhtml TEST_PowhegPythia8.yoda
```

The plots live in subfolders in the folder "rivet-plots" .
They can be also seen as a webpage, maybe for now (and to practise the scp command) copy them to your local computer:

```
scp -r USERNAME@lxplus.cern.ch:PATH/TO/YOUR/FOLDER/rivet-plots . 
```
Obviously you have to change the USERNAME and the PATH/TO/YOUR/FOLDER to reflect your own username and folder.
If you are not sure what the full path of the folder is that you are in, just type "pwd".

If you just want to look at a plot quickly without copying everything over, you can just do:
```
display rivet-plots/TEST_ROUTINE/nJets.png
scp -r USERNAME@lxplus.cern.ch:PATH/TO/YOUR/FOLDER/rivet-plots .
```

## Now: Make first changes to the code

Your first task	would be to add	the event selection relevant for your analysis to the .cc code,	recompile and run again locally. 

ALWAYS test the code first locally on a test sample before you submit it to the grid later!!!


There are some existing rivet routines with unfolded data from ATLAS and CMS here that you can use as examples: https://rivet.hepforge.org/analyses.html
Each of the "ATLAS_" strings there points to one rivet routine, and by clicking on the link you can see the C++ code directly in the browser.


## After making the changes: put this code in your own git repository (not into mine, otherwise the three of you will interfere)!

## More plotting options can be found here:

https://gitlab.com/hepcedar/yoda/-/blob/release-2-1-x/doc/PlotConfig.md?ref_type=heads

If you implement a new plot and want to add a caption to the axes, please look at TEST_ROUTINE.plot to see how this works!

==============================================================================================================
==============================================================================================================

## Submission to the grid and more: will be added next week ;) 
```
```


## Samples to use when running on the grid (still under development)

The three of you need different samples for your studies, so ignore the ones that are not listed for you:

# David:
mc16_13TeV.560102.MGPy8EG_Toponium_1L.merge.EVNT.e8562_e8455
mc16_13TeV.802380.Py8EG_Toponium_2L.merge.EVNT.e8562_e8455
mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.merge.EVNT.e6337_e5984


# Mike:
mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.evgen.EVNT.e8514
mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.evgen.EVNT.e8514
mc23_13p6TeV.604482.PhPy8EG_A14_ttbar_hdamp258p75_recToTop_singleLep.evgen.EVNT.e8589
mc23_13p6TeV.604483.PhPy8EG_A14_ttbar_hdamp258p75_recoilToTop_dil.evgen.EVNT.e8589


# Isabel
mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.merge.EVNT.e6337_e5984
mc16_13TeV.601284.PhPy8EG_A14_ttbar_nospin_hdamp258p75_SingleLep.merge.EVNT.e8448_e7400
mc16_13TeV.601403.PhPy8EG_A14_ttbar_nospin_hdamp258p75_dil.merge.EVNT.e8448_e7400