import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
from chimera.tkgui import saveReplyLog

# change to folder with data files
os.chdir("/Users/wojciechpotrzebowski/Documents/masaxs/CaMComplex")
ref_names = [rn for rn in os.listdir("refs2") if rn.endswith(".pdb")]
file_names = [fn for fn in os.listdir("pdbs") if fn.endswith(".pdb")]

for ref in ref_names:
    os.chdir("/Users/wojciechpotrzebowski/Documents/masaxs/CaMComplex/refs2")
    rc("open " + ref)
    for fn in file_names:
        #replyobj.status("Processing " + fn) # show what file we're working on
        os.chdir("/Users/wojciechpotrzebowski/Documents/masaxs/CaMComplex/pdbs")
        rc("open " + fn)
        rc('matchmaker #0 #1 iter 15.0')
        rc('close #1')
    rc('close #0')
saveReplyLog("entire.log")
# uncommenting the line below will cause Chimera to exit when the script is done
#rc("stop now")
