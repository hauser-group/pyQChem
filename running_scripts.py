# pyQchem - Input/Output-Tools for Q-Chem
# Copyright (C) 2014  Matthew Goldey
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.

#####################################################################
#                                                                   #
#                 pyQchem - Scripts for running Q-Chem              #
#                                                                   #
#####################################################################

# modified my AWH Feb 2014: global method 'run()' is now hidden to reduce confusion,
# inputfile and multifile objects have their own 'run()' method which calls '_run()'.

from input_classes import *

def _run(inp_file,name='',loc53='',qchem='',nt=1,np=1,timestamp=False):
    """This is a script for running Q-Chem inside of iPython given an input file object, name, 53.0 file, location, and number of threads or processors"""
    #tag it all with current time for safety's sake
    import time,os
    curtime=time.strftime("%Y%m%d%H%M%S")
    if name=='':
        name=inp_file.runinfo.name
    if name=='':
        name=curtime
    if timestamp!=True:
        name=name.replace(".in","")
    else: 
        if name.endswith(".in"):
            name.replace(".in",curtime)
        else:
        	name=name+curtime
    if loc53=='':
        loc53=inp_file.runinfo.loc53
    if qchem=='':
        qchem=inp_file.runinfo.qchem
    if nt==1:
        nt=inp_file.runinfo.nt
    if np==1:
        np=inp_file.runinfo.np
    


    #make script file
    scr=name+".sh"
    scr_out=open(scr,'w')

    #source appropriate Q-Chem
    if qchem!='':
        scr_out.write("source ~/"+qchem+"\n")

    #copy 53.0 if restarting
    if loc53!='':
        inp_file.rem.scf_guess('read')
        scr_out.write("mkdir $QCSCRATCH/"+name+".dir\n")
        if loc53.endswith(".in.53.0"):
            scr_out.write('cp '+loc53+' $QCSCRATCH/'+name+'.dir/53.0\n')
        else:
            scr_out.write('cp '+loc53+'/53.0 $QCSCRATCH/'+name+'.dir/\n')

    inp_file.write(name+".in")

    #write qchem command to script file
    if nt>1:
    	scr_out.write("qchem -nt "+str(nt)+" "+name+".in "+name+".out "+name+".dir "+"\n")
    elif np>1:
    	scr_out.write("qchem -np "+str(nt)+" "+name+".in "+name+".out "+name+".dir "+"\n")
    else:
    	scr_out.write("qchem "+name+".in "+name+".out "+name+".dir "+"\n")

    #close and run
    scr_out.close()
    os.popen("bash "+scr).read()
    return

def queue(joblist,num_workers=1):
	"""This is a simple queue for running through a list of jobs (based on the content of each runinfo object).  
	When in doubt, set num_workers to the number of cores on the machine.
	Advanced options are currently not supported."""
	import Queue
	import threading
	q_in = Queue.Queue(maxsize=0)
	q_out = Queue.Queue(maxsize=0)
	# process that each worker thread will execute until the Queue is empty
	def worker():
	    while True:
	        # get item from queue, do work on it, let queue know processing is done for one item
	        item = q_in.get()
	        _run(item)
	        q_out.put(item.runinfo.name)
	        q_in.task_done()

	# another queued thread we will use to print output
	def printer():
	    while True:
	        # get an item processed by worker threads and print the result. Let queue know item has been processed
	        item = q_out.get()
	        print "Completed ", item
	        q_out.task_done()

	# launch all of our queued processes
	# Launches a number of worker threads to perform operations using the queue of inputs
	for i in range(num_workers):
	     t = threading.Thread(target=worker)
	     t.daemon = True
	     t.start()

	# launches a single "printer" thread to output the result (makes things neater)
	t = threading.Thread(target=printer)
	t.daemon = True
	t.start()

	# put items on the input queue
	for item in joblist:
		q_in.put(item)

	# wait for two queues to be emptied (and workers to close)
	q_in.join()       # block until all tasks are done
	q_out.join()

	print "Processing Complete"