#! /usr/bin/env python

import glob
import os
import subprocess

manta='/hpc/cog_bioinf/ridder/users/cshneider/Manta-1.1.0_110517/manta-1.1.0.centos5_x86_64/bin/configManta.py'
#print manta
#'/hpc/local/CentOS7/cog_bioinf/manta_1.0.3/bin/configManta.py'
indir='/home/cog/cshneider/hpc_cshneider/Breast_Cancer_Pilot_outside_IAP'
#'/hpc/cog_bioinf/kloosterman/users/mroosmalen/ovarian/OC9/bams' # CHANGE THIS LINE
#print 'indir:', indir
genome='/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta'
#print genome

#Delly artifact
#chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']

path = os.path.dirname(os.path.realpath(__file__))
#print "path:", path

#tumor_bams = [] #list()
#normal_bams = [] #list()
#hold_jids = [] #list()

def qsub( threads, time, mem, mail, jobID, logdir, hold_jids, batch_file ):

	qsub = "qsub -pe threaded "+threads+" -l h_rt="+time+" -l h_vmem="+mem+" -m a -M "+mail+" -N "+jobID+" -o "+logdir+"/"+jobID+".log -e "+logdir+"/"+jobID+".err"
	if len(hold_jids) > 0:
		qsub += " -hold_jid "+",".join(hold_jids)
        #print "qsub:", qsub
	qsub += " "+batch_file
    #print "qsub:", qsub
	subprocess.call(qsub, shell=True)


def qsub_launcher(direct):

    tumor_bams = [] #list()
    normal_bams = [] #list()
    hold_jids = [] #list()

    direct = str(direct)
    #print 'direct:', direct

    #print "glob.glob(indir+direct+'/*.bam'):", glob.glob(indir+'/'+direct+'/*.bam')

    for bam in glob.glob(indir+'/'+direct+'/*.bam'):

    	filename = bam.split("/")[-1]
        #print 'filename:', filename

    	if filename.find("_B") != -1: #filename.startswith("B"):# CHANGE THIS LINE
    		normal_bams.append(bam)
    	else:
    		tumor_bams.append(bam)

    #print "normal_bams:", normal_bams
    #print "tumor_bams:", tumor_bams


    for normal in normal_bams:
    	n = normal.split("/")[-1].split("_")[0]
        #print "n:", n
    	for tumor in tumor_bams:
    		t = tumor.split("/")[-1].split("_")[0]
            #print "t:", t

    		cmd = manta+" --normalBam "+normal+" --tumorBam "+tumor+" --referenceFasta "+genome+" --runDir "+path+"/"+n+"_"+t+"\n\n"
            #print "cmd with config file:", cmd
    		cmd += path+"/"+n+"_"+t+"/runWorkflow.py -m local -j 8\n"
            #print "cmd with workflow script:", cmd

    		text = """
    #! /bin/bash

    echo `date` \": Running on \"`uname -n`

    """+cmd+"""

    echo `date`: Done
    """
    		jobID = n+"_"+t+"_manta"
            #print "jobID:", jobID
    		bashfile = path+"/"+jobID+".sh"
            #print "bashfile:", bashfile
    		bash = open( bashfile,'w' )
    		bash.write( text )
    		bash.close()
    		qsub( '8', '4:0:0', '10G', 'EMAIL@EMAIL', jobID, path, hold_jids, bashfile )


#Have manually formed this directory_list which contained all directories present in pwd w/o subdirectories; had to remove initial_prep manually
#find -maxdepth 1 -type d | cut -c3- > directory_list.txt #can automate this line if only have folders that strictly required
###directory_file = open(os.path.abspath('directory_list.txt'), 'r')
###print 'len(directory_file.readlines()):', len(directory_file.readlines())
###directory_file.close()

directory_file = open(os.path.abspath('directory_list.txt'), 'r')
#counter = 0
for direct in directory_file.readlines():
    direct = direct.strip()
    #print 'direct:', direct
    qsub_launcher(direct)
    #counter += 1
    #print counter
    #break
directory_file.close()
