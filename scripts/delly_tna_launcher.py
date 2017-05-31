#! /usr/bin/env python

import glob
import os
import subprocess
from os.path import realpath

delly='/hpc/cog_bioinf/ridder/users/cshneider/Delly-0.7.7_140517/delly_v0.7.7_parallel_linux_x86_64bit'
#/delly_v0.7.7_CentOS5.4_x86_64bit'
#'/hpc/local/CentOS7/cog_bioinf/delly_v0.7.6/delly'
indir='/home/cog/cshneider/hpc_cshneider/Breast_Cancer_Pilot_outside_IAP'
#'/hpc/cog_bioinf/kloosterman/users/mroosmalen/ovarian/OC9/bams'
genome='/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta'

chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']

#path = os.path.dirname(os.path.realpath(__file__))
path = os.path.dirname(os.path.realpath(__file__))+'/Delly'
print "path:", path

#tumor_bams = list()
#normal_bams = list()
svtypes = [ "DEL", "INS", "INV", "DUP", "BND" ]

def qsub( threads, time, mem, mail, jobID, logdir, hold_jids, batch_file ):

	qsub = "qsub -pe threaded "+threads+" -l h_rt="+time+" -l h_vmem="+mem+" -m a -M "+mail+" -N "+jobID+" -o "+logdir+"/"+jobID+".log -e "+logdir+"/"+jobID+".err"
	if len(hold_jids) > 0:
		qsub += " -hold_jid "+",".join(hold_jids)
	qsub += " "+batch_file

	subprocess.call(qsub, shell=True)


def qsub_launcher_delly(direct,subdir):

    direct = str(direct)
    print 'direct:', direct

    tumor_bams = [] #list()
    normal_bams = [] #list()

    print "glob.glob(indir+'/'+direct+'/*.bam'):", glob.glob(indir+'/'+direct+'/*.bam')

    for bam in glob.glob(indir+'/'+direct+'/*.bam'):
    	filename = bam.split("/")[-1]
    	if filename.find("_B") != -1: #if filename.startswith("B"):   #CHANGE HERE TOO!
    		normal_bams.append(bam)
    	else:
    		tumor_bams.append(bam)

    for normal in normal_bams:
    	n = normal.split("/")[-1].split(".")[0] #.split("_")[0]
        n_realp = realpath(normal).split("/")[-1].split("_")[0]
    	for tumor in tumor_bams:
    		t = tumor.split("/")[-1].split(".")[0] #.split("_")[0]
            t_realp= realpath(tumor).split("/")[-1].split("_")[0]
    		tsvfile = path+"/"+n+"_"+t+'.tsv'
    		tsv = open( tsvfile, 'w' )
    		tsv.write( n_realp+"\t"+'control'+"\n")
    		tsv.write( t_realp+"\t"+'tumor' )
    		tsv.close()

    		for svtype in svtypes:
    			if svtype == 'BND':
    				for chr in chrs:
    					for chr2 in chrs:
    						if chrs.index(chr2) <= chrs.index(chr):
    							continue

    						excl = open( path+"/"+chr+"-"+chr2+".excl",'w' )
    						for chr3 in chrs:
    							if chr3 != chr and chr3 != chr2:
    								excl.write(chr3+"\n")
    						excl.close()

    						hold_jids = [] #list()

    						cmd = delly+" call -t "+svtype+" -x "+path+"/"+chr+"-"+chr2+".excl"+" -o "+path+"/"+n+"_"+t+"_"+svtype+"_"+chr+"-"+chr2+".bcf -g "+genome+" "+normal+" "+tumor
    						text = """
    #! /bin/bash

    echo `date` \": Running on \"`uname -n`

    """+cmd+"""

    echo `date`: Done
    """
    						jobID = n+"_"+t+"_"+svtype+"_"+chr+"-"+chr2+"_call"

    						bashfile = path+"/"+jobID+".sh"
    						bash = open( bashfile,'w' )
    						bash.write( text )
    						bash.close()

    						qsub( '8', '4:0:0', '10G', 'EMAIL@EMAIL', jobID, path, hold_jids, bashfile )
    						hold_jids.append(jobID)

    						jobID = n+"_"+t+"_"+svtype+"_"+chr+"-"+chr2+"_filter"

    						cmd = delly+" filter -t "+svtype+" -f somatic -o "+path+"/"+subdir+"/"+n+"_"+t+"_"+svtype+"_"+chr+"-"+chr2+".pre.bcf -s "+tsvfile+" "+path+"/"+n+"_"+t+"_"+svtype+"_"+chr+"-"+chr2+".bcf"
    						text = """
    #! /bin/bash

    echo `date` \": Running on \"`uname -n`

    """+cmd+"""

    echo `date`: Done
    """
    						bashfile = path+"/"+jobID+".sh"
    						bash = open( bashfile,'w' )
    						bash.write( text )
    						bash.close()

    						qsub( '8', '3:0:0', '10G', 'EMAIL@EMAIL', jobID, path, hold_jids, bashfile )
    						hold_jids.append(jobID)

    			else:
    				for chr in chrs:
    					excl = open( path+"/"+chr+".excl",'w' )
    					for chr2 in chrs:
    						if chr2 != chr:
    							excl.write(chr2+"\n")
    					excl.close()

    					hold_jids = [] #list()


    					cmd = delly+" call -t "+svtype+" -x "+path+"/"+chr+".excl"+" -o "+path+"/"+n+"_"+t+"_"+svtype+"_"+chr+".bcf -g "+genome+" "+normal+" "+tumor
    					text = """
    #! /bin/bash

    echo `date` \": Running on \"`uname -n`

    """+cmd+"""

    echo `date`: Done
    """
    					jobID = n+"_"+t+"_"+svtype+"_"+chr+"_call"

    					bashfile = path+"/"+jobID+".sh"
    					bash = open( bashfile,'w' )
    					bash.write( text )
    					bash.close()

    					qsub( '8', '4:0:0', '10G', 'EMAIL@EMAIL', jobID, path, hold_jids, bashfile )
    					hold_jids.append(jobID)

    					jobID = n+"_"+t+"_"+svtype+"_"+chr+"_filter"

    					cmd = delly+" filter -t "+svtype+" -f somatic -o "+path+"/"+subdir+"/"+n+"_"+t+"_"+svtype+"_"+chr+".pre.bcf -s "+tsvfile+" "+path+"/"+n+"_"+t+"_"+svtype+"_"+chr+".bcf"
    					text = """
    #! /bin/bash

    echo `date` \": Running on \"`uname -n`

    """+cmd+"""

    echo `date`: Done
    #"""
    					bashfile = path+"/"+jobID+".sh"
    					bash = open( bashfile,'w' )
    					bash.write( text )
    					bash.close()

    					qsub( '8', '3:0:0', '10G', 'EMAIL@EMAIL', jobID, path, hold_jids, bashfile )
    					hold_jids.append(jobID)

directory_file = open(os.path.abspath('directory_list.txt'), 'r') #so this script is in the same directory as this txt file.
counter = 0
for direct in directory_file.readlines():
    direct = direct.strip()
    print 'direct:', direct
    subdir = direct+"_Delly_filtered2_SVs" #n+"_"+t
    print subdir
    print "mkdir" + " " + path + "/" + subdir
    os.system("mkdir" + " " + path + "/" + subdir)
    qsub_launcher_delly(direct,subdir)
    counter += 1
    print counter
    break #NOTEE: #when launch with the break, need to delete first entry, MMC01021, from directory_list.txt!!!
directory_file.close()

#"""insert break into function run to make sure that runs only on one patient initially and turn off qsub for now!"""
