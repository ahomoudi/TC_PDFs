#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##################################################################################################
# Routine to download in a systematically way the CORDEX data from
# https://esgf-data.dkrz.de/search/cordex-dkrz/ 
#
# @author: Ana Aguilar - ana.aguilar@mailbox.tu-dresden.de
# September 2020
##################################################################################################

'''==================================== Load Libraries ======================================='''

import os #for comunicating with the bash
#import subprocess #for comunicating with the bash
import shutil #to move file
import glob

from pyesgf.search import SearchConnection
from pyesgf.logon import LogonManager 
from myproxy.client import MyProxyClient #pip install myproxyclient from console
from OpenSSL import SSL
MyProxyClient.SSL_METHOD = SSL.TLSv1_2_METHOD


'''================================= Preparations ============================================='''

#Set Working Directory
#src = '/mnt/Shared/000Documentos/001_Master_HSE/Thesis/CORDEX/Wget/'
src = '/scratch/ws/1/anag999a-DownPTY/CORDEX/'
os.chdir(src) 

conn = SearchConnection('http://esgf-data.dkrz.de/esg-search', distrib=False) #connect to a server

# ===========================================================================================
# The first time this need to be run
# Note: This never worked in the HPC, so I ran it in my computer and then copied the entire 
# content of my ~/esg folder to the same folder in the home in the HPC
#
#OPENID = 'https://esgf-data.dkrz.de/esgf-idp/openid/anaguilarc'
#lm.logon_with_openid(openid=OPENID, password=None, bootstrap=True)
#lm.is_logged_on()
#
# ===========================================================================================

#loggin to ESGF
lm = LogonManager()
#lm.logoff()
#lm.logon_with_openid('https://esgf-data.dkrz.de/esgf-idp/openid/XXXXX', password= 'XXXXXXXX') 
#lm.is_logged_on()


var6 = ['ta300'
        ,'ta500'
        ,'ta700'
        ,'ta850'
        ,'ua850'
        ,'ua300'
        ,'va850'
        ,'va300'
        ,'uas'
        ]

var3 = ['pr'
        ,'psl'
        ]

'''================================= WGET - Bash creation ====================================='''
# ===========================================================================================
# The following functions will:
#    1. Create a list of variables that need to be downloaded, depending on the frequency 
#       available in CORDEX under the rest of the giving constrains
#    2. Generate the WGET scrips
#    3. Generate the txt scrip for downloding in the HPC (can be done, but not necessary now)
#
# ===========================================================================================

def Wget_hist_6hrs():
    for v in var6:
        ctx = conn.new_context(
            project="CORDEX", 
            institute="CLMcom-ETH",
            driving_model="NCC-NorESM1-M",
            experiment="historical", 
            time_frequency="6hr",
            domain="WAS-22", 
            ensemble="r1i1p1", 
            variable= v
            )
        hits = ctx.hit_count
        for i in range(hits):
            ds = ctx.search()[i]
            fc = ds.file_context()
            wget_script_content = fc.get_download_script()
            script_path = src + 'wget_hist_'+ v + '.sh'
            with open(script_path, "w") as writer:
                writer.write(wget_script_content)
            #script_path2 = src + 'wget_'+ v + '.txt'
            #HPC_script = (
'''#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH -J "CORDEX"   # job name
#SBATCH --mail-user=ana.aguilar@mailbox.tu-dresden.de   # email address
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT,TIME_LIMIT_90
#SBATCH -A p_climate_panama
#SBATCH --output=outp-m-%j.out
#SBATCH --error=erro-m-%j.out

srun bash'''
           # )
           # with open (script_path2, 'w+') as rsh:
            #    rsh.write(HPC_script + ' ' + script_path + ' -H -d ')

def Wget_hist_3hrs():
    for v in var3:
        ctx = conn.new_context(
            project="CORDEX", 
            institute="CLMcom-ETH",
            driving_model="NCC-NorESM1-M",
            experiment="historical", 
            time_frequency="3hr",
            domain="WAS-22", 
            ensemble="r1i1p1", 
            variable= v
            )
        hits = ctx.hit_count
        for i in range(hits):
            ds = ctx.search()[i]
            fc = ds.file_context()
            wget_script_content = fc.get_download_script()
            script_path = src + 'wget_hist_'+ v + '.sh'
            with open(script_path, "w") as writer:
                writer.write(wget_script_content)
            
def Wget_RCP85_6hrs():
    for v in var6:
        ctx = conn.new_context(
            project="CORDEX", 
            institute="CLMcom-ETH",
            driving_model="NCC-NorESM1-M",
            experiment="rcp85", 
            time_frequency="6hr",
            domain="WAS-22", 
            ensemble="r1i1p1", 
            variable= v
            )
        hits = ctx.hit_count
        for i in range(hits):
            ds = ctx.search()[i]
            fc = ds.file_context()
            wget_script_content = fc.get_download_script()
            script_path = src + 'wget_rcp85_'+ v + '.sh'
            with open(script_path, "w") as writer:
                writer.write(wget_script_content)
                
def Wget_RCP85_3hrs():
    for v in var3:
        ctx = conn.new_context(
            project="CORDEX", 
            institute="CLMcom-ETH",
            driving_model="NCC-NorESM1-M",
            experiment="rcp85", 
            time_frequency="3hr",
            domain="WAS-22", 
            ensemble="r1i1p1", 
            variable= v
            )
        hits = ctx.hit_count
        for i in range(hits):
            ds = ctx.search()[i]
            fc = ds.file_context()
            wget_script_content = fc.get_download_script()
            script_path = src + 'wget_rcp85_'+ v + '.sh'
            with open(script_path, "w") as writer:
                writer.write(wget_script_content)

            
            
Wget_hist_6hrs()
Wget_hist_3hrs()
Wget_RCP85_6hrs()
Wget_RCP85_3hrs()

'''======================================== Download ==========================================='''
#list the nc files
files = glob.glob('*.sh')

#downloadn them
for f in files:
    os.system('bash '+ f)
    

#lm.logoff()
    
'''=================================== Post Processing ========================================='''    

#Delete files before 1990 and after 2050
fileH = list(range(1950,1990))
fileR = list(range(2050,2100))

for h in fileH:
    pattern = glob.glob('*'+str(h)+'*.nc')
    for p in pattern:
        os.remove(p)

for r in fileR:
    pattern = glob.glob('*'+str(r)+'*.nc')
    for p in pattern:
        os.remove(p)        

#Move historical and rcp8.5 to its respective folders
hist = glob.glob('*historical*.nc')
rcp = glob.glob('*rcp85*.nc')

for h in hist:
    shutil.move(src + h,
                src + 'Historical/' + h)

for r in rcp:
    shutil.move(src + r,
                src + 'RCP_85/' + r)    

#Create folders from 1990-2049 and Move files to its yearly folder
folderhist = src + 'Historical/'
foldercp = src + 'RCP_85/'
folder1 = list(range(1990,2006))
folder2 = list(range(2006,2050))


os.chdir(folderhist)
for f in folder1:
    os.mkdir(str(f))
       
for f in folder1:
    for y in glob.glob('*'+str(f)+'123118.nc'):
        if str(f) in y:
            shutil.move(y, src + 'Historical/' + str(f) + '/')
 
os.chdir(foldercp)
for f in folder2:
    os.mkdir(str(f))

for f in folder2:
    for y in glob.glob('*'+str(f)+'123118.nc'):
        if str(f) in y:
            shutil.move(y, src + 'RCP_85/' + str(f) + '/')
