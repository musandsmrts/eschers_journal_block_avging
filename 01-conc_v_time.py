#!/usr/bin/env python

import math
import numpy as np
import matplotlib
import time
import sys # for sys.exit()
matplotlib.use('agg')

params = {'mathtext.default': 'regular',
          'font.family': 'Times New Roman'
         }
matplotlib.rcParams.update(params)

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import gridspec

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
gs.update(wspace=0,hspace=0)
ax = fig.add_subplot(gs[0,0])

################################################################################
# Set testing state ############################################################
################################################################################
#  0 : full data set  ##########################################################
#  1 : 4 data points  ##########################################################
################################################################################

test_state = 0

################################################################################
# Load data

print('Loading data...')
if test_state==0:
    data01 = np.loadtxt('tvDvZ_-120_120_D1CSVR_20200930.dat')
elif test_state==1:
    data01 = np.loadtxt('tvDvZ_testset.004.dat')
else:
    print('Invalid test state')
    sys.exit()
print('Data loaded')
################################################################################

sampleCnt = len(data01)
data_len = len(data01[:,0])

# We have a dataset with N samples
print('Let\'s say we have a data set, with N={} samples.'.format(data_len))

# Calculate average concentration at Z=-30 Angstrom
mean_neg30 = np.mean(data01[1:sampleCnt,91])
#                             ^ Upper bound is up-to but NOT including.
#                               If it was inclusive, sampleCnt-1.
print('We can look at its sample mean.')
print('Sample mean [O2] at {} Ang. = {}'.format(data01[0,91],mean_neg30))

# Calculate sample standard deviation at Z=-30 Angstrom
print('But does that tell the whole story?')
std_neg30 = np.std(data01[1:sampleCnt,91])
print('The sample standard deviation [O2] at {} Ang. = {}'.format(data01[0,91],std_neg30))

# Naively calculate the standard error of the mean (SEM)
naive_sem_neg30 = std_neg30/math.sqrt(sampleCnt-1)
print('Naively we can calculate a standard error of mean: {}'.format(naive_sem_neg30))

####
# Plot [O2] at Z=-30 Angstrom vs time
####
ax.plot(data01[1:sampleCnt,0]/100,data01[1:sampleCnt,91],alpha=0.85,color='r',linewidth=0.5)
ax.plot(data01[1:sampleCnt,0]/100,np.full(sampleCnt-1,mean_neg30),alpha=0.85,color='m',linewidth=0.5)

# Title
ax.set_title('[O$_2$] at Z=-30$\AA{}$ vs Time (M?)',fontname='Times New Roman',fontweight='bold',fontsize=14)

# X Axis Styling
ax.set_xlabel('Time (ns)', fontname='Times New Roman', fontweight='bold', fontsize=12)

ax.tick_params(axis='x',labelsize=12)
for tick in ax.get_xticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

ax.set_xlim([0,1128.23])

# Y Axis Styling
ax.set_ylabel('[O$_2$] at Z=-30', fontname='Times New Roman', fontweight='bold', fontsize=14)

#plt.yticks(np.arange(-100,101,10))
ax.tick_params(axis='y',labelsize=10)
for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

ax.set_ylim([-0.001,0.2])
  
# Grid
plt.grid(1,which='both',dashes=(1,1),linewidth=0.5)

fig.set_size_inches(7,4)
fig.savefig('O_2_conc_v_time_Z=-30.png',dpi=500,transparent=True,bbox_inches='tight')

################################################################################
# Plot [O2] at Z=-30 Angstrom vs time with blocks containing two frames
####

# Calculate blocks
blk_2 = np.full((data_len-1)/2,-1.)

print('Initial blk_2 values: {} {} ...'.format(blk_2[0],blk_2[1]))

for it1 in range(0,(data_len-1)/2):
    blk_2[it1] = np.mean(data01[(it1*2)+1:(it1*2)+1+1+1,91])
#                                     ^       ^ ^ ^ exclusive range
#                                     |       | |
#                                     |       | +-- for block size 2
#                                     |       | 
#                                     +-------+-- +1 because header
    print('{} {}'.format(blk_2[it1],np.mean(data01[(it1*2)+1:(it1*2)+1+1+1,91])))
    #print(blk_2[it1])

print('Computed blk_2 values: {} {} ...'.format(blk_2[0],blk_2[1]))

# Create time data for blocks
time_2 = np.full((data_len-1)/2,-1)
for it1 in range(0,len(time_2)):
    time_2[it1] = 2*it1

print('time_2 values: {} {} ...'.format(time_2[0],time_2[1]))

# Clear figure
plt.cla()

# Plot figure

ax.plot(time_2/100.,blk_2,alpha=0.85,color='g',linewidth=0.5)

# Title
ax.set_title('[O$_2$] at Z=-30$\AA{}$ vs Time Blocksize=2',fontname='Times New Roman',fontweight='bold',fontsize=14)

# X Axis Styling
ax.set_xlabel('Time (ns)', fontname='Times New Roman', fontweight='bold', fontsize=12)

ax.tick_params(axis='x',labelsize=12)
for tick in ax.get_xticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

ax.set_xlim([0,1128.23])

# Y Axis Styling
ax.set_ylabel('[O$_2$] at Z=-30', fontname='Times New Roman', fontweight='bold', fontsize=14)

ax.tick_params(axis='y',labelsize=10)
for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

ax.set_ylim([-0.001,0.2])
  
# Grid
plt.grid(1,which='both',dashes=(1,1),linewidth=0.5)

fig.set_size_inches(7,4)
fig.savefig('O_2_conc_v_time_Z=-30.2blk.png',dpi=500,transparent=True,bbox_inches='tight')

################################################################################
# Average from full trajectory
####

# Initialize array of averages
z_slices = len(data01[0,:])-1
print(z_slices)

avgs = np.zeros(z_slices)

for z0 in range(0,z_slices):
    avgs[z0] = np.mean(data01[1:data_len,z0+1])

# Clear figure
plt.cla()

ax.plot(avgs[:],data01[0,1:z_slices+1],alpha=0.85,color='r',linewidth=1.,label='z=0 avg')

# Title

ax.set_title('O2 Conc. vs Z for 5 $O_2$ D1CSVR',fontname='Times New Roman',fontweight='bold',fontsize=14)

# X Axis Styling

ax.set_xlabel('[O2]', fontname='Times New Roman', fontweight='bold', fontsize=12)

ax.tick_params(axis='x',labelsize=12)
for tick in ax.get_xticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

plt.xscale('log')
ax.set_xlim([0.00025,0.042])

# Y Axis Styling

ax.set_ylabel('Z', fontname='Times New Roman', fontweight='bold', fontsize=14)

plt.yticks(np.arange(-100,101,10))
ax.tick_params(axis='y',labelsize=10)
for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

ax.set_ylim([-102.,102.])
  
# Grid
plt.grid(1,which='both',dashes=(1,1),linewidth=0.5)

fig.set_size_inches(4,7)
fig.savefig('O2conc_v_Z.png',dpi=500,transparent=True,bbox_inches='tight')

sys.exit()
