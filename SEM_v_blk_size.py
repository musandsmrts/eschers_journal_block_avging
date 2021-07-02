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

print('Loading data...')
data01 = np.loadtxt('01-tvDvZ_-120_120_D1CSVR_20200930.dat')
print('Data loaded')

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

####
# Plot [O2] at Z=-30 Angstrom vs time with blocks containing two frames
####

# Calculate blocks
blck_2 = np.full((data_len-1)/2,-1)

print('blck_2 values: {} {} {}...'.format(blck_2[0],blck_2[1],blck_2[2]))

# Create time data for blocks
time_2 = np.full((data_len-1)/2,-1)
for it1 in range(0,len(time_2)):
    time_2[it1] = 2*it1

print('time_2 values: {} {} {}...'.format(time_2[0],time_2[1],time_2[2]))


# Clear figure
plt.cla()

# Plot figure

ax.plot(time_2/100.,blck_2,alpha=0.85,color='g',linewidth=0.5)

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

#plt.yticks(np.arange(-100,101,10))
ax.tick_params(axis='y',labelsize=10)
for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

ax.set_ylim([-1.1,0.2])
  
# Grid
plt.grid(1,which='both',dashes=(1,1),linewidth=0.5)

fig.set_size_inches(7,4)
fig.savefig('O_2_conc_v_time_Z=-30.2blk.png',dpi=500,transparent=True,bbox_inches='tight')

sys.exit()
################################################################################
################################################################################

# Once the script is deemed complete, delete remaining code.

# Initialize array of averages
z_slices = len(data01[0,:])-1
print(z_slices)

avgs = np.zeros(z_slices)
print('avgs')
print(avgs)
print('data01[0,:]')
print(data01[0,:])
print('data01[0,1:z_slices+1]')
print(data01[0,1:z_slices+1])
#               ^ the colon range operator seems exclusive at the upper limit 
#                  so we need to add +1


# Block Averaging
blk_size = 10000 # 100 ns, 100 frms/ns
#blk_size = 5000 # 50 ns, 100 frms/ns
blk_lim = (data_len-1)/blk_size
print('Block size: {}'.format(blk_size))
print('Data length: {}'.format(data_len-1))
print('Block limit: {}'.format(blk_lim))
blk_means = np.zeros((z_slices,blk_lim))
blk_avgs = np.zeros((z_slices,))
conc_blk_SEMs = np.zeros((z_slices,))
print(blk_means.shape)

for z0 in range(0,z_slices):
    for blk in range(0,blk_lim):
        lim_up = blk_size*(blk+1)
        lim_lo = blk_size*blk

	blk_means[z0,blk] = np.mean(data01[lim_lo+1:lim_up+1,z0+1])

    print('blk STD:{} {}'.format(z0,np.std(blk_means[z0])))
    conc_blk_SEMs[z0] = np.std(blk_means[z0])/math.sqrt(blk_lim)

    #avgs[z0] = np.mean(data01[1:blk_size*blk_lim+1,z0+1])
     # For comparison to blk_avg. For identical data set, identical averages.
    blk_avgs[z0] = np.mean(blk_means[z0])
    avgs[z0] = np.mean(data01[1:data_len,z0+1])
#                             ^ we start at the 1st row (not the 0th),
#                               discarding the header row, which contains the
#                               z-value of the slice

    if 0: #For print statements
        print('len(data01[:,z0+1])')
        print(len(data01[:,z0+1]))
        print('data01[:,z0+1]')
        print(data01[:,z0+1])
    
        print('len(data01[1:data_len,z0+1])')
        print(len(data01[1:data_len,z0+1]))
        print('data01[1:data_len,z0+1]')
        print(data01[1:data_len,z0+1])
print(avgs)

for i0 in range(0,z_slices):
    print('Full avgs:{}, {}'.format(data01[0,i0+1],avgs[i0]))
    print('Blk avgs: {}, {}; SEM: {}'.format(data01[0,i0+1],np.mean(blk_means[i0]),conc_blk_SEMs[i0]))
    print('blk_avgs: {}'.format(blk_avgs[i0]))

#print('Z:')
#print(data01[0,1:z_slices+1])
#print('Avg Conc:')
#print(avgs[:])

# Average from full trajectory
#ax.plot(avgs[:],data01[0,1:z_slices+1],alpha=0.85,color='r',linewidth=1.,label='z=0 avg')

# Average from blocks
ax.plot(blk_avgs[:],data01[0,1:z_slices+1],alpha=0.85,color='r',linewidth=1.,label='z=0 avg')
ax.fill_betweenx(data01[0,1:z_slices+1],blk_avgs[:]-2*conc_blk_SEMs[:],blk_avgs[:]+2*conc_blk_SEMs[:],alpha=0.33,color='r',linewidth=0.,label='z=0 avg')


################################################################################
# Title
################################################################################
#ax.set_title('$O_2$ Concentration Probability Distribution at Z=30$\AA$',fontname='Times New Roman',fontweight='bold',fontsize=16)
ax.set_title('SEM vs Block size for 5 $O_2$ D1CSVR',fontname='Times New Roman',fontweight='bold',fontsize=14)

################################################################################
# X Axis Styling
################################################################################

ax.set_xlabel('Block size', fontname='Times New Roman', fontweight='bold', fontsize=12)

ax.tick_params(axis='x',labelsize=12)
for tick in ax.get_xticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

plt.xscale('log')
ax.set_xlim([0.00025,0.042])

################################################################################
# Y Axis Styling
################################################################################

ax.set_ylabel('SEM', fontname='Times New Roman', fontweight='bold', fontsize=14)

plt.yticks(np.arange(-100,101,10))
ax.tick_params(axis='y',labelsize=10)
for tick in ax.get_yticklabels():
    tick.set_fontname('Times New Roman')
    tick.set_fontweight('bold')

#plt.yscale('log')
ax.set_ylim([-102.,102.])
  
#ax.legend(loc='center',fontsize=8)

################################################################################
# Grid
################################################################################
plt.grid(1,which='both',dashes=(1,1),linewidth=0.5)

fig.set_size_inches(4,7)
#fig.savefig('61-Conc_v_Z_20200508_dimCheck.png',dpi=500)
fig.savefig('SEM_v_blk_size.png',dpi=500,transparent=True,bbox_inches='tight')

