#!/usr/bin/python

import jump_interface
import cluster_parameter
import numpy as np
import brewer2mpl as brew
from collections import OrderedDict
from pylab import *


COLOR_MAP = (brew.get_map('Set3', 'qualitative', 10).mpl_colors +
	brew.get_map('Set2', 'qualitative', 6).mpl_colors)



NUM_CLUSTERS = len(LEGEND) - 1





dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
ji = jump_interface.JumpInterface(dbPath)
clusters = OrderedDict()
included_clusters = OrderedDict()
jumps = OrderedDict()
signal_indices = OrderedDict()
slopes = OrderedDict()
rats = ['V6','V4']

for r in rats:
	l = len(ji.get_jumps(r))
	clusters[r] = -1 * np.ones(l)
	#by default don't include 1 to 2 or 2 to 1
	included_clusters[r] = range(1, NUM_CLUSTERS-1) 
	jumps[r] = ji.get_jumps(r)
	signal_indices[r] = ji.jget_signal_indices(r)
	slopes[r] = []
	clusters[r], slopes[r] = cluster_parameter.main(
		jumps[r], included_clusters[r])	

fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width      # height in inches
fig_size =  [fig_width,fig_height]
params = {'savefig.format': 'png',
         'axes.labelsize': 10,
         'text.fontsize': 10,
         'legend.fontsize': 10,
         'xtick.labelsize': 8,
         'ytick.labelsize': 8,
         'text.usetex': True,
         'figure.figsize': fig_size,
         'figure.dpi': 100,
         'legend.markerscale':2}
pylab.rcParams.update(params)

LEGEND = OrderedDict([(-1,(r"unclustered", COLOR_MAP[2])), (0,(r"$2\rightarrow1$", COLOR_MAP[0])), 
(1,(r"$3\rightarrow2$", COLOR_MAP[15])), (2,(r"$4\rightarrow3$", COLOR_MAP[3])), 
(3,(r"$5\rightarrow4$", COLOR_MAP[4])), (4,(r"$6\rightarrow5$", COLOR_MAP[5])), 
(5,(r"$5\rightarrow6$", COLOR_MAP[6])), (6,(r"$4\rightarrow5$", COLOR_MAP[7])), 
(7,(r"$3\rightarrow4$", COLOR_MAP[8])), (8,(r"$2\rightarrow3$", (131./255,39./255,39./255))), 
(9,(r"$1\rightarrow2$", COLOR_MAP[10]))])
rat = 'V4'
fig_jumps = plt.figure()
ax_jumps = fig_jumps.add_subplot(111, aspect='equal')
f = np.arange(20e3, 80e3, 100)
ax_jumps.set_xlim(f[0], f[-1])
ax_jumps.set_ylim(f[0], f[-1])
ax_jumps.set_xlabel('f1')
ax_jumps.set_ylabel('f2')
fig_jumps.tight_layout()

#plot clustered jumps
for c in included_clusters[rat]:
	idx = np.where(clusters[rat]==c)[0]
	print(str(c)+':'+str(len(idx)))
	print(LEGEND[c])
	line, = ax_jumps.plot(jumps[rat][idx][:,0], 
		jumps[rat][idx][:,1], 'o', markersize=2.1, 
		color = LEGEND[c][1], label = LEGEND[c][0])
	if len(slopes[rat]) != 0:
		line, = ax_jumps.plot(f, f*slopes[rat][c], 
			color = LEGEND[c][1])

handles, labels = ax_jumps.get_legend_handles_labels()
legend = ax_jumps.legend(handles[::-1], labels[::-1],loc='lower right', numpoints=1, shadow=True)
fig_jumps.canvas.draw()
fig_jumps.savefig('/home/matthew/work/writing/jump_paper/'+rat)

rat = 'V6'

fig_jumps = plt.figure()
ax_jumps = fig_jumps.add_subplot(111, aspect='equal')
f = np.arange(20e3, 80e3, 100)
ax_jumps.set_xlim(f[0], f[-1])
ax_jumps.set_ylim(f[0], f[-1])
ax_jumps.set_xlabel('f1')
ax_jumps.set_ylabel('f2')
fig_jumps.tight_layout()

#plot clustered jumps
for c in included_clusters[rat]:
	idx = np.where(clusters[rat]==c)[0]
	print(str(c)+':'+str(len(idx)))
	print(LEGEND[c])
	line, = ax_jumps.plot(jumps[rat][idx][:,0], 
		jumps[rat][idx][:,1], 'o', markersize=2.1, 
		color = LEGEND[c][1], label = LEGEND[c][0])
	if len(slopes[rat]) != 0:
		line, = ax_jumps.plot(f, f*slopes[rat][c], 
			color = LEGEND[c][1])

		
fig_jumps.canvas.draw()

fig_jumps.savefig('/home/matthew/work/writing/jump_paper/'+rat)
