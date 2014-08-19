import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
import jump_interface
import minimize_cost
import os
from IPython import embed
from collections import OrderedDict
from math import sqrt
import brewer2mpl as brew
from matplotlib.ticker import AutoMinorLocator

COLOR_MAP = (brew.get_map('Set3', 'qualitative', 10).mpl_colors +
	brew.get_map('Accent', 'qualitative', 6).mpl_colors)

DB_PATH = "/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db"
RAT_CLUSTERS = {'V1':[3,6], 'V2':[3,6], 'V3':[3,6], 'V4':[1,2,3,4,5,6,7,8],
        'V5':[3,6], 'V6':[1,2,3,4,5,6,7,8], 'V15':[1,3,6,8], 'V17':[3,6], 
        'V18':[3,6], 'V31':[1,3,6,8], 'V32': [3,6]}

def get_file_list():
	root_path = '/home/matthew/Dropbox/Work/vocalization_analysis/resampled_theta'
	file_list = []
	for root, _, filenames in os.walk(root_path):
		for filename in filenames:
			file_list.append(os.path.join(root, filename))

	return file_list

def set_rc_parameters():
	fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
	inches_per_pt = 1.0/72.27               # Convert pt to inch
	golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	fig_width = fig_width_pt*inches_per_pt  # width in inches
	fig_height = fig_width*golden_mean      # height in inches
	fig_size =  [fig_width,fig_height]
	params = {'savefig.format': 'png',
	         'axes.labelsize': 10,
	         'text.fontsize': 10,
	         'legend.fontsize': 10,
	         'xtick.labelsize': 8,
	         'ytick.labelsize': 8,
	         'text.usetex': True,
	         'figure.figsize': fig_size,
	         'figure.dpi': 100}
	plt.rcParams.update(params)

def confidence_interval(measurements, confidence=0.025):
	num_bins = int(np.sqrt(len(measurements)))
	#plt.hist(measurements, bins=num_bins)
	density, edges = np.histogram(measurements, bins=num_bins)
	density = density / float(sum(density))
	distribution = np.cumsum(density)

	a = sum(distribution < confidence)
	a = edges[a]
	b = len(edges) - sum(distribution > (1-confidence))
	b = edges[b]

	return [a,b]

def plot_results():
	ji = jump_interface.JumpInterface(DB_PATH)
	rats = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V15', 'V17', 'V18', 
	'V31', 'V32']
	labels = ['rat 1', 'rat 2', 'rat 3', 'rat 4', 'rat 5', 'rat 6', 'rat 7', 'rat 8', 'rat 9', 
	'rat 10', 'rat 11', 'edge tone', 'shallow cavity', 'deep cavity', 'hole tone']
	whistles = [('edge tone',-0.54), ('shallow cavity',-0.17), ('deep cavity',0.08), ('hole tone',0.25)]
	whistles = OrderedDict(whistles)
	initials = [[.61,0],[0.38,0],[.36,0],[0.2,0],[0.28,0],[0.43,0],[0.5,0],[.25,0],[0.4,0],[0.39,0],[0.64,0]]
	#initials = [[-.25,0] for _ in range(12)]

	fig_t, ax_t = plt.subplots()
	ax_t.set_ylim([0,len(labels)+1])
	ax_t.set_yticks(range(1,len(labels)+1))
	ax_t.set_yticklabels(labels)
	ax_t.set_xlabel('$\gamma$')
	ax_t.set_xlim([-0.6,0.8])
	#ax_t.set_ylabel('rat')
	ax_t.grid(True)
	minorLocator   = AutoMinorLocator()


	"""
	fig_b, ax_b = plt.subplots()
	ax_b.set_ylim([0,12])
	ax_b.set_yticks(range(1,12))
	ax_b.set_yticklabels(rats)
	ax_b.set_xlabel('b')
	ax_b.set_ylabel('rat')
	fig_b.tight_layout()
	"""

	i = 1
	for rat, initial in zip(rats, initials):
		print(rat)
		f = '/home/matthew/Dropbox/Work/vocalization_analysis/resampled_theta/'+ rat + '_measurements.npy'
		data = np.load(f)
		jumps = ji.get_jumps(rat)
		included_clusters = RAT_CLUSTERS[rat]
		theta = data[:,0]
		b = data[:,1]
		theta_min, b_min = minimize_cost.main(jumps, initial, included_clusters)
		ax_t.plot(1*theta_min, i, color='blue', marker='o' )
		#ax_b.plot(b_min, i, 'bo')
		theta_error = confidence_interval(theta)
		b_error = confidence_interval(b)
		ax_t.plot(theta_error, [i,i], color='blue', marker='|')
		#ax_b.plot(b_error, [i,i], color='blue', marker='|')
		i+=1

	for w, g in whistles.iteritems():
		ax_t.plot(g, i, color='red', marker='^')
		i+=1

	ax_t.xaxis.set_minor_locator(minorLocator)
	fig_t.tight_layout()
	fig_t.savefig('/home/matthew/work/writing/jump_paper/theta_error.png')
	plt.show()

if __name__ == '__main__':
	set_rc_parameters()
	plot_results()

