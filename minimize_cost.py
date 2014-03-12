import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import find
import sqlite3 as lite
import jump_interface
from IPython import embed
from functools import partial
import brewer2mpl as brew
from collections import OrderedDict
from scipy.optimize import fmin_ncg
import poly_cluster
from matplotlib import animation

def slope(n1, n2, theta):
	return float(n2-theta) / (n1-theta)

def dslope(n1, n2, theta):
	return float(n1-n2) / (n1 - theta)**2

SLOPES = OrderedDict([(0, partial(slope,2,1)), (1, partial(slope,3,2)), 
	(2, partial(slope,4,3)), (3, partial(slope,5,4)), (4,partial(slope,6,5)), 
	(5, partial(slope,5,6)), (6,partial(slope,4,5)), (7,partial(slope,3,4)), 
	(8,partial(slope,2,3)), (9, partial(slope,1,2))])

DSLOPES = OrderedDict([(0, partial(dslope,2,1)), (1, partial(dslope,3,2)), 
	(2, partial(dslope,4,3)), (3, partial(dslope,5,4)), (4,partial(dslope,6,5)), 
	(5, partial(dslope,5,6)), (6,partial(dslope,4,5)), (7,partial(dslope,3,4)), 
	(8,partial(dslope,2,3)), (9, partial(dslope,1,2))])

COLOR_MAP = (brew.get_map('Set3', 'qualitative', 10).mpl_colors +
	brew.get_map('Accent', 'qualitative', 6).mpl_colors)

LEGEND = OrderedDict([(-1,("unclustered", COLOR_MAP[2])), (0,("2 to 1", COLOR_MAP[0])), 
(1,("3 to 2", COLOR_MAP[15])), (2,("4 to 3", COLOR_MAP[3])), 
(3,("5 to 4", COLOR_MAP[4])), (4,("6 to 5", COLOR_MAP[5])), 
(5,("5 to 6", COLOR_MAP[6])), (6,("4 to 5", COLOR_MAP[7])), 
(7,("3 to 4", COLOR_MAP[8])), (8,("2 to 3", (131./255,39./255,39./255))), 
(9,("1 to 2", COLOR_MAP[10]))])

REVERSE_LEGEND = OrderedDict([("unclustered", -1), ("2 to 1", 0), ("3 to 2", 1), ("4 to 3", 2), 
("5 to 4", 3), ("6 to 5", 4), ("5 to 6", 5), ("4 to 5", 6), ("3 to 4", 7), ("2 to 3", 8),
("1 to 2", 9)])

NUM_CLUSTERS = len(LEGEND) - 1

#Error of all the clusters. We're trying to minimize this.
def dcost(x, *args):
	theta = x[0]
	b = x[1]
	jumps = args[0]
	included_clusters = args[1]

	h = [(n,SLOPES[n](theta)) for n in included_clusters]
	h = OrderedDict(h)
	dh = [(n, DSLOPES[n](theta)) for n in included_clusters]
	dh = OrderedDict(dh)
	
	cluster = poly_cluster.main(jumps, h, included_clusters)

	dcost_theta = 0
	dcost_b = 0
	for c in h.keys(): 
		idx = find(cluster==c)
		l = len(jumps[idx,:])
		for j in jumps[idx,:]:
			x = (j[0] / h[c] + j[1] - b) / (h[c] + 1 / h[c])
			y = h[c] * x + b 
			
			dx_theta = ((h[c]**2 + 1) * (j[1] - b) * dh[c] - 
				2 * (j[0] + (j[1] - b) * h[c]) * h[c] * dh[c]) / (h[c]**2 + 1)**2
			
			dy_theta = ((h[c]**2 + 1) * (j[0] * dh[c] + 2 * (j[1] - b) * h[c] * dh[c]) 
				- 2 * (j[0] * h[c] + (j[1] - b) * h[c]**2) * h[c] * dh[c]) / (h[c]**2 + 1)**2
			dx_b = -h[c] / (h[c]**2 + 1)
			dy_b = -h[c]**2 / (h[c]**2 + 1)

			dcost_theta += (2 * (j[0] - x) * dx_theta 
				+ 2 * ((j[1] - b) - y) * dy_theta) / l

			dcost_b += (2 * (j[0] - x) * dx_b 
				+ 2 * ((j[1] - b) - y) * dy_b) / l

	return np.array([dcost_theta, dcost_b])

def cost(x, *args):
	theta = x[0]
	b = x[1]
	jumps = args[0]
	included_clusters = args[1]
	h = [(n,SLOPES[n](theta)) for n in included_clusters]
	h = OrderedDict(h)

	cluster = poly_cluster.main(jumps, h, included_clusters)

	cost = 0 
	for c in h.keys(): 
		idx = find(cluster==c)
		l = len(jumps[idx,:])
		for j in jumps[idx,:]:
			x = (j[0] / h[c] + j[1] - b) / (h[c] + 1 / h[c])
			y = h[c] * x + b 
			cost += (((j[0] - x)**2 + (j[1]-y)**2)) / l 

	return cost

def cost_animator(n, *args):
	line = args[0]
	theta = args[1]
	error = args[2]
	line.set_data(theta[0:n], error[0:n])
	line.set_marker('o')
	line.set_linestyle('')
	line.set_color('red')

	return line,


def update_cluster_animator(n, line, jumps, cluster_frame, h_frame, 
	included_clusters):
	cluster = cluster_frame[n]
	h = h_frame[n]
	f1 = np.arange(20e3,80e3,100)
	for i, c  in enumerate(included_clusters):
		idx = find(cluster==c)
		line[i].set_data(jumps[idx,0],jumps[idx,1])
		line[i+l/2].set_data(f1, h[c]*f1)


	return line

def animate_clustering(jumps, included_clusters, trajectory):
	"""
	fig = plt.figure()
	ax = plt.axes(xlim=(0, 1), ylim=(0, 4e7))
	line, = ax.plot([], [], lw=2)
	"""

	"""
	theta = np.arange(0,1,0.01)
	error = np.array([cost([t,10], jumps, included_clusters) for t in theta])
	ax.plot(theta, error)
	"""

	fig = plt.figure()
	ax = plt.axes(xlim=(20e3, 80e3), ylim=(20e3, 80e3))
	line = []
	l = 2*len(included_clusters)
	for i in range(l):
		line += ax.plot([],[])
	
	for i in range(l/2):
		c = included_clusters[i]
		line[i].set_color(LEGEND[c][1])
		line[i].set_linestyle('')
		line[i].set_marker('o')
		line[i].set_markersize(2)
		line[i+l/2].set_color(LEGEND[c][1])
		line[i+l/2].set_linestyle('-')

	theta = trajectory[:,0]
	h_frame = []
	cluster_frame = []
	for t in theta:
		h = OrderedDict([(n,SLOPES[n](t)) for n in included_clusters])
		h_frame += [h]
		cluster_frame += [poly_cluster.main(jumps, h, included_clusters)]
	#error = np.array([cost([t,10], jumps, included_clusters) for t in theta])
	anim = animation.FuncAnimation(fig, update_cluster_animator, 
		fargs = (line, jumps, cluster_frame, h_frame, included_clusters), 
		init_func=None, frames=len(theta), interval=1000, blit=True)

	
	"""
	anim = animation.FuncAnimation(fig, cost_animator, 
		fargs = (line, theta, error), 
		init_func=None, frames=len(trajectory)+1, interval=400, blit=True)
	"""

	plt.show()

def main(jumps, included_clusters = range(NUM_CLUSTERS)):
	mins, trajectory = fmin_ncg(f=cost, x0=[0.36,0], fprime=dcost,
	 args=(jumps, included_clusters), avextol=1e-5, disp=0, callback=None, 
	 retall=True)

	theta_min = mins[0]
	b_min = mins[1]
	trajectory = np.array(trajectory)
	#embed()
	#animate_clustering(jumps, included_clusters, trajectory)
	print(theta_min)
	print(b_min)
	
	return [theta_min, b_min]

if __name__ == "__main__":
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	ji = jump_interface.JumpInterface(dbPath)
	"""
	rat_clusters={'V1':[3,6], 'V2':[3,6], 'V3':[3,6], 'V4':[1,2,3,4,5,6,7,8],
		'V5':[3,6], 'V6':[1,2,3,4,5,6,7,8], 'V15':[1,3,6,8], 'V17':[3,6], 
		'V18':[3,6], 'V31':[1,3,6,8], 'V32': [3,6]}
	"""
	rat_clusters={'V6':[1,2,3,4,5,6,7,8]}
	for rat, included_clusters in rat_clusters.iteritems():
		jumps = ji.get_jumps(rat)[:,0:2]
		main(jumps, included_clusters)

