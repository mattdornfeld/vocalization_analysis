import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import find
from matplotlib.path import Path
import matplotlib.patches as patches
import sqlite3 as lite
import jump_interface
import point_browser
from math import atan, tan, sqrt
from IPython import embed
from functools import partial
import brewer2mpl as brew
from collections import OrderedDict
from scipy.optimize import fmin_ncg
import poly_cluster

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
	
	cluster = poly_cluster.main(h)

	dcost_theta = 0
	dcost_b = 0
	for c in h.keys(): 
		idx = find(cluster==c)
		l = len(jumps)
		for j in jumps:
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

	cluster = poly_cluster.main(h)

	cost = 0 
	for c in h.keys(): 
		idx = find(cluster==c)
		l = len(jumps)
		for j in jumps:
			x = (j[0] / h[c] + j[1] - b) / (h[c] + 1 / h[c])
			y = h[c] * x + b 
			cost += (((j[0] - x)**2 + (j[1]-y)**2)) / l 
	
	#print('b='+str(b))
	#print('theta='+str(theta))
	#print('cost='+str(cost))
	return cost



def main(jumps, included_clusters = range(NUM_CLUSTERS)):
	"""
	fig, ax = plt.subplots()
	ax.plot(jumps[:,0], jumps[:,1], 'o')
	plt.show()
	"""
	

	theta_min, b_min  = fmin_ncg(f=cost, x0=[0,0], fprime=dcost,
	 args=(jumps, included_clusters), avextol=1e-5, disp=0)
	""" 
	slopes = [(n,SLOPES[n](r_min)) for n in included_clusters]
	slopes = OrderedDict(slopes)
	bisecting_slopes= find_bisector(slopes)
	vertices = polygon_vertices(bisecting_slopes, included_clusters)
	cluster = poly_cluster(jumps[:,0:2], vertices)
	cluster = shuffle(jumps[:,0:2], cluster, slopes)
	"""

	#plot_jumps(jumps[:,0:2], cluster, slopes)
	
	#plot_polygons(vertices, slopes)
	
	"""
	fig_error = plt.figure()
	ax_error = fig_error.add_subplot(111)
	ax_error.plot(rationals, e)
	ax_error.set_xlabel('theta')
	ax_error.set_ylabel('error')
	ax_error.set_title("theta = " + str(r_min))
	fig_error.show()
	"""
	
	#print(slopes)

	return theta_min, b_min
	"""
	for jump_index in jump_indices:
		ji.update_cluster(jump_index, cluster)
	#point_browser.PointBrowser(ji)
	"""
	

def plot_jumps(jumps, cluster, slopes):
	num_clusters = len(slopes)
	fig1 = plt.figure(figsize=(12, 12), dpi=100)

	ax1 = fig1.add_subplot(111, aspect='equal')
	#ax1.set_title('After Jump Frequency Vs Before Jump Frequency')
	ax1.set_xlabel('f1', size=12)
	ax1.set_ylabel('f2', size=12)
	#transitions = [str(n+1)+' to ' +str(n) for n in range(2,6)]
	#transitions = transitions+[str(n)+' to ' +str(n+1) for n in range(5,1,-1)]

	f = np.arange(20e3, 80e3, 100)
	for c in slopes.keys():
		#plot the data set a second time and color code them by cluster
		idx = find(cluster==c)
		if len(idx)==0: continue
		minf = np.amin(jumps[idx,:][:,0])
		maxf = np.amax(jumps[idx,:][:,0])
		line, = ax1.plot(jumps[idx,:][:,0],jumps[idx,:][:,1], 'o', 
			markersize=2.3, color = LEGEND[c][1], label=LEGEND[c][0])
		line, = ax1.plot(f, slopes[c]*f, color = LEGEND[c][1])

	#create legend ax_jumps.legend()
	handles, labels = ax1.get_legend_handles_labels()
	rect = [0.73, 0.758, 0.1, 0.15] #l,b,w,h
	ax_legend = fig1.add_axes(rect)
	ax_legend.get_xaxis().set_visible(False)
	ax_legend.get_yaxis().set_visible(False)
	ax_legend.set_axis_bgcolor((0.75, 0.75, 0.75, 1.0))
	ax_legend.set_frame_on(False)
	ax_legend.legend(handles[::-1], labels[::-1],
	 markerscale=3, numpoints=1)

	#set axis limits

	ax1.set_xlim(f[0], f[-1])
	ax1.set_ylim(f[0], f[-1])
	plt.show()
	#fig1.savefig(rat+'.png', bbox_inches='tight', pad_inches=0)

def plot_polygons(vertices, slopes):
	fig_error = plt.figure()
	ax_error = fig_error.add_subplot(111) #, aspect='equal')
	ax_error.set_title('Illustration of Polygon Clustering Technique')
	ax_error.set_xlabel('F1')
	ax_error.set_ylabel('F2')
	x = np.arange(20e3, 85e3, 100)
	num_clusters=len(vertices)
	codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
	
	for m in slopes:
		ax_error.plot(x, m*x)
	
	#embed()
	for n,v in vertices.iteritems():
		print(n)
		ax_error.plot(v[:,0], v[:,1], 'ro')
		poly = Path(v, codes, closed = True)
		print(plt.get_cmap('jet')((float(n)+1)/(num_clusters-1)))
		patch = patches.PathPatch(poly, facecolor = LEGEND[n][1], lw=2)
		ax_error.add_patch(patch)	

if __name__ == "__main__":
	#Define command line arguments accepted
	#parser = argparse.ArgumentParser()
	#parser.add_argument("-f", "--fileName", help="Data path")
	#parser.add_argument("-theta", "--rat", help="Rat Name")
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	ji = jump_interface.JumpInterface(dbPath)
	#args = parser.parse_args()
	"""
	rat_clusters={'V1':[3,6], 'V2':[3,6], 'V3':[3,6], 'V4':[1,2,3,4,5,6,7,8],
		'V5':[3,6], 'V6':[1,2,3,4,5,6,7,8], 'V15':[1,3,6,8], 'V17':[3,6], 
		'V18':[3,6], 'V31':[1,3,6,8], 'V32': [3,6]}
	"""
	rat_clusters={'V6':[1,2,3,4,5,6,7,8]}
	for rat, included_clusters in rat_clusters.iteritems():
		jumps = ji.get_jumps(rat)[:,0:2]
		main(jumps, included_clusters)

