import argparse
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

def slope(n1, n2, r):
	return float(n2-r) / (n1-r)

SLOPES = OrderedDict([(0, partial(slope,2,1)), (1, partial(slope,3,2)), (2, partial(slope,4,3)), (3, partial(slope,5,4)),
(4,partial(slope,6,5)), (5, partial(slope,5,6)), (6,partial(slope,4,5)), 
(7,partial(slope,3,4)), (8,partial(slope,2,3)), (9, partial(slope,1,2))])

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

#linear regression forced through origin
def regression(x,y):
	x = np.mat(x)
	y = np.mat(y)
	x.shape=(1,x.size)
	y.shape=(1,y.size)
	m= y * x.T * np.linalg.inv(x*x.T)
	return float(m)

#Error of all the clusters. We're trying to minimize this.
def cost_function(jumps, cluster, slopes):
	cost = 0 
	for c in slopes.keys(): 
		idx = find(cluster==c)
		l = len(jumps[idx,0:2])
		for j in jumps[idx,0:2]:
			x = (j[0] / slopes[c] + j[1]) / (slopes[c] + 1 / slopes[c])
			y = slopes[c] * x
			cost += (((j[0] - x)**2 + (j[1]-y)**2)) / l 
	
	return cost

#takes an array of slopes and returns the slopes of the bisectors of these lines sorted along with the original slopes	
def find_bisector(slopes):
	keys = slopes.keys()
	bisecting_slopes = np.zeros(len(slopes)-1)
	for i, c in enumerate(keys):
		if c == keys[-1]: break
		m1 = slopes[c] 
		m2 = slopes[keys[i+1]]
		bisecting_slopes[i] = tan(atan(m1)/2 + atan(m2)/2)
	
	return bisecting_slopes

#takes slopes of lines, min, and max x values. returns vertices of quadrigon that lie on lines at min and max values. 
def polygon_vertices(slopes, included_clusters, minf = 20e3, maxf = 85e3):
	#vertices = [0 for x in range(len(slopes)+1)]
	vertices = [(n,[]) for n in included_clusters]
	vertices = OrderedDict(vertices)
	
	for i, m1 in enumerate(slopes):
		c = included_clusters[i] 
		if i == 0:
			vertices[c] = np.array([(minf, 0), (minf, m1*minf), (maxf, m1*maxf),
			 (maxf, 0), (minf, 0)])
		else:
			m2 = slopes[i-1]
			vertices[c] = np.array([(minf, minf*m2), (maxf, m2*maxf), 
				(maxf, m1*maxf), (minf, m1*minf), (minf, minf*m2)]			)
		
	m1 = slopes[-1]
	vertices[included_clusters[-1]] = np.array([(minf, m1*minf),
	 (minf, m1*maxf), (maxf, m1*maxf), (minf, m1*minf), (minf, m1*minf)])
			
	return vertices

#creates polygons based on the lines in slopes and checks to see which jumps are in which polygon
def poly_cluster(jumps, vertices, codes = 
	[Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]):

	cluster = -1*np.ones(len(jumps))
	num_clusters = len(vertices)
	for c, v in vertices.iteritems():
		poly = Path(v, codes, closed = True)
		for ij, j in enumerate(jumps):
			if poly.contains_point(j)==1:
				cluster[ij] = c
	
	return cluster

#shuffle the points around between neighboring clusters to minimize error
def shuffle(jumps, cluster, slopes):
	num_clusters = len(slopes)
	keys = slopes.keys()
	for ij, j in enumerate(jumps):
		#embed()
		c = int(cluster[ij])
		if c==-1: continue
		idx = keys.index(c) 
		e1 = j[1]-slopes[c]*j[0] #original error
		#handle edge cases first 
		if (c==keys[0] and e1<0) or (c==keys[-1] and e1>0): continue
		elif c==keys[0] and e1>0:
			nc = keys[idx+1] #next cluster
			e2 = j[1]-slopes[nc]*j[0]
			if abs(e2)<abs(e1): cluster[ij]=c+1
		elif c==keys[-1] and e1<0:
			pc = keys[idx-1] #prev cluster
			e2 = j[1]-slopes[pc]*j[0]
			if abs(e2)<abs(e1): cluster[ij]=c-1
		else:
			#if above the line and new error is less than original bring it up one cluster
			#if below the line bring it down one
			if e1>0:
				nc = keys[idx+1] #next cluster
				e2 = j[1]-slopes[nc]*j[0]
				if abs(e2)<abs(e1): cluster[ij]=c+1
			if e1<0:
				pc = keys[idx-1] #prev cluster
				e2 = j[1]-slopes[pc]*j[0]
				if abs(e2)<abs(e1): cluster[ij]=c-1
	
	return cluster
		
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

def insert_jumps(ji, jumps, signal_indices, cluster):
	#put the bad jumps on top
	q1 = np.ones((len(jumps)))
	q0 = ji.get_jump_indices_quality(0)
	bad_jumps = ji.get_jumps()[q0]
	bad_signal_indices = ji.get_signal_indices()[q0]
	bad_cluster = -1 * np.ones(len(bad_jumps))
	jumps = np.vstack((jumps, bad_jumps))
	#embed()
	signal_indices = np.concatenate((signal_indices, bad_signal_indices))
	cluster = np.concatenate((cluster, bad_cluster))
	quality = np.concatenate((q1,q0))
	ji.insert_jumps(jumps, signal_indices, cluster, quality)


def main(jumps, included_clusters = range(NUM_CLUSTERS)):
	#Organize input from db and command line
	#dbPath = args.fileName
	#rat = args.rat
	#jumps = ji.get_jumps(rat)
	#signal_indices = ji.jget_signal_indices(rat)
	#jump_indices = ji.get_jump_indices(rat)
	
	#Brute force through different values of r(theta)
	#Find error of each on
	rationals = np.arange(-1,1,0.01)
	rTrue = 0
	e = []
	for r in rationals:
		print(r)
		slopes = [(n,SLOPES[n](r)) for n in included_clusters]
		slopes = OrderedDict(slopes)
		bisecting_slopes= find_bisector(slopes)
		vertices = polygon_vertices(bisecting_slopes, included_clusters)
		cluster = poly_cluster(jumps[:,0:2], vertices)
		cluster = shuffle(jumps[:,0:2], cluster, slopes)
		e.append(cost_function(jumps[:,0:2], cluster, slopes))
	
	#take r with min eror
	e = np.array(e)
	r_min = rationals[np.where(e==min(e))][0]
	print('r_min=' + str(r_min)) 
	slopes = [(n,SLOPES[n](r_min)) for n in included_clusters]
	slopes = OrderedDict(slopes)
	bisecting_slopes= find_bisector(slopes)
	vertices = polygon_vertices(bisecting_slopes, included_clusters)
	cluster = poly_cluster(jumps[:,0:2], vertices)
	cluster = shuffle(jumps[:,0:2], cluster, slopes)

	plot_jumps(jumps[:,0:2], cluster, slopes)
	
	#plot_polygons(vertices, slopes)
	
	fig_error = plt.figure()
	ax_error = fig_error.add_subplot(111)
	ax_error.plot(rationals, e)
	ax_error.set_xlabel('theta')
	ax_error.set_ylabel('error')
	ax_error.set_title("r = " + str(r_min))
	fig_error.show()
	
	#print(slopes)

	return cluster, slopes
	"""
	for jump_index in jump_indices:
		ji.update_cluster(jump_index, cluster)
	#point_browser.PointBrowser(ji)
	"""
	

if __name__ == "__main__":
	#Define command line arguments accepted
	#parser = argparse.ArgumentParser()
	#parser.add_argument("-f", "--fileName", help="Data path")
	#parser.add_argument("-r", "--rat", help="Rat Name")
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
		jumps = ji.get_jumps(rat)
		main(jumps, included_clusters)

