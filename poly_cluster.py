import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from math import atan, tan
from collections import OrderedDict

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

def main(jumps, slopes, included_clusters):
	bisecting_slopes = find_bisector(slopes)
	vertices = polygon_vertices(bisecting_slopes, included_clusters)
	cluster = poly_cluster(jumps, vertices)

	return cluster

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