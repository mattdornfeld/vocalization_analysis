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

#linear regression forced through origin
def regression(x,y):
	x = np.mat(x)
	y = np.mat(y)
	x.shape=(1,x.size)
	y.shape=(1,y.size)
	m= y * x.T * np.linalg.inv(x*x.T)
	return float(m)

#Error of all the clusters. We're trying to minimize this.
def costFunction(jumps, cluster, slopes):
	cost = 0
	for i,c in enumerate(np.unique(cluster[cluster!=-1])): 
		idx = find(cluster==c)
		l = len(jumps[idx,0:2])
		for j in jumps[idx,0:2]:
			x = (j[0] / slopes[i] + j[1]) / (slopes[i] + 1 / slopes[i])
			y = slopes[i] * x
			cost += (((j[0] - x)**2 + (j[1]-y)**2)) / l 
	
	return cost


def costFunction2(jumps, cluster, slopes):
	cost = 0
	for i,c in enumerate(np.unique(cluster[cluster!=-1])): 
		idx = find(cluster==c)
		l = len(jumps[idx,0:2])
		cost=0
		for j in jumps[idx,0:2]:
			cost += abs(j[0]*slopes[i] - j[1])/l	

		print(cost)
	
	return cost

#takes an array of slopes and returns the slopes of the bisectors of these lines sorted along with the original slopes	
def findBisector(slopes):
	bisectingSlopes = np.zeros(len(slopes)-1)
	for i, m1 in enumerate(slopes):
		if i == len(slopes)-1:
			break
		m2 = slopes[i+1]
		bisectingSlopes[i] = tan(atan(m1)/2+atan(m2)/2)
	
	return bisectingSlopes

#takes slopes of lines, min, and max x values. returns vertices of quadrigon that lie on lines at min and max values. 
def polygonVertices(slopes, minf = 20e3, maxf = 85e3):
	vertices = [0 for x in range(len(slopes)+1)]	
	
	for i, m1 in enumerate(slopes):
		if i == 0:
			vertices[i] = np.array([(minf, 0), (minf, m1*minf), (maxf, m1*maxf), (maxf, 0), (minf, 0)])
		else:
			m2 = slopes[i-1]
			vertices[i] = np.array([(minf, minf*m2), (maxf, m2*maxf), (maxf, m1*maxf), (minf, m1*minf), (minf, minf*m2)]			)
		
	m1 = slopes[-1]
	vertices[-1] = np.array([(minf, m1*minf), (minf, m1*maxf), (maxf, m1*maxf), (minf, m1*minf), (minf, m1*minf)])
			
	return vertices

#creates polygons based on the lines in slopes and checks to see which jumps are in which polygon
def polyCluster(jumps, vertices, codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]):
	cluster = np.zeros(len(jumps))
	num_clusters = len(vertices)
	for iv, v in enumerate(vertices):
		poly = Path(v, codes, closed = True)
		for ij, j in enumerate(jumps):
			if poly.contains_point(j)==1:
				cluster[ij] = iv + 1
	
	return cluster

#shuffle the points around between neighboring clusters to minimize error
def shuffle(jumps, cluster, slopes):
	num_clusters = int(np.amax(cluster)+1)
	for ij, j in enumerate(jumps):
		c = cluster[ij]
		if c==-1: continue
		e1 = j[1]-slopes[c]*j[0] #original error
		#handle edge cases first 
		if (c==0 and e1<0) or (c==num_clusters-1 and e1>0): 
			continue
		elif c==0 and e1>0:
			e2 = j[1]-slopes[c+1]*j[0]
			if abs(e2)<abs(e1): cluster[ij]=c+1
		elif c==num_clusters-1 and e1<0:
			e2 = j[1]-slopes[c-1]*j[0]
			if abs(e2)<abs(e1): cluster[ij]=c-1
		else:
			#if above the line and new error is less than original bring it up one cluster
			#if below the line bring it down one
			if e1>0:
				e2 = j[1]-slopes[c+1]*j[0]
				if abs(e2)<abs(e1): cluster[ij]=c+1
			if e1<0:
				e2 = j[1]-slopes[c-1]*j[0]
				if abs(e2)<abs(e1): cluster[ij]=c-1
	
	return cluster
		
def plotJumps(jumps, cluster, slopes):
	num_clusters = len(slopes)
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, aspect='equal')
	ax1.set_title('After Jump Frequency Vs Before Jump Frequency')
	ax1.set_xlabel('F1')
	ax1.set_ylabel('F2')
	transitions = [str(n+1)+' to ' +str(n) for n in range(2,6)]
	transitions = transitions+[str(n)+' to ' +str(n+1) for n in range(5,1,-1)]


	for n in xrange(0,num_clusters):
		#plot the data set a second time and color code them by cluster
		idx = find(cluster==n)
		if len(idx)==0: continue
		minf = np.amin(jumps[idx,:][:,0])
		maxf = np.amax(jumps[idx,:][:,0])
		f = np.arange(minf, maxf, 100)
		line, = ax1.plot(jumps[idx,:][:,0],jumps[idx,:][:,1], 'o', markersize=2, color = plt.get_cmap('jet')((float(n)+1)/(num_clusters)))
		line, = ax1.plot(f, slopes[n]*f, color = plt.get_cmap('jet')((float(n)+1)/(num_clusters)), linewidth=1.3, label=transitions[n])
	ax1.legend()

def plotPolygons(vertices, slopes):
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111) #, aspect='equal')
	ax2.set_title('Illustration of Polygon Clustering Technique')
	ax2.set_xlabel('F1')
	ax2.set_ylabel('F2')
	x = np.arange(20e3, 85e3, 100)
	num_clusters=len(slopes)
	codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
	
	#for m in slopes:
	#	ax2.plot(x, m*x)
	
	for n,v in enumerate(vertices):
		ax2.plot(v[:,0], v[:,1], 'ro')
		poly = Path(v, codes, closed = True)
		print(plt.get_cmap('jet')((float(n)+1)/(num_clusters-1)))
		patch = patches.PathPatch(poly, facecolor = plt.get_cmap('jet')((float(n)+1)/(num_clusters)), lw=2)
		ax2.add_patch(patch)	

def findSlopes(r):
	slopes = [(n+1-r)/(n-r) for n in range(2,6)]
	slopes[len(slopes):]=([(n-r)/(n+1-r) for n in range(5,1,-1)])
	slopes = np.array(slopes)
	slopes = np.sort(slopes)

	return slopes

#sort jumps first based on cluster then based on distance from the origin	
def sort_jumps(jumps, signal_indices, cluster):
	sorted_jumps = np.zeros(jumps.shape)
	sorted_signal_indices = np.zeros(len(signal_indices))
	distance = np.sqrt(jumps[:,0]**2+jumps[:,1]**2)
	idx = cluster.argsort() 
	cluster = cluster[idx] #cluster is now sorted
	jumps = jumps[idx,:]
	distance = distance[idx]
	signal_indices = signal_indices[idx]  
	for c in np.unique(cluster):
		idx1 = find(cluster == c)
		j = jumps[idx1,:]
		idx2 = distance[idx1].argsort()
		sorted_jumps[idx1,:] = jumps[idx1,:][idx2,:]
		sorted_signal_indices[idx1] = signal_indices[idx1][idx2]
	
	return sorted_jumps, sorted_signal_indices, cluster

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


def main(ji, rat):
	#Organize input from db and command line
	#dbPath = args.fileName
	#rat = args.rat
	jumps = ji.get_jumps(rat)
	signal_indices = ji.jget_signal_indices(rat)
	jump_indices = ji.get_jump_indices(rat)
	
	#Brute force through different values of r(theta)
	#Find error of each one
	rationals = np.arange(0,1,0.01)
	rTrue = 0
	e = []
	for r in rationals:
		print(r)
		slopes = findSlopes(r)
		bisectingSlopes= findBisector(slopes)
		vertices = polygonVertices(bisectingSlopes)
		cluster = polyCluster(jumps[:,0:2], vertices)
		cluster = shuffle(jumps[:,0:2], cluster, slopes)
		e.append(costFunction(jumps[:,0:2], cluster, slopes))
	
	#take r with min eror
	e = np.array(e)
	rmin = rationals[np.where(e==min(e))][0]
	print('rmin=' + str(rmin)) 
	slopes = findSlopes(rmin)
	bisectingSlopes= findBisector(slopes)
	vertices = polygonVertices(bisectingSlopes)
	cluster = polyCluster(jumps[:,0:2], vertices)

	#plotJumps(jumps[:,0:2], cluster, slopes)
	#plotPolygons(vertices, slopes)
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	ax2.plot(rationals, e)
	ax2.set_xlabel('theta')
	ax2.set_ylabel('error')
	fig2.show()
	#plt.show()

	return cluster
	"""
	for jump_index in jump_indices:
		ji.update_cluster(jump_index, cluster)
	#point_browser.PointBrowser(ji)
	"""
	

if __name__ == "__main__":
	#Define command line arguments accepted
	parser = argparse.ArgumentParser()
	#parser.add_argument("-f", "--fileName", help="Data path")
	parser.add_argument("-r", "--rat", help="Rat Name")
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	ji = jump_interface.JumpInterface(dbPath)
	args = parser.parse_args()
	main(ji, args.rat)
