from matplotlib.widgets import Slider
import numpy as np
import argparse
from customTools import jumpDatabaseAPI
import matplotlib.pyplot as plt
from matplotlib import mlab
from math import atan, tan
from matplotlib.path import Path
import matplotlib.patches as patches
from IPython import embed

def costFunction(jumps, cluster, slopes):
	cost = 0
	for i,c in enumerate(np.unique(cluster[cluster!=-1])): 
		idx = mlab.find(cluster==c)
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
def polyCluster(jumps, bisectingSlopes):
	cluster = -1*np.ones(len(jumps))
	codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
	vertices = polygonVertices(bisectingSlopes)
	nClusters = len(bisectingSlopes)+1
	for iv, v in enumerate(vertices):
		poly = Path(v, codes, closed = True)
		for ij, j in enumerate(jumps):
			if poly.contains_point(j)==1:
				cluster[ij] = iv
	
	return cluster
	
def plotJumps(jumps, cluster, clusterSlopes):
	nClusters = len(clusterSlopes)
	ax1 = fig1.add_subplot(111, aspect='equal')
	ax1.set_title('After Jump Frequency Vs Before Jump Frequency')
	ax1.set_xlabel('F1')
	ax1.set_ylabel('F2')
	for n in xrange(0,nClusters):
		#plot the data set a second time and color code them by cluster
		idx = mlab.find(cluster==n)
		minf = np.amin(jumps[idx,:][:,0])
		maxf = np.amax(jumps[idx,:][:,0])
		f = np.arange(minf, maxf, 100)
		line, = ax1.plot(jumps[idx,:][:,0],jumps[idx,:][:,1], 'o', markersize=2, color = plt.get_cmap('jet')((float(n)+1)/(nClusters)))
		line, = ax1.plot(f, clusterSlopes[n]*f, color = plt.get_cmap('jet')((float(n)+1)/(nClusters)))
	
	plt.draw()

def makePlot(r):
	ax1.cla()
	slopes = [(n+1-r)/(n-r) for n in range(2,4)+[5]]
	slopes[len(slopes):]=([(n-r)/(n+1-r) for n in [5]+range(3,1,-1)])
	slopes = np.array(slopes)
	slopes = np.sort(slopes)
	bisectingSlopes= findBisector(slopes)
	cluster = polyCluster(jumps, bisectingSlopes)
	print('\r'+'r='+str(r))
	costFunction(jumps, cluster, slopes)
	plotJumps(jumps, cluster, slopes)
	plt.draw()


if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--dbPath", help="Data path")
	parser.add_argument("-r", "--rat", help="Rat Name")
	args = parser.parse_args()
	
	jdb = jumpDatabaseAPI.jumpDatabaseAPI(args.dbPath, args.rat)
	jumps = jdb.getJumps()
	valid = np.load('/home/matthew/Dropbox/Work/Python/pointBrowsers/indices.npz')
	valid = valid['yes'][:,1].astype(int)
	jumps = jumps[valid,:]
	#embed()

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, aspect='equal')
	plt.subplots_adjust(bottom=0.25)
	makePlot(r=0.25)
	axr = fig1.add_axes([0.25, 0.05, 0.65, 0.03], axisbg='lightgoldenrodyellow')
	r = Slider(axr, 'r', -1, 1, valinit=0.25, dragging = True)
	r.on_changed(makePlot)
	plt.show()		
