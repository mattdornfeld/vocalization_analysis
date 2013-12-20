import jump_interface
import cluster_parameter
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons, LassoSelector, RadioButtons 
from matplotlib.path import Path
import numpy as np
import libtfr as tfr
from IPython import embed
import brewer2mpl as brew
from operator import not_
from collections import OrderedDict


"""
Class for investigating jumps. It plots the after vs before jump freqs.You 
can click on them to display the reassinged spectrogram in another window.
"""

COLOR_MAP = (brew.get_map('Set3', 'qualitative', 10).mpl_colors +
	brew.get_map('Set2', 'qualitative', 6).mpl_colors)

LEGEND = OrderedDict([(-1,("unclustered", COLOR_MAP[1])), (0,("2 to 1", COLOR_MAP[0])), 
(1,("3 to 2", COLOR_MAP[2])), (2,("4 to 3", COLOR_MAP[3])), 
(3,("5 to 4", COLOR_MAP[4])), (4,("6 to 5", COLOR_MAP[5])), 
(5,("5 to 6", COLOR_MAP[6])), (6,("4 to 5", COLOR_MAP[7])), 
(7,("3 to 4", COLOR_MAP[8])), (8,("2 to 3", COLOR_MAP[9])), 
(9,("1 to 2", COLOR_MAP[10]))])

REVERSE_LEGEND = OrderedDict([("unclustered", -1), ("2 to 1", 0), ("3 to 2", 1), ("4 to 3", 2), 
("5 to 4", 3), ("6 to 5", 4), ("5 to 6", 5), ("4 to 5", 6), ("3 to 4", 7), ("2 to 3", 8),
("1 to 2", 9)])

CLUSTERS = OrderedDict([(0, {'n1':2, 'n2':1}), (1, {'n1':3, 'n2':2}), 
	(2, {'n1':4, 'n2':3}), (3, {'n1':5, 'n2':4}), (4, {'n1':6, 'n2':5}), 
	(5, {'n1':5, 'n2':6}), (6, {'n1':4, 'n2':5}), (7, {'n1':3, 'n2':4}), 
	(8, {'n1':2, 'n2':3}), (9, {'n1':1, 'n2':2})])

NUM_CLUSTERS = len(LEGEND) - 1

class PointBrowser:
	def __init__(self, rat, ji):
		self.ji = ji
		self.rats = ji.get_rats()
		if rat not in self.rats:
			raise Exception("That rat doesn't exist.")

		self.clusters = OrderedDict()
		self.included_clusters = OrderedDict()
		self.jumps = OrderedDict()
		self.signal_indices = OrderedDict()
		self.slopes = OrderedDict()
		for r in self.rats:
			l = len(ji.get_jumps(r))
			self.clusters[r] = -1 * np.ones(l)
			#by default don't include 1 to 2 or 2 to 1
			self.included_clusters[r] = range(1, NUM_CLUSTERS-1) 
			self.jumps[r] = ji.get_jumps(r)
			self.signal_indices[r] = ji.jget_signal_indices(r)
			self.slopes[r] = []

		self.rat_index = np.where(self.rats == rat)[0][0]
		self.rat = rat
		self.setup_figures()
		plt.show()

		
	def change_rat(self, event):
		if event.key not in ('up', 'down'): return
		if event.key == 'up':
			self.rat_index = (self.rat_index + 1)%len(self.rats)
		else: 
			self.rat_index = (self.rat_index - 1)%len(self.rats)
		
		self.rat = self.rats[self.rat_index]
		self.setup_figures()

	
	def reload_defaults(self, event):
		if event.key != "ctrl+" + 'R': return
		self.jumps[self.rat] = ji.get_jumps(self.rat)
		self.signal_indices[self.rat] = ji.jget_signal_indices(self.rat)
		self.clusters[self.rat] = -1 * np.ones(len(self.jumps[self.rat])) 

		self.setup_figures()

	def setup_figures(self):
		self.jump_index = 0
		print(self.rat)
		if 'fig_jumps' not in self.__dict__.keys():
			#create figures
			self.fig_jumps = plt.figure()
			self.fig_mgram = plt.figure()
			self.fig_hist = plt.figure()
			self.fig_hist_phi = plt.figure()
			
			#Connect click and press events to callback functions
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.reload_defaults)
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.on_press)
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.change_rat)
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.cluster_data)
			self.fig_mgram.canvas.mpl_connect('key_press_event', self.on_press)
			self.picker = self.fig_jumps.canvas.mpl_connect('pick_event', self.on_pick)
			
			#create axes
			self.ax_jumps = self.fig_jumps.add_subplot(111, aspect='equal')
			self.ax_mgram = self.fig_mgram.add_subplot(111)
			self.ax_hist = self.fig_hist.add_subplot(111)
			self.ax_hist_phi = self.fig_hist_phi.add_subplot(111)

			#create cluster checker
			rect = [0.00, 0.7, 0.1, 0.25] #l,b,w,h
			self.ax_checker = self.fig_jumps.add_axes(rect)

			#create lasso/point selector toggles
			rect = [0.00, 0.4, 0.1, 0.25] #l,b,w,h
			self.ax_toggle = self.fig_jumps.add_axes(rect)
			labels = ["Point Selector", "Exclusive Lasso", "Inclusive Lasso"]
			self.toggle = RadioButtons(self.ax_toggle, labels, 0)
			self.toggle.on_clicked(self.selection_mode)


		else:
			self.ax_jumps.cla()
			self.ax_mgram.cla()
			self.ax_hist.cla()
			self.ax_checker.cla()

		#update cluster checker
		labels = tuple(LEGEND[n][0] for n in range(NUM_CLUSTERS))
		actives = tuple(n in self.included_clusters[self.rat] 
			for n in range(NUM_CLUSTERS))
		self.check = CheckButtons(self.ax_checker, labels, actives)
		self.check.on_clicked(self.on_check) 

		#set axes properties
		title = (str(self.rat) + " After Jump Frequency Vs Before Jump"  
			" Frequency")
		self.ax_jumps.set_title(title)
		self.ax_jumps.set_xlabel('F1')
		self.ax_jumps.set_ylabel('F2')

		#plot jumps invisibly for click selector
		line,=self.ax_jumps.plot(self.jumps[self.rat][:,0], 
			self.jumps[self.rat][:,1],'o', 
			markersize=2, visible = False , picker=5)

		f = np.arange(20e3, 80e3, 100)
		
		#plot unclustered jumps first
		n = -1
		idx = np.where(self.clusters[self.rat]==n)[0]
		line, = self.ax_jumps.plot(self.jumps[self.rat][idx][:,0], 
				self.jumps[self.rat][idx][:,1],'o',markersize=2.3, 
				color = LEGEND[n][1], label = LEGEND[n][0])
		#plot clustered jumps
		for c in self.included_clusters[self.rat]:
			idx = np.where(self.clusters[self.rat]==c)[0]
			print(str(c)+':'+str(len(idx)))
			print(LEGEND[c])
			line, = self.ax_jumps.plot(self.jumps[self.rat][idx][:,0], 
				self.jumps[self.rat][idx][:,1], 'o', markersize=2.3, 
				color = LEGEND[c][1], label = LEGEND[c][0])
			if len(self.slopes[self.rat]) != 0:
				line, = self.ax_jumps.plot(f, f*self.slopes[self.rat][c], 
					color = LEGEND[c][1])

		#plot histogram
		hist, bins = np.histogram(self.jumps[self.rat][:,1]/self.jumps[self.rat][:,0], bins=50)
		width = 0.7 * (bins[1] - bins[0])
		center = (bins[:-1] + bins[1:]) / 2
		self.ax_hist.bar(center, hist, align='center', width=width)

		"""
		for c in self.included_clusters[self.rat]:
			n1 = CLUSTERS[c]['n1']
			n2 = CLUSTERS[c]['n2']
			idx = np.where(clutsers)
			phi = ((n1*self.jumps[self.rat][:,1] - n2*self.jumps[self.rat][:,0]) 
				/ (self.jumps[self.rat][:,1]-self.jumps[self.rat][:,0]))
			hist, bins = np.histogram(phi, bins=50)
			width = 0.7 * (bins[1] - bins[0])
			center = (bins[:-1] + bins[1:]) / 2
			self.ax_hist_phi.bar(center, hist, align='center', width=width)
		"""



		#create legend self.ax_jumps.legend()
		handles, labels = self.ax_jumps.get_legend_handles_labels()
		rect = [0.1, 0.8, 0.1, 0.15] #l,b,w,h
		self.ax_legend = self.fig_jumps.add_axes(rect)
		self.ax_legend.get_xaxis().set_visible(False)
		self.ax_legend.get_yaxis().set_visible(False)
		self.ax_legend.set_axis_bgcolor((0.75, 0.75, 0.75, 1.0))
		self.ax_legend.set_frame_on(False)
		self.ax_legend.legend(handles, labels, markerscale=3)

		#plot red dot
		self.selected,= self.ax_jumps.plot([self.jumps[self.rat][0,0]], 
			[self.jumps[self.rat][0,1]], 'o', color='red', visible=False)
		
		#set axis limits
		self.ax_jumps.set_xlim(f[0], f[-1])
		self.ax_jumps.set_ylim(f[0], f[-1])
		
		#draw stuff
		self.fig_jumps.canvas.draw()
		self.fig_mgram.canvas.draw()
		self.fig_hist.canvas.draw()

	def selection_mode(self, label):
		labels = ["Point Selector", "Exclusive Lasso", "Inclusive Lasso"]
		if label == labels[0]:
			if 'lasso' in self.__dict__.keys():
				self.lasso.disconnect_events()
			self.fig_jumps.canvas.mpl_connect('pick_event', self.on_pick)
		elif label == labels[1]:
			self.lasso = LassoSelector(self.ax_jumps, self.exclusive_on_lasso)
			self.fig_jumps.canvas.mpl_disconnect(self.picker)
		elif label == labels[2]:
			self.lasso = LassoSelector(self.ax_jumps, self.inclusive_on_lasso)
			self.fig_jumps.canvas.mpl_disconnect(self.picker) 

	def exclusive_on_lasso(self, verts):
		p = Path(verts)
		ind = p.contains_points(self.jumps[self.rat][:,0:2])
		ind = map(not_, ind)
		ind = np.where(ind)[0]
		self.jumps[self.rat] = self.jumps[self.rat][ind,:]
		self.signal_indices[self.rat] = self.signal_indices[self.rat][ind]
		self.clusters[self.rat] = self.clusters[self.rat][ind]
		self.setup_figures()

	def inclusive_on_lasso(self, verts):
		p = Path(verts)
		ind = p.contains_points(self.jumps[self.rat][:,0:2])
		ind = np.where(ind)[0]
		self.jumps[self.rat] = self.jumps[self.rat][ind,:]
		self.signal_indices[self.rat] = self.signal_indices[self.rat][ind]
		self.clusters[self.rat] = self.clusters[self.rat][ind]
		self.setup_figures()

	def on_check(self, label):
		cluster = REVERSE_LEGEND[label]
		if cluster in self.included_clusters[self.rat]:
			self.included_clusters[self.rat].remove(cluster)
		else:
			self.included_clusters[self.rat].insert(cluster, cluster)

		self.included_clusters[self.rat].sort()	
		print(self.included_clusters[self.rat])

		self.fig_jumps.canvas.draw()
	
	def cluster_data(self, event):
		if event.key != "ctrl+" + 'C': return
		self.clusters[self.rat], self.slopes[self.rat] = cluster_parameter.main(
			self.jumps[self.rat], self.included_clusters[self.rat])
		self.setup_figures()
		for n, m in zip(self.included_clusters[self.rat], self.slopes[self.rat]):
			label = LEGEND[n][0]
			top = self.ax_hist.get_ylim()[1]
			self.ax_hist.annotate(label, (m, top))

	
	def on_press(self, event):
		#Show next plot if your press, previous plot if your press p
		if self.jump_index is None: return
		if event.key not in ('n', 'p'): return
		if event.key=='n': inc = 1
		else:  inc = -1
		
		#increment jump index  
		self.jump_index += inc
		self.jump_index = int(np.clip(self.jump_index, 0, len(self.jumps[self.rat])-1))
		
		#draw yellow circle on graph
		self.selected.set_data(self.jumps[self.rat][self.jump_index,0], 
			self.jumps[self.rat][self.jump_index,1])
		self.selected.set_visible(True)
		self.fig_jumps.canvas.draw()

		#update signal data and graph mgram
		self.get_signal_data()
		self.graph_mgram()
	
	def on_pick(self, event):		
		#event.ind is the list of all jump points within a certain
		#distance of selected point 
		N = len(event.ind)
		if not N: return True
		
		#find closest jump point to selected location
		x = event.mouseevent.xdata
		y = event.mouseevent.ydata
		distances = np.hypot(x-self.jumps[self.rat][event.ind,0], 
			y-self.jumps[self.rat][event.ind,1])
		m = distances.argmin()
		self.jump_index = int(event.ind[m])	
		if self.jump_index is None: 

			return

		#draw yellow circle
		self.selected.set_data(self.jumps[self.rat][self.jump_index,0], 
		 self.jumps[self.rat][self.jump_index,1])
		self.selected.set_visible(True)
		self.fig_jumps.canvas.draw()
		
		#update signal data and graph mgram
		self.get_signal_data()
		self.graph_mgram()

	def get_signal_data(self):
		#self.signal_index = self.ji.get_signal_index(self.jump_index) 
		signal_index = self.signal_indices[self.rat][self.jump_index]
		self.signal = self.ji.get_signal(signal_index) 
		self.t, self.s = np.split(self.ji.get_curve(signal_index),[1],1)
		self.jump = self.jumps[self.rat][self.jump_index]
		self.fs = self.ji.get_fs(signal_index)
		print(self.jump)
		
	def graph_mgram(self):
		self.ax_mgram.cla()
	
		#define constants
		nf = 512
		tstep = 30
		yskip = 20
		tskip = 100
		nf2=nf/2+1 #number of display points
		yconv=int(self.fs/(2*nf2))

		#calculate mgram
		mgram = tfr.tfr_spec(self.signal, nf, tstep, nf)
		mgram[mgram==0] = np.amin(mgram[mgram>0])

		#setup labels
		self.ax_mgram.set_title('Reassinged Spectrogram With Curve Fit')
		self.ax_mgram.set_xlabel('Time (s)')
		self.ax_mgram.set_ylabel('Frequency (Hz)')

		#setup time axis coordinates 
		t_axis_coords = np.arange(0, mgram.shape[1] , tskip)
		t_data_coords = np.arange(0, mgram.shape[1]*tstep/self.fs, 
			tstep*tskip/self.fs)
		t_data_coords = np.array([round(elem, 2) for elem in t_data_coords])
		self.ax_mgram.set_xticks(t_axis_coords)
		self.ax_mgram.set_xticklabels(t_data_coords)

		#setup frequency axis coordintes
		f_axis_coords = np.arange(0, nf2, yskip)
		f_data_coords = np.arange(0, self.fs/2, yskip*yconv)
		self.ax_mgram.set_yticks(f_axis_coords)
		self.ax_mgram.set_yticklabels(f_data_coords)

		#plot mgram
		self.ax_mgram.imshow(np.log(mgram), vmin=np.mean(np.log(mgram)), 
			origin='lower')

		"""
		plot jumps and fitted curve in axis (reassinged) coordinates
		the data coordinates are for display only
		"""
		self.ax_mgram.plot(self.t*self.fs/tstep, self.s/yconv, color='green', 
			label='Fitted Curve')
		self.ax_mgram.plot(self.jump[2]*self.fs/tstep, self.jump[0]/yconv, 
			'mo', label='Jump Points')
		self.ax_mgram.plot(self.jump[3]*self.fs/tstep, self.jump[1]/yconv,'mo')

		self.fig_mgram.canvas.draw()


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

if __name__ == '__main__':
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	rat = 'V6'
	ji = jump_interface.JumpInterface(dbPath)
	PointBrowser(rat, ji)