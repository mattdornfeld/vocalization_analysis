import jump_interface
import cluster_parameter
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import numpy as np
import libtfr as tfr
from IPython import embed
import brewer2mpl as brew

"""
Class for investigating jumps. It plots the after vs before jump freqs.You 
can click on them to display the reassinged spectrogram in another window.
"""
COLOR_MAP = brew.get_map('Set3', 'qualitative', 10).mpl_colors
LEGEND = {-1:("unclustered", COLOR_MAP[1]), 0:("3 to 2", COLOR_MAP[0]), 
1:("4 to 3", COLOR_MAP[2]), 2:("5 to 4", COLOR_MAP[3]), 
3:("6 to 5", COLOR_MAP[4]), 4:("5 to 6", COLOR_MAP[5]), 
5:("4 to 5", COLOR_MAP[6]), 6:("3 to 4", COLOR_MAP[7]), 
7:("2 to 3", COLOR_MAP[8])}

class PointBrowser:
	def __init__(self, rat, ji):
		self.ji = ji
		self.rats = ji.get_rats()
		if rat not in self.rats:
			raise Exception("That rat doesn't exist.")

		self.clusters = dict()
		self.slopes = dict()
		for r in self.rats:
			l = len(ji.get_jumps(r))
			self.clusters[r] = -1 * np.ones(l) 
			self.slopes[r] = []

		self.rat_index = np.where(self.rats == rat)[0][0]
		self.rat = rat
		self.initialize()
		plt.show()

		
	def change_rat(self, event):
		if event.key not in ('up', 'down'): return
		if event.key == 'up':
			self.rat_index = (self.rat_index + 1)%len(self.rats)
		else: 
			self.rat_index = (self.rat_index - 1)%len(self.rats)
		
		self.rat = self.rats[self.rat_index]
		self.initialize()

	
	def initialize(self):
		self.jump_index = 0
		print(self.rat)

		#load data from db
		self.jumps = ji.get_jumps(self.rat)
		self.signal_indices = ji.jget_signal_indices(self.rat)

		self.setup_figures()

	def setup_figures(self):
		if 'fig_jumps' not in self.__dict__.keys():
			#create figures
			self.fig_jumps = plt.figure()
			self.fig_mgram = plt.figure()
			
			#Connect click and press events to callback functions
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.on_press)
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.change_rat)
			self.fig_jumps.canvas.mpl_connect('key_press_event', self.cluster_data)
			self.fig_mgram.canvas.mpl_connect('key_press_event', self.on_press)
			self.fig_jumps.canvas.mpl_connect('pick_event', self.on_pick)
			
			#create axes
			self.ax_jumps = self.fig_jumps.add_subplot(111)
			self.ax_mgram = self.fig_mgram.add_subplot(111)

			#create cluster checker
			rect = [0.00, 0.8, 0.1, 0.15]
			button_axis = self.fig_jumps.add_axes(rect)
			labels = tuple(LEGEND[n][0] for n in range(8))
			actives = tuple(False for n in range(8))
			check = CheckButtons(button_axis, labels, actives)
			check.on_clicked(self.on_check)  


		else:
			self.ax_jumps.cla()
			self.ax_mgram.cla()

		#set axes properties
		title = (str(self.rat) + " After Jump Frequency Vs Before Jump"  
			" Frequency")
		self.ax_jumps.set_title(title)
		self.ax_jumps.set_xlabel('F1')
		self.ax_jumps.set_ylabel('F2')

		#plot jumps invisibly for click selector
		line,=self.ax_jumps.plot(self.jumps[:,0],self.jumps[:,1],'o', 
			markersize=2, visible = False , picker=5)

		f = np.arange(20e3, 80e3, 100)
		
		#plot unclustered juumps first
		n = -1
		idx = np.where(self.clusters[self.rat]==n)[0]
		line, = self.ax_jumps.plot(self.jumps[idx][:,0], 
				self.jumps[idx][:,1],'o',markersize=2.3, 
				color = LEGEND[n][1], label = LEGEND[n][0])
		#plot clustered jumps
		for n in range(8):
			idx = np.where(self.clusters[self.rat]==n)[0]
			print(str(n)+':'+str(len(idx)))
			line, = self.ax_jumps.plot(self.jumps[idx][:,0], 
				self.jumps[idx][:,1],'o',markersize=2.3, 
				color = LEGEND[n][1], label = LEGEND[n][0])
			if len(self.slopes[self.rat]) != 0:
				line, = self.ax_jumps.plot(f, f*self.slopes[self.rat][n], 
					color = LEGEND[n][1])

		self.ax_jumps.legend()
		
		self.selected,= self.ax_jumps.plot([self.jumps[0,0]], 
			[self.jumps[0,1]], 'o', color='red', visible=False)
		self.ax_jumps.set_xlim(f[0], f[-1])
		self.ax_jumps.set_ylim(f[0], f[-1])
		self.fig_jumps.canvas.draw()
		self.fig_mgram.canvas.draw()

		#embed()
	
	def on_check(self, label):
		print(label)
		print(1)
		self.fig_jumps.canvas.draw()
	
	def cluster_data(self, event):
		if event.key != "ctrl+" + 'C': return
		self.clusters[self.rat], self.slopes[self.rat] = cluster_parameter.main(self.ji, self.rat)
		self.setup_figures()

	
	def on_press(self, event):
		#Show next plot if your press, previous plot if your press p
		if self.jump_index is None: return
		if event.key not in ('n', 'p'): return
		if event.key=='n': inc = 1
		else:  inc = -1
		
		#increment jump index  
		self.jump_index += inc
		self.jump_index = int(np.clip(self.jump_index, 0, len(self.jumps)-1))
		
		#draw yellow circle on graph
		self.selected.set_data(self.jumps[self.jump_index,0], 
			self.jumps[self.jump_index,1])
		self.selected.set_visible(True)
		self.fig_jumps.canvas.draw()

		#update signal data and graph mgram
		self.get_signal_data()
		self.graph_mgram()
	
	def on_pick(self, event):	
		print(event)	
		#event.ind is the list of all jump points within a certain
		#distance of selected point 
		N = len(event.ind)
		if not N: return True
		
		#find closest jump point to selected location
		x = event.mouseevent.xdata
		y = event.mouseevent.ydata
		distances = np.hypot(x-self.jumps[event.ind,0], 
			y-self.jumps[event.ind,1])
		m = distances.argmin()
		self.jump_index = int(event.ind[m])	
		if self.jump_index is None: 

			return

		#draw yellow circle
		self.selected.set_data(self.jumps[self.jump_index,0], 
		 self.jumps[self.jump_index,1])
		self.selected.set_visible(True)
		self.fig_jumps.canvas.draw()
		
		#update signal data and graph mgram
		self.get_signal_data()
		self.graph_mgram()

	def get_signal_data(self):
		#self.signal_index = self.ji.get_signal_index(self.jump_index) 
		signal_index = self.signal_indices[self.jump_index]
		self.signal = self.ji.get_signal(signal_index) 
		self.t, self.s = np.split(self.ji.get_curve(signal_index),[1],1)
		self.jump = self.jumps[self.jump_index]
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

if __name__ == '__main__':
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	rat = 'V6'
	ji = jump_interface.JumpInterface(dbPath)
	PointBrowser(rat, ji)