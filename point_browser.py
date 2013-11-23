import matplotlib.pyplot as plt
import jump_interface
import numpy as np
import libtfr as tfr
from IPython import embed
import brewer2mpl as brew

class PointBrowser:
	def __init__(self, jdb):
		self.jdb = jdb
		self.jump_index = 0
		
		#load data from db
		self.jumps = jdb.get_jumps()
		self.clusters = jdb.get_clusters()
		self.signal_indices = jdb.get_signal_indices()
		
		#Only retain jumps of quality 1
		
		quality = jdb.get_jump_indices_quality(1)
		self.jumps = self.jumps[quality]
		self.clusters = self.clusters[quality]
		self.signal_indices = self.signal_indices[quality]
		


		self.__setup_figures()

	def __setup_figures(self):
		num_clusters = np.amax(self.clusters)+1

		#create figures
		self.fig_jumps = plt.figure()
		self.fig_mgram = plt.figure()
		
		#create axes
		self.ax_jumps = self.fig_jumps.add_subplot(111, aspect= 'equal')
		self.ax_mgram = self.fig_mgram.add_subplot(111) 
		
		#set axes properties
		self.ax_jumps.set_title('After Jump Frequency Vs Before Jump \
			Frequency')
		self.ax_jumps.set_xlabel('F1')
		self.ax_jumps.set_ylabel('F2')
		
		#plot self.jumps
		line,=self.ax_jumps.plot(self.jumps[:,0],self.jumps[:,1],'o', \
			markersize=2, visible = False , picker=5)
	
		color_map = brew.get_map('Set2', 'qualitative',\
		 num_clusters).mpl_colors
		for n in xrange(0, num_clusters):
			#plot the data set a second time and color code them by cluster
			idx = np.where(self.clusters==n)[0]
			line, = self.ax_jumps.plot(self.jumps[idx][:,0], \
				self.jumps[idx][:,1],'o',markersize=2.3, \
				color = color_map[n])
			#line, = self.ax_jumps.plot(self.jumps[idx][:,0],self.jumps[idx][:,1],'o',markersize=2, color = plt.get_cmap('Paired')((float(n)+1)/(num_clusters-1)))
		
		self.selected,= self.ax_jumps.plot([self.jumps[0,0]], \
		[self.jumps[0,1]], 'o', ms=12, alpha=0.4, color='red', \
		visible=False)

		#Connect click and press events to callback functions
		self.fig_jumps.canvas.mpl_connect('key_press_event', self.on_press)
		self.fig_mgram.canvas.mpl_connect('key_press_event', self.on_press)
		self.fig_jumps.canvas.mpl_connect('pick_event', self.on_pick)
		#embed()

		plt.show()
		
								  
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
		self.selected.set_data(self.jumps[self.jump_index,0], self.jumps[self.jump_index,1])
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
		distances = np.hypot(x-self.jumps[event.ind,0], y-self.jumps[event.ind,1])
		m = distances.argmin()
		self.jump_index = int(event.ind[m])	
		if self.jump_index is None: return

		#draw yellow circle
		self.selected.set_data(self.jumps[self.jump_index,0], self.jumps[self.jump_index,1])
		self.selected.set_visible(True)
		self.fig_jumps.canvas.draw()
		
		#update signal data and graph mgram
		self.get_signal_data()
		self.graph_mgram()

	def get_signal_data(self):
		#self.signal_index = self.jdb.get_signal_index(self.jump_index) 
		signal_index = self.signal_indices[self.jump_index]
		self.signal = self.jdb.get_signal(signal_index) 
		self.t, self.s = np.split(self.jdb.get_curve(signal_index),[1],1)
		self.jump = self.jumps[self.jump_index]
		self.jump.shape = (1,4)
		print(self.jump)
		
	def graph_mgram(self):
		self.ax_mgram.cla()
	
		#define constants
		Fs = 300e3
		nf = 512
		tstep = 30
		yskip = 20
		tskip = 100
		nf2=nf/2+1 #number of display points
		yconv=int(Fs/(2*nf2))

		#calculate mgram
		mgram = tfr.tfr_spec(self.signal, nf, tstep, nf)
		mgram[mgram==0] = np.amin(mgram[mgram>0])

		#setup labels
		self.ax_mgram.set_title('Reassinged Spectrogram With Curve Fit')
		self.ax_mgram.set_xlabel('Time (s)')
		self.ax_mgram.set_ylabel('Frequency (Hz)')

		#setup time axis coordinates 
		t_axis_coords = np.arange(0, mgram.shape[1] , tskip)
		t_data_coords = np.arange(0, mgram.shape[1]*tstep/Fs, tstep*tskip/Fs)
		t_data_coords = np.array([round(elem, 2) for elem in t_data_coords])
		self.ax_mgram.set_xticks(t_axis_coords)
		self.ax_mgram.set_xticklabels(t_data_coords)
		
		#setup frequency axis coordintes
		f_axis_coords = np.arange(0, nf2, yskip)
		f_data_coords = np.arange(0, Fs/2, yskip*yconv)
		self.ax_mgram.set_yticks(f_axis_coords)
		self.ax_mgram.set_yticklabels(f_data_coords)

		#plot mgram
		self.ax_mgram.imshow(np.log(mgram), vmin=np.mean(np.log(mgram)), origin='lower')
		
		#plot jumps and fitted curve in axis (reassinged) coordinates
		#the data coordinates are for display only
		self.ax_mgram.plot(self.t*Fs/tstep, self.s/yconv, color='green', label='Fitted Curve')
		self.ax_mgram.plot(self.jump[:,2]*Fs/tstep, self.jump[:,0]/yconv, 'mo', label='Jump Points')
		self.ax_mgram.plot(self.jump[:,3]*Fs/tstep, self.jump[:,1]/yconv,'mo')
	
		self.fig_mgram.canvas.draw()

if __name__ == '__main__':
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	rat = 'V6'
	jdb = jump_interface.JumpInterface(dbPath, rat)
	PointBrowser(jdb)