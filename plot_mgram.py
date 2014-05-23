import numpy as np
import matplotlib.pyplot as plt
import libtfr as tfr
import jump_interface
from math import sqrt

def set_rc_parameters():
	fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
	inches_per_pt = 1.0/72.27               # Convert pt to inch
	golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	fig_width = fig_width_pt*inches_per_pt  # width in inches
	fig_height = fig_width*0.4       # height in inches
	fig_size =  [fig_width,fig_height]
	params = {'savefig.format': 'png',
			 'font.size': '6',
	         'axes.labelsize': 10,
	         'text.fontsize': 10,
	         'legend.fontsize': 10,
	         'xtick.labelsize': 8,
	         'ytick.labelsize': 8,
	         'text.usetex': True,
	         'figure.figsize': fig_size,
	         'figure.dpi': 100,
	         'image.aspect': 2}
	plt.rcParams.update(params)

def get_signal_data(ji, jump_index):
	rat = 'V6'
	signal_index = ji.get_signal_index(jump_index) 
	jump_indices = ji.get_jump_indices(signal_index)
	jumps = np.array([ji.get_jump(i) for i in jump_indices])
	signal = ji.get_signal(signal_index) 
	t, s = np.split(ji.get_curve(signal_index),[1],1)
	fs = ji.get_fs(signal_index)
	
	return signal, s, t, jumps, fs

def low_h(mt, thresh):
	pt = sum(mt,0)
	pn = mt / pt
	valid = -sum(pn*np.log2(pn),0) #entropy as a function of time
	average_h = sum(valid)/len(valid) #time average of H
	valid[valid<thresh] = 1
	valid[valid>thresh] = 0

	return valid, average_h

def graph_mgram(signal, s, t, jumps, fs):
	#define constants
	nf = 512
	tstep = 30
	tskip = 250
	nf2=nf/2+1 #number of display points
	yconv = int(fs/(2*nf2))
	fskip = 20e3/yconv

	#calculate mgram
	mgram = tfr.tfr_spec(signal, nf, tstep, nf)
	mgram[mgram==0] = np.amin(mgram[mgram>0])
	valid, _ = low_h(mgram, thresh=6.9)
	valid.shape = (len(valid),1)
	s = s * valid

	#setup labels
	#ax_mgram.set_title('Reassinged Spectrogram With Curve Fit')
	fig_mgram = plt.figure()
	ax_mgram = fig_mgram.add_subplot(111)
	ax_mgram.set_xlabel('Time (s)')
	ax_mgram.set_ylabel('Frequency (Hz)')

	#setup time axis coordinates 
	t_axis_coords = np.arange(0, mgram.shape[1] , tskip)
	t_data_coords = np.arange(0, mgram.shape[1]*tstep/fs, 
		tstep*tskip/fs)
	t_data_coords = np.array([round(elem, 2) for elem in t_data_coords])
	ax_mgram.set_xticks(t_axis_coords)
	ax_mgram.set_xticklabels(t_data_coords)

	#setup frequency axis 
	f_axis_coords = np.arange(0, nf2, fskip)
	f_data_coords = np.arange(0, fs/2, fskip*yconv)
	ax_mgram.set_yticks(f_axis_coords)
	ax_mgram.set_yticklabels(f_data_coords)
	ax_mgram.tick_params(axis='both', which='major', labelsize=20)

	for label in (ax_mgram.get_xticklabels() + ax_mgram.get_yticklabels()):
		label.set_fontsize(13)

	"""
	plot jumps and fitted curve in axis (reassinged) coordinates
	the data coordinates are for display only
	"""
	ax_mgram.plot(t*fs/tstep, s/yconv, color='orange', linewidth='0.5')
	ax_mgram.plot(jumps[:,2]*fs/tstep, jumps[:,0]/yconv, 'o', 
		color='#98E71A', markersize=4)
	ax_mgram.plot(jumps[:,3]*fs/tstep, jumps[:,1]/yconv, 'o', 
		color='#98E71A', markersize=4)

	#plot mgram
	im = ax_mgram.imshow(np.log(mgram), vmin=np.mean(np.log(mgram)), 
		origin='lower', aspect='auto')
	im.set_cmap('Purples')

	ax_mgram.set_ylim([0,80e3/yconv])
	fig_mgram.tight_layout()
	fig_mgram.savefig('/home/matthew/work/writing/jump_paper/specgram.png')

if __name__ == '__main__':
	set_rc_parameters()
	jump_index = 156
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	ji = jump_interface.JumpInterface(dbPath)
	signal, s, t, jumps, fs = get_signal_data(ji,jump_index)
	graph_mgram(signal, s, t, jumps, fs)
	plt.show()

