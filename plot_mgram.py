import numpy as np
import matplotlib.pyplot as pyplot
import libtfr as tfr
import jump_interface

def set_rc_parameters():
	fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
	inches_per_pt = 1.0/72.27               # Convert pt to inch
	golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	fig_width = fig_width_pt*inches_per_pt  # width in inches
	fig_height = fig_width*golden_mean      # height in inches
	fig_size =  [fig_width,fig_height]
	params = {'savefig.format': 'png',
	         'axes.labelsize': 10,
	         'text.fontsize': 10,
	         'legend.fontsize': 10,
	         'xtick.labelsize': 8,
	         'ytick.labelsize': 'medium',
	         'ytick.major.size': 2,
	         'text.usetex': True,
	         'figure.figsize': fig_size,
	         'figure.dpi': 100}
	plt.rcParams.update(params)

def get_signal_data(jump_index):
	rat = 'V6'
	#signal_index = ji.get_signal_index(jump_index) 
	signal_index = signal_indices[rat][jump_index]
	signal = ji.get_signal(signal_index) 
	t, s = np.split(ji.get_curve(signal_index),[1],1)
	jump = jumps[rat][jump_index]
	fs = ji.get_fs(signal_index)
	
	return signal, s, t, jump, fs

def low_h(mt, thresh):
	pt = sum(mt,0)
	pn = mt / pt
	valid = -sum(pn*np.log2(pn),0) #entropy as a function of time
	average_h = sum(valid)/len(valid) #time average of H
	valid[valid<thresh] = 1
	valid[valid>thresh] = 0

	return valid, average_h

def graph_mgram(signal, s, t, jump, fs):
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
	fig_mgram, ax_mgram = plt.subplots()
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

	#plot mgram
	im = ax_mgram.imshow(np.log(mgram), vmin=np.mean(np.log(mgram)), 
		origin='lower')
	im.set_cmap('gray')

	"""
	plot jumps and fitted curve in axis (reassinged) coordinates
	the data coordinates are for display only
	"""
	ax_mgram.plot(t*fs/tstep, s/yconv, color='orange', 
		label='Fitted Curve', linewidth='3')
	ax_mgram.plot(jump[2]*fs/tstep, jump[0]/yconv, 
		label='Jump Points', marker='o', color='#6599FF', markersize=8)
	ax_mgram.plot(jump[3]*fs/tstep, jump[1]/yconv, 
		marker='o', color='#6599FF', markersize=8)

	ax_mgram.set_ylim([0,80e3/yconv])
	fig_mgram.tight_layout()

if __name__ == '__main__':
	set_rc_parameters()
	jump_index = 100
	signal, s, t, jump, fs = get_signal_data(jump_index)
	graph_mgram(signal, s, t, jump, fs)

