import libtfr as tfr
import argparse
import numpy as np
from matplotlib.mlab import find
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import label
from scipy.stats import linregress
from customtools import arrayStack
import jump_interface
from IPython import embed

global Fs, tstep, nf, nf2
Fs=250e3 #sampling frequency 
tstep=30	#number of time steps in original coordinates in one time step of reassigned coordinate
nf=512 # number of frequency points in spectrogram
nf2=nf/2+1 #number of display points
yconv=int(Fs/(2*nf2))

def mgram(signal):
	mt = tfr.tfr_spec(signal,nf,tstep,nf)
	mt[mt==0] = np.amin(mt[mt>0])
	return mt
	
#returns low entropy time points and the time averaged entropy
def lowH(mt,thresh):
	pt=sum(mt,0)
	pn=mt/pt
	valid=-sum(pn*np.log2(pn),0) #entropy as a function of time
	averageH=sum(valid)/len(valid) #time average of H
	valid[valid<thresh]=1
	valid[valid>thresh]=0

	return valid, averageH

#yskip=number of frequency steps in reassigned coodinates to skip between ticks when plotting mgrams
#tskip=number of time steps in reassigned coordinates to skip between ticks when plotting mgrams
def makePlot(s, t, signal, jumps, yskip = 20, tskip = 100):
	signal = np.reshape(signal,signal.size)
	mt = mgram(signal)
	plt.imshow(np.log(mt), vmin=np.mean(np.log(mt)), origin='lower')
	fdisp = np.arange(0,Fs/2,yskip*yconv) 
	plt.yticks(np.arange(0,nf2,yskip),fdisp)
	tdisp = np.arange(0,(mt.shape[1]*tstep)/Fs,tstep*tskip/Fs)
	tdisp = np.array([round(elem, 2) for elem in tdisp])
	plt.xticks(np.arange(0,mt.shape[1],tskip),tdisp)
	plt.plot(t*Fs/tstep,s/yconv, color='red',label='Fitted Curve')
	
	if jumps!=[]:
		jumps.shape=(jumps.size/4,4)		
		plt.plot(jumps[:,2]*Fs/tstep,jumps[:,0]/yconv,'go',label='Jump Points')
		plt.plot(jumps[:,3]*Fs/tstep,jumps[:,1]/yconv,'go')
	
	plt.xlabel('Time (s)')
	plt.ylabel('Frequency (Hz)')
	plt.title('Spectrogram With Curve Fit')

def fitCurve(mt, thresh = 6.9):
	realFreqs = np.arange(int(10e3/yconv), int(100e3/yconv)) #indices of realistic frequency range 10e3 to 100e3 hz
	tr = np.shape(mt)[1]
	s=np.zeros(tr)
	ts = np.arange(0, np.shape(mt)[1]*tstep,tstep)/Fs
	valid, averageH=lowH(mt,thresh) #valid are time points with entropy below thresh. This function also returns the time average of H.
	#print valid		
	for t in range(0,tr):
		if valid[t]==0:
			continue
		
		#m = np.argmax(mt[realFreqs, t], axis=0)
		m = find(mt[:,t]==max(mt[realFreqs,t]))[0]
		
		try:
			pt=sum(mt[m-3:m+4,t])
		except:
			continue
		
		s[t]=sum(mt[m-3:m+4, t]*np.arange(m-3,m+4)*yconv/pt)
	
	return s, ts, valid, averageH
	
#input curves and valid points returns jump times and frequencies 
def jumpFrequencies(s,t,valid):
	split_ds = 5e3
	split_dt = 1e-3
	ds = np.diff(s)/split_ds
	dt = np.diff(t)/split_dt
	dr = np.sqrt(np.power(ds,2) + np.power(dt,2)) #compute distance of jump from origin
	window = np.ones(t.size-1)	#Discard points within 10% of edge
	window[0:0.1*t.size] = 0
	window[0.9*t.size:t.size] = 0
	dr = window*dr
	l = label(dr<1)	#jump points have dr<1	
	index = find(l[0]==0)	#indices of jump points	
	remove = np.sort(np.append(find(np.diff(index)<5), find(np.diff(index)<5)+1))	#remove points whose indices are within 5 places of each other
	index = np.delete(index, remove)
	
	remove = np.array([])    
	for i in range(0,index.size):
		if dr[index[i]]<2:
			remove = np.append(remove,i)
	
	remove.dtype = int
	index = np.delete(index,remove)
	jumps = np.zeros((100,4))
	
	for i in range(0,len(index)):
		#This if statement checks several conditions. If any of them are satisified the jump point is rejected. 1) Is jump point in high entropy region, 2) Is the before or the after jump frequency to high or too low, 
		#3) Is the the difference between before and after frequencies too high or too low
		if find(valid[index[i]-15:index[i]+15]==0)!=[] or s[index[i]]<15e3 or s[index[i]+1]<15e3 or s[index[i]]>80e3 or s[index[i]+1]>80e3 or abs(s[index[i]]-s[index[i]+1])<7e3 or abs(s[index[i]+1]-s[index[i]])>25e3:
			continue
		else:
			#not sure what this does
			if i != len(index)-1 and abs(index[i+1]-index[i])<10:				
				continue
			if i!=0 and abs(index[i-1]-index[i])<10:				
				continue
		jumps[i,:] = np.array([s[index[i]],s[index[i]+1],t[index[i]],t[index[i]+1]])
	
	jumps = np.delete(jumps,find(jumps==0))
	jumps.shape = (len(jumps)/4, 4)
	
	return jumps
	
#does linear regressions on four points before jump and four points after jump. 
#It calculates the new jump frequencies as the values of the regession lines at the midpoint of the jump times
def refineJump(jump, s, ts):
	i1, i2 = jump[2]*Fs/tstep, jump[3]*Fs/tstep #before after jump time indices
	m1, b1 = linregress(ts[i1-3:i1+1], s[i1-3:i1+1])[0:2]  
	m2, b2 = linregress(ts[i2:i2+4], s[i2:i2+4])[0:2]	
	tmid = (jump[3]+jump[2])/2
	f1 = b1+m1*tmid # refined before jump freq 
	f2 = b2+m2*tmid # refined after jump freq
	
	jump[0] = f1
	jump[1] = f2
	
	return jump

def main(dbPath):
	ji = jump_interface.JumpInterface(dbPath)
	#signalIndices = ji.getSignalIndices()
	#signalIndices = np.unique(signalIndices)
	
	rats = ji.get_rats()
	for rat in rats:
		jumps = arrayStack.arrayStack(5)
		rat = str(rat)
		ji.change_rat(rat)
		signal_indices = ji.get_signal_indices_rat()
		for signal_index in signal_indices:
			signal = ji.get_signal(signal_index)
			mt = mgram(signal)
			s, t, valid, averageH = fitCurve(mt)
			signalJumps = jumpFrequencies(s, t, valid)
			curve = np.append(s,t)
			curve = curve.reshape((len(s),2))
			ji.insert_curve(curve) 
			
			for j, jump in enumerate(signalJumps):
				jump = refineJump(jump, s, t)
				if jump != []:
					jumps.push(np.append(jump, signal_index))
		
		jumps = jumps.returnArray()
		cluster = -1 * np.ones(len(jumps))
		quality = 1 * np.ones(len(jumps))
		ji.insert_jumps(jumps[:,0:4], jumps[:,4], cluster, quality)  
		
		"""	
		fig = plt.figure()
		ax = fig.add_subplot(111)
		print(jumps.arrayShape())
		print(signalJumps.shape)
		print(jumps.returnArray())
		makePlot(s, t, signal, signalJumps)
		plt.show()
		"""
		
		

		
		

if __name__ == "__main__":
	#Define command line arguments accepted
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--dbPath", help="Data path")
	#parser.add_argument("-r", "--rat", help="ratid")
	args = parser.parse_args()
	
	main(args.dbPath)
