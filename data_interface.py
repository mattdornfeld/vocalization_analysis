import h5py
import numpy as np
import os
from IPython import embed
import jump_interface

class DataInterface:
	def __init__(self):
		pass
	
	def connect(self, path):
		self.f = h5py.File(path)

	def get_Fs(self):
		return self.f['USV_info']['soundFs'][0][0]

	def get_quality(self, mic_num, signal_num):
		ref_USV = self.f['USV'][mic_num][0] 
		group_USV = self.f.__getitem__(ref_USV)
		q_ref = group_USV['q'][signal_num][0]
		q = group_USV.__getitem__(q_ref)[0][0]
		
		return chr(q)

	def get_num_calls(self, mic_num):
		obj = self.f['USVclip'][mic_num][0]
		group = self.f.__getitem__(obj)
		
		return len(group['sound'])

	def get_signal(self, signal_num, mic_num):
		obj = self.f['USVclip'][mic_num][0]
		group = self.f.__getitem__(obj)
		signal = np.array(group.__getitem__(group['sound'][signal_num][0])).flatten() 

		return signal
	
	def get_rat_ids(self):
		group = self.f['USV_info']['stages']
		data_set0 = np.array(group.__getitem__(group['rats'][0][0])).flatten()
		data_set1 = np.array(group.__getitem__(group['rats'][1][0])).flatten()
		rat_ids = [''.join(str(chr(i)) for i in data_set0)]
		rat_ids.append(''.join(str(chr(i)) for i in data_set1))
		
		return rat_ids

def get_file_list():
	root_path = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/Calls/V1-V6'
	file_list = []
	for root, _, filenames in os.walk(root_path):
		for filename in filenames:
			file_list.append(os.path.join(root, filename))

	return file_list

def main():
	db_path = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	ji = jump_interface.JumpInterface(db_path)
	null = chr(0) + chr(0)
	di = DataInterface()
	file_list = get_file_list()
	for file in file_list:
		print(file)
		di.connect(file)
		Fs = di.get_Fs()
		rat_ids = di.get_rat_ids()
		embed()
		for mic_num, id in enumerate(rat_ids):
			if id == null or id == 'V6':
				break

			ji.change_rat(id)

			num_calls = di.get_num_calls(mic_num)
			for signal_num in range(num_calls):
				q = di.get_quality(mic_num, signal_num)
				if q not in ['v', ' ']:
					break 
				signal = di.get_signal(signal_num, mic_num)
				#.insert_signal(signal, Fs)

if __name__ == "__main__":
	main()
