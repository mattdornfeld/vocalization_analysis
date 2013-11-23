#! /usr/bin/env python

import sqlite3 as lite
import numpy as np
import io
from IPython import embed

"""
This a class for interacting with jump.db, thee database that stores rat
calls, fitted curves, and jump points.
"""

class JumpInterface:
	def __init__(self, db_path, rat=None):			
		lite.register_adapter(np.ndarray, self.__array_to_text)
		lite.register_converter("ARRAY", self.__text_to_array)
		self.con = lite.connect(db_path, detect_types=lite.PARSE_DECLTYPES)

		self.rats = self.get_rats()
		if rat != None:
			self.change_rat(rat)

	def get_rats(self):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT rat FROM signals")
			rats = np.array(cur.fetchall()).flatten()
			
			return np.unique(rats)
		
	def change_rat(self, rat):
		if rat[0] != 'V':
			raise Exception("Invalid rat name.")
		if unicode(rat) not in self.rats:
			raise Warning("This rat doesn't exist yet.")

		self.rat = rat

	def get_jump(self, jump_index):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT f1, f2, t1, t2 FROM jumps WHERE rat=? AND \
				jump_index=?", (self.rat, jump_index))
			jump = np.array(cur.fetchone()).flatten()
			
			return jump
	
	def get_jumps(self, quality=1):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT f1, f2, t1, t2 FROM jumps WHERE rat=? AND \
			 quality=? ", (self.rat, quality))
			jumps = np.array(cur.fetchall())
			
			return jumps
		
	def get_clusters(self, quality=1):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT cluster FROM jumps WHERE rat=? AND \
				quality=?", (self.rat, quality))
			
			cluster = np.array([c[0] for c in cur.fetchall()])					
			
			return cluster
	
	def get_signal(self, signal_index):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT signal FROM signals WHERE signal_index = ? ", (signal_index,))
			
			signal = cur.fetchone()
			if signal == None:
				raise Exception("There is no signal stored at signal_index = " + str(signal_index))
			
			return signal[0]
			
	def get_curve(self, signal_index):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT curve FROM curves WHERE rat = ? AND signal_index = ? ", (self.rat, signal_index))
			
			curve = cur.fetchone()
			if curve == None:
				raise Exception("There is no curve stored at signal_index = " + str(signal_index))

			return curve[0]
	
	def get_jump_indices(self, signal_index):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT jump_index FROM jumps WHERE rat = ? AND signal_index = ?", (self.rat, signal_index))          
			jump_indices = np.array(cur.fetchall())
			jump_indices.shape = (len(jump_indices),)
			if len(jump_indices) == 0:
				raise Exception("There are no jumps associated with signal_index = " + str(signal_index))
			
			return jump_indices 

	def get_jump_indices_quality(self, quality):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT jump_index FROM jumps WHERE rat = ? AND quality = ?", (self.rat, quality))          

			jump_indices = np.array(cur.fetchall())
			jump_indices.shape = (len(jump_indices),)
			
			return jump_indices

	def get_signal_index(self, jump_index):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT signal_index FROM jumps WHERE rat = :rat AND jump_index = :jump_index", {"rat": self.rat, "jump_index": jump_index})			
			
			signal_index = cur.fetchone()
			if signal_index == None:
				return signal_index
			else:
				return signal_index[0]

	def get_signal_indices_rat(self):
		with self.con:
			cur = self.con.cursor()
			cur.execute("SELECT signal_index FROM signals WHERE rat=?", \
				(self.rat,))

			return np.array(cur.fetchall()).flatten()
			
	def get_signal_indices(self, quality=1):
		with self.con:
			cur = self.con.cursor()
			
			cur.execute("SELECT signal_index FROM jumps WHERE rat = ? \
				AND quality=?", (self.rat, quality))			
			
			signal_indices = np.array([i[0] for i in cur.fetchall()])							
		
			return signal_indices
	
	def get_column_names(self):
		with self.con:
			cur = self.con.cursor()      
			cur.execute("PRAGMA table_info(jumps)")
			
			tableInfo = cur.fetchall()
			print('Table Name: jumps')
			for i in tableInfo:
				print i[0], i[1], i[2]
			
			cur.execute("PRAGMA table_info(signals)")
			tableInfo = cur.fetchall()
			print('\nTable Name: signals')
			for i in tableInfo:
				print i[0], i[1], i[2]
			
			cur.execute("PRAGMA table_info(curves)")	
			tableInfo = cur.fetchall()
			print('\nTable Name: curves')
			for i in tableInfo:
				print i[0], i[1], i[2]
	
	def insert_jumps(self, jumps, signal_indices, clusters, quality):
		with self.con:
			cur = self.con.cursor()	
			for index, j in enumerate(jumps):
				i = int(signal_indices[index])
				c = clusters[index]
				q = quality[index]
				cur.execute("INSERT INTO jumps VALUES(? ,?, ?, ?, ?, ?, ?, \
					?)", (self.rat, j[0], j[1], j[2], j[3], i, c, q))
	
	def insert_curve(self, curve):
		if curve.shape != (len(curve),2):
			raise Exception("Array curve must have shape (len(curve),2).")

		with self.con:
			cur = self.con.cursor()
			
			cur.execute("INSERT INTO curves(rat, curve) VALUES(?,?)", \
				(self.rat, curve))
	
	def insert_signal(self, signal, Fs):
		if signal.shape != (len(signal),):
			raise Exception("Array signal must have shape (len(signal),).")

		with self.con:
			cur = self.con.cursor()
			cur.execute("INSERT INTO signals(rat, Fs, signal) \
				VALUES(?,?,?)", (self.rat, Fs, signal))

	def delete_signal(self, signal_index):
		with self.con:
			jump_indices = self.get_jump_indices(signal_index)
			for jump_index in jump_indices:
				self.delete_jump(jump_index)

	def delete_jump(self, jump_index):
		with self.con:
			signal_index = self.get_signal_index(jump_index)
			if signal_index == None:
				raise Exception("I don't think that jump exists.")

			jump_indices = self.get_jump_indices(signal_index)

			cur = self.con.cursor()
			#delete jump and deincrement indices
			cur.execute("DELETE FROM jumps WHERE rat=? AND jump_index=?", (self.rat, jump_index))
			#if that jump is the only one from a given signal, delete the signal and curves
			if len(jump_indices) == 1:
				cur.execute("DELETE FROM signals WHERE rat=? AND signal_index=?", (self.rat, signal_index))
				cur.execute("DELETE FROM curves WHERE rat=? AND signal_index=?", (self.rat, signal_index))				
				
	#tells sqlite how to store numpy arrays
	def __array_to_text(self, arr):
		out = io.BytesIO()
		np.save(out, arr)
		out.seek(0)
		
		return buffer(out.read())

	#tells sqlite how to read numpy arrays
	def __text_to_array(self, text):
		out = io.BytesIO(text)
		out.seek(0)

		return np.load(out)

	def close(self):
		self.con.close()

if __name__ == '__main__':
	dbPath = '/media/matthew/1f84915e-1250-4890-9b74-4bfd746e2e5a/jump.db'
	rat = 'V6'
	if 'jdb' in locals():
		jdb.close()
	jdb = JumpInterface(dbPath, rat)
	
	

