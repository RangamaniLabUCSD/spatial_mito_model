import os, sys,time
import subprocess as sp

'''
This script was used from a terminal to run all the spatial simulations in MCell
GCG
8.11.23
'''
python_path ='mcell4'
run_file = 'Scene.main.mdl'
tdir = './'
NUM_SEEDS = 10

with open("sim_id.txt", "w+") as f:
	f.write("simulation id \n")

for i in range(1, int(NUM_SEEDS) + 1):
	Proc = sp.Popen([python_path,'-seed',str(i),run_file],cwd = tdir, close_fds=True)
	pid_ = Proc.pid
	with open("sim_id.txt", "a") as f:
		f.write(str(pid_)+ "\t\t" + str(i) + "\n")
