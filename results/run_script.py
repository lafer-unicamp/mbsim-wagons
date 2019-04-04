#!/usr/bin/python

import multiprocessing
import subprocess

def worker(i,bush):

	name = 'freqtest-' + bush + str(i) + 'Hz'

	# Read in the file
	with open('../input2.in', 'r') as file :
		filedata = file.read()

	# Replace the target string
	filedata = filedata.replace('SYSTEM_NAME=test', 'SYSTEM_NAME=' + name)
	filedata = filedata.replace('FREQUENCY=10', 'FREQUENCY=' + str(i))
	if bush == '3d':
		filedata = filedata.replace('BOLSTER_BUSHINGS=1', 'BOLSTER_BUSHINGS=1')
	else:
		filedata = filedata.replace('BOLSTER_BUSHINGS=1', 'BOLSTER_BUSHINGS=0')

	# Write the file out again
	with open('../' + name + '.in', 'w') as file:
		file.write(filedata)

	subprocess.call(['../bin/test', '../' + name + '.in'])

	print('Fim\n')

   	return

if __name__ == '__main__':
    jobs = []
    freqs = range(1,30,5)
    for i in freqs:
        p1 = multiprocessing.Process(target=worker,args=(i,'3d'))
	p2 = multiprocessing.Process(target=worker,args=(i,'1d'))
        jobs.append([p1,p2,])
        p1.start()
	p2.start()
