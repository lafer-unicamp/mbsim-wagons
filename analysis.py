#!-*-coding: utf8 -*-

from subprocess import call
import os
import shutil

def writeInputFile(om,t0,amp,name):
	# opens file
	f = open('input.in','w')
	print f

	# writes overall controls
	f.write("%======================== OVERALL CONTROLS ===================================%\n")
	f.write("SIMULATION_TIME="+str(10)+"\n")
	f.write("SYSTEM_NAME="+name+"\n")
	# writes motion parameters
	f.write("%======================== MOTION INPUT =======================================%\n")
	f.write("AMPLITUDE="+str(amp)+"\n")
	f.write("ANGULAR_VELOCITY="+str(om)+"\n")
	f.write("MOVEMENT_DELAY="+str(t0)+"\n")
	# writes geometric parameter
	f.write("%======================== WAGON GEOMETRIC DATA ===============================%\n")
	f.write("CAR_WHEEL_BASE="+str(12.3)+"\n")
	f.write("TRUCK_SEMI_WHEEL_BASE="+str(1.575)+"\n")
	#writes inertial parameters
	f.write("BOX_MASS="+str(65134.)+"\n")
	f.write("BOX_Ixx="+str(52.8e3)+"\n")	# longitudinal
	f.write("BOX_Iyy="+str(1284.e3)+"\n")	# vertical
	f.write("BOX_Izz="+str(1290.e3)+"\n")	# transversal
	f.closed
	print "The file has been closed"


for i in range(4):
	executable = "/home/leonardo/Projects/mbsim/mbsim-friction-wedge/main"
	pathName = "Wagon"+str(i)			# simulation name
	rootDir = "./"
	destinationDir = "./"+pathName+"/"
	# creates analysis directory
	try:
		os.mkdir(destinationDir)
	except OSError:
		# if the directory already existis, asks the user to delete it
		print "O diretório "+destinationDir+"já existe.\n"
		userInput = raw_input("Apagar o existente e criar um novo? [y/N]")
		if (userInput == "y" or userInput == "Y"):
			shutil.rmtree(destinationDir)
			os.mkdir(destinationDir)
		else:
			# if the user choose not to delete the directory, exits the routine
			break
	os.chdir(destinationDir)
	writeInputFile(15.9*(1+i)/4,2,0.01,pathName)
	call(executable)
	os.chdir("../")
