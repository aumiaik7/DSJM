#!/bin/bash
>Testing_script/testResultsSloRlf.csv
iterations=5

	for file in data/*.mtx
	do
	
	 # do something on "$file"
		#file name is now like data/Matrix.mtx
		#the following lines are used to remove .mtx and /data	  
		#split file name using '.' delimeter 
		IFS='.'
		read -a direc <<< "${file}"
		#split file name using '.' delimeter 
		IFS='/'
		read -a name <<< "${direc}"


			examples/gcolor -i "$file" -m slo_rlf >>Testing_script/testResultsSloRlf.csv

	done
	echo "Done!!"

	#echo $1
	
	
 

