#!/bin/bash
#>testSdoNew
iterations=1

	#for file in ~/Desktop/data/benchmark/*.mtx
	for file in ~/Desktop/gpia/*.mtx
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


			examples/gcolor -i "$file" -m sdo >>testSdoOld
			
			

	done
	echo "Done!!"

	#echo $1
	
	
 

