#!/bin/bash
	#This script runs for every .mtx file in "data_exact" folder
	#For each file it creates a fileName.out file and prints the results
	#of the exact paritioning algorithms (four tie-breaking strategies)	
	for file in data_exact/*.mtx
	do
	
	    # do something on "$file"
		#file name is now like data_heuristics/Matrix.mtx
		#the following lines are used to remove .mtx and /data_heuristics	  
		#split file name using '.' delimeter 
		IFS='.'
		read -a direc <<< "${file}"
		#split file name using '.' delimeter 
		IFS='/'
		read -a name <<< "${direc}"

		>results_exact/${name[1]}.out
		
			#echo -e "Output format: File name, no. of rows, no. of columns, no. of non-zero elements, method name, ordering time, coloring time, colors/partitions \n\n">>results_exact/${name[1]}.out 
			#SIMPLE
			echo -e "######################## Tie-breaking method: SIMPLE">>results_exact/${name[1]}.out  
			../examples/gcolor -i "$file" -m exact -t "1">>results_exact/${name[1]}.out 
			echo -e "\n">>results_exact/${name[1]}.out
			#SEWELL 
			echo -e "######################## Tie-breaking method: SEWELL">>results_exact/${name[1]}.out  
			../examples/gcolor -i "$file" -m exact -t "2">>results_exact/${name[1]}.out 
			echo -e "\n">>results_exact/${name[1]}.out
			#SEGUNDO
			echo -e "######################## Tie-breaking method: SEGUNDO">>results_exact/${name[1]}.out  
			../examples/gcolor -i "$file" -m exact -t "3">>results_exact/${name[1]}.out 
			echo -e "\n">>results_exact/${name[1]}.out
			#NEW
			echo -e "######################## Tie-breaking method: NEW">>results_exact/${name[1]}.out  
			../examples/gcolor -i "$file" -m exact -t "4">>results_exact/${name[1]}.out 
			echo -e "\n">>results_exact/${name[1]}.out
			
	done
	echo "Done!!"

	
	
 

