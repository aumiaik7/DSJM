#!/bin/bash
	#This script runs for every .mtx file in "data_heuristics" folder
	#For each file it creates a fileName.out file and prints the results
	#of the heuristic ordering and paritioning algorithms (time and color)	
	for file in data_heuristics/*.mtx
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

		>results_heuristics/${name[1]}.out
		
			echo -e "Output format: File name, no. of rows, no. of columns, no. of non-zero elements, method name, ordering time, coloring time, colors/partitions \n\n">>results_heuristics/${name[1]}.out 
			#LFO
			echo -e "Method: LFO">>results_heuristics/${name[1]}.out  
			../examples/gcolor -i "$file" -m lfo>>results_heuristics/${name[1]}.out 
			echo -e "\n">>results_heuristics/${name[1]}.out
			#SLO 
			echo -e "Method: SLO">>results_heuristics/${name[1]}.out  
			../examples/gcolor -i "$file" -m slo>>results_heuristics/${name[1]}.out 
			echo -e "\n">>results_heuristics/${name[1]}.out
			#IDO
			echo -e "Method: IDO">>results_heuristics/${name[1]}.out  
			../examples/gcolor -i "$file" -m ido>>results_heuristics/${name[1]}.out 
			echo -e "\n">>results_heuristics/${name[1]}.out
			#SDO
			echo -e "Method: SDO">>results_heuristics/${name[1]}.out  
			../examples/gcolor -i "$file" -m sdo>>results_heuristics/${name[1]}.out 
			echo -e "\n">>results_heuristics/${name[1]}.out
			#RLF
			echo -e "Method: RLF">>results_heuristics/${name[1]}.out  
			../examples/gcolor -i "$file" -m rlf>>results_heuristics/${name[1]}.out 
			echo -e "\n">>results_heuristics/${name[1]}.out

	done
	echo "Done!!"

	
	
 

