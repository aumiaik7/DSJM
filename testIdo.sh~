#!/bin/bash
>Testing_script/testResultsAvg.csv
>testTrack
iterations=1

	for file in /dfs/Data_MM/*.mtx
	#for file in ~/Data_DSJM/*.mtx
	do
	echo "-----------------File Name: $file -------------------">>testTrack
	
	 # do something on "$file"
		#file name is now like data/Matrix.mtx
		#the following lines are used to remove .mtx and /data	  
		#split file name using '.' delimeter 
		IFS='.'
		read -a direc <<< "${file}"
		#split file name using '.' delimeter 
		IFS='/'
		read -a name <<< "${direc}"
>tempTest
##################  LFO
#		for ((i=0;i<iterations;i++));
#       		do
#			examples/gcolor -i "$file" -m lfo >> tempTest 	    
#    		done 
#		total_ordering_time=0
#		total_coloring_time=0
#		total_time=0
#		while read -r line1
#		do 
#			#split line using '='
#			IFS=','
#			read -a split_res <<< "${line1}"
#	
#			ordering_time=${split_res[5]}
#			coloring_time=${split_res[6]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
#			total_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}+${ordering_time}; exit(0)}")
#			total_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
#		done < tempTest
#		avg_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}/${iterations}; exit(0)}")
#		avg_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}/${iterations}; exit(0)}")
#		total_time=$(awk "BEGIN {printf \"%.6f\",${avg_ordering_time}+${avg_coloring_time}; exit(0)}")
#		echo "$file"",""${split_res[1]}"",""${split_res[2]}"",""${split_res[3]}"",""${split_res[4]}"","$avg_ordering_time","$avg_coloring_time","$total_time",""${split_res[7]}">>Testing_script/testResultsAvg.csv
>tempTest
##################  SLO
		for ((i=0;i<iterations;i++));
       		do
			examples/gcolor -i "$file" -m slo >> tempTest 	  
			echo "SLO Iteration $i">>testTrack
    		done 
		total_ordering_time=0
		total_coloring_time=0
		total_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[5]}
			coloring_time=${split_res[6]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}/${iterations}; exit(0)}")
		total_time=$(awk "BEGIN {printf \"%.6f\",${avg_ordering_time}+${avg_coloring_time}; exit(0)}")
		echo "$file"",""${split_res[1]}"",""${split_res[2]}"",""${split_res[3]}"",""${split_res[4]}"","$avg_ordering_time","$avg_coloring_time","$total_time",""${split_res[7]}">>Testing_script/testResultsAvg.csv

		#echo "-------------------Done SLO for $file">>testTrack
>tempTest
##################  IDO
		for ((i=0;i<iterations;i++));
       		do
			examples/gcolor -i "$file" -m ido >> tempTest 	
			echo "IDO Iteration $i">>testTrack    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		total_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[5]}
			coloring_time=${split_res[6]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}/${iterations}; exit(0)}")
		total_time=$(awk "BEGIN {printf \"%.6f\",${avg_ordering_time}+${avg_coloring_time}; exit(0)}")
		echo "$file"",""${split_res[1]}"",""${split_res[2]}"",""${split_res[3]}"",""${split_res[4]}"","$avg_ordering_time","$avg_coloring_time","$total_time",""${split_res[7]}">>Testing_script/testResultsAvg.csv

		#echo "-------------------Done IDO for $file">>testTrack

>tempTest
##################  RLF
		for ((i=0;i<iterations;i++));
       		do
			examples/gcolor -i "$file" -m rlf >> tempTest 	 
			echo "RLF Iteration $i">>testTrack   
    		done 
		total_ordering_time=0
		total_coloring_time=0
		total_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[5]}
			coloring_time=${split_res[6]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}/${iterations}; exit(0)}")
		total_time=$(awk "BEGIN {printf \"%.6f\",${avg_ordering_time}+${avg_coloring_time}; exit(0)}")
		echo "$file"",""${split_res[1]}"",""${split_res[2]}"",""${split_res[3]}"",""${split_res[4]}"","$avg_ordering_time","$avg_coloring_time","$total_time",""${split_res[7]}">>Testing_script/testResultsAvg.csv

		#echo "-------------------Done RLF for $file">>testTrack

>tempTest
##################  SDO
		for ((i=0;i<iterations;i++));
       		do
			examples/gcolor -i "$file" -m sdo >> tempTest 
			echo "SDO Iteration $i">>testTrack	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		total_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[5]}
			coloring_time=${split_res[6]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.6f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.6f\",${total_coloring_time}/${iterations}; exit(0)}")
		total_time=$(awk "BEGIN {printf \"%.6f\",${avg_ordering_time}+${avg_coloring_time}; exit(0)}")
		echo "$file"",""${split_res[1]}"",""${split_res[2]}"",""${split_res[3]}"",""${split_res[4]}"","$avg_ordering_time","$avg_coloring_time","$total_time",""${split_res[7]}">>Testing_script/testResultsAvg.csv

		#echo "-------------------Done SDO for $file">>testTrack
		echo -e "-------------Done Ordering and Coloring for $file -----------------\n">>testTrack

	done
	echo "Done!!"

	#echo $1
	
	
 

