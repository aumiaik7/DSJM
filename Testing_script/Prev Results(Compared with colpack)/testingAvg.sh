#!/bin/bash
>Testing_script/testResultsAvg.csv
iterations=5

	for file in Graphs/*.mtx
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
################# Column
>tempTest
##################  Natural
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" NATURAL COLUMN_PARTIAL_DISTANCE_TWO >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest

		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv
>tempTest
#################  LFO
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" LARGEST_FIRST COLUMN_PARTIAL_DISTANCE_TWO  >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv
>tempTest
#################  SLO
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" SMALLEST_LAST COLUMN_PARTIAL_DISTANCE_TWO  >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv

>tempTest
#################  IDO
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" INCIDENCE_DEGREE COLUMN_PARTIAL_DISTANCE_TWO  >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv

################# Row
>tempTest
##################  Natural
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" NATURAL ROW_PARTIAL_DISTANCE_TWO >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest

		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv
>tempTest
#################  LFO
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" LARGEST_FIRST ROW_PARTIAL_DISTANCE_TWO  >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv
>tempTest
#################  SLO
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" SMALLEST_LAST ROW_PARTIAL_DISTANCE_TWO  >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv

>tempTest
#################  IDO
		for ((i=0;i<iterations;i++));
       		do
			./ColPack "$file" INCIDENCE_DEGREE ROW_PARTIAL_DISTANCE_TWO  >> tempTest 	    
    		done 
		total_ordering_time=0
		total_coloring_time=0
		while read -r line1
		do 
			#split line using '='
			IFS=','
			read -a split_res <<< "${line1}"
	
			ordering_time=${split_res[3]}
			coloring_time=${split_res[4]}
			#echo $ordering_time
			#echo $coloring_time
			#cumilitive sum of fopt and time 
			total_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}+${ordering_time}; exit(0)}")
			total_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}+${coloring_time}; exit(0)}")
		    #echo "Two lines $line1 and $line2 $total_fopt"	
		done < tempTest
		avg_ordering_time=$(awk "BEGIN {printf \"%.2f\",${total_ordering_time}/${iterations}; exit(0)}")
		avg_coloring_time=$(awk "BEGIN {printf \"%.2f\",${total_coloring_time}/${iterations}; exit(0)}")
		echo "$file"",""${split_res[1]}"","$avg_ordering_time","$avg_coloring_time",""${split_res[5]}">>Testing_script/testResultsAvg.csv


		#./ColPack "$file" LARGEST_FIRST COLUMN_PARTIAL_DISTANCE_TWO >> testResults.csv 
		#./ColPack "$file" SMALLEST_LAST COLUMN_PARTIAL_DISTANCE_TWO >> testResults.csv
		#./ColPack "$file" INCIDENCE_DEGREE COLUMN_PARTIAL_DISTANCE_TWO >> testResults.csv
		#./ColPack "$file" NATURAL ROW_PARTIAL_DISTANCE_TWO >> testResults.csv 
		#./ColPack "$file" LARGEST_FIRST ROW_PARTIAL_DISTANCE_TWO >> testResults.csv 
		#./ColPack "$file" SMALLEST_LAST ROW_PARTIAL_DISTANCE_TWO >> testResults.csv
		#./ColPack "$file" INCIDENCE_DEGREE ROW_PARTIAL_DISTANCE_TWO >> testResults.csv
		#Put gcolor result and profilings in appropriate folders
		
	done
	echo "Done!!"

	#echo $1
	
	
 

