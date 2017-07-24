#! /bin/bash
ARRAY1=(5. 10. 20. 30. 40. 50. 60. 70. 80. 90. 100. 125. 150. 175. 200. 250. 300. 350.) # peak jet
ARRAY2=(0. 0.1 0.2 1.) # cvis in Smagorinsky model


ELEMENTS1=${#ARRAY1[@]} # elements in first array
ELEMENTS2=${#ARRAY2[@]} # elements in second array

for (( i=0;i<$ELEMENTS1;i++)); do
	for (( j=0;j<$ELEMENTS2;j++)); do
			# Runs with the hm process switched on:
	 		echo ${ARRAY1[${i}]} ${ARRAY2[${j}]} 
			sed -e "s/u_jet=50./u_jet=${ARRAY1[${i}]}/" namelist.in > /tmp/namelist.tmp
 			sed -e "s/cvis=0.2/cvis=${ARRAY2[${j}]}/" /tmp/namelist.tmp > namelist.run
 			
 			#mpiexec -n 64 ./main.exe namelist.run > std.out
 
 			#mv /tmp/output.nc /tmp/output_${i}_${j}.nc 		 			
 			
	done
done

