#Load Gromos Settings
. ../../gromos_settings.sh


for x in $( ls ${GROMDIR}/testbuild* -d);
do 
	if [[ ${x} == *"testbuild"*"mpi"* ]];
	then
		echo ${x};
		binDIR="${x}/bin"
		dirP=$(basename $x)
		dirP=${dirP/build}
		dirP="out_${dirP}"

		echo ${dirP}
		rm -rf ${dirP}
		cp -r templateDir $dirP

		cd ${dirP}
		./submit_all.sh ${binDIR}  || echo "Failed: ${x} " >> Failed.log
		cd ..
	fi
done
