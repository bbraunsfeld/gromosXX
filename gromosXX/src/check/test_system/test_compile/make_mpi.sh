NAME="testbuild_mpi"

. ../gromos_settings.sh

echo "GROMPATH: ${GROMSRC}"
cd ${GROMSRC}

./Config.sh


mkdir ${NAME} -p
cd ${NAME} 

echo "${PWD}"


../configure --enable-mpi CC=mpiCC CXX=mpiCC > ${NAME}.log  || $(echo "OHOH! ->conf" && exit 1);

make -j8 > ${NAME}.log || $(echo "OHOH! -> MAKE" && exit 1);
make install  > ${NAME}.log  || $(echo "OHOH! -> MAKE install" && exit 1);
