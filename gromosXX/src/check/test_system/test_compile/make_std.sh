NAME="testbuild_standard"

. ../gromos_settings.sh


cd ${GROMSRC}

mkdir ${NAME}

cd ${NAME} 

../configure  > ${NAME}.log  || $(echo "OHOH! ->conf" && exit 1);

make -j8  > ${NAME}.log  || $(echo "OHOH! -> MAKE" && exit 1);
make install  > ${NAME}.log  || $(echo "OHOH! -> MAKE install" && exit 1);
