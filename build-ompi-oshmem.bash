#!/bin/bash

git clone https://gitlab.com/hjelmn/xpmem.git 
git clone https://github.com/openucx/ucx
git clone https://github.com/open-mpi/ompi ompi-pmix
git clone https://github.com/libevent/libevent libevent-pmix
git clone https://github.com/pmix/pmix pmix

WD=`pwd`
PREFIX=${WD}/install.shmem-full

cd xpmem; ./autogen.sh; ./configure --prefix=$PREFIX; make -j; make install

cd $WD/ucx; git checkout v1.2.x; ./autogen.sh; ./configure --prefix=$PREFIX --disable-numa --with-xpmem=$PREFIX; make -j; make install

cd $WD/libevent-pmix; ./autogen.sh; ./configure --prefix=$PREFIX; make -j; make install

cd $WD/pmix; git checkout v2.0.2; ./autogen.pl; ./configure --prefix=$PREFIX --with-libevent=$PREFIX; make -j; make install

cd $WD/ompi-pmix; git checkout v3.0.0; ./autogen.pl; ./configure --prefix=$PREFIX --with-ucx=$PREFIX --with-libevent=$PREFIX --with-pmix=$PREFIX; make -j; make install

cd $WD

echo "PATH=${PREFIX}/bin:$PATH" > source-shmem.sh
echo "LD_LIBRARY_PATH=${PREFIX}/lib:$LD_LIBRARY_PATH" >> source-shmem.sh
