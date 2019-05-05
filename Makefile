all :	seamcarver_pthreads seamcarver_serial seamcarver_SHMEM seamcarver_int64
seamcarver_pthreads:	seamcarver_pthreads.cpp
		g++ -std=c++11 seamcarver_pthreads.cpp -o seamcarver_pthreads `pkg-config --cflags lqr-1` `freetype-config --cflags`\
		-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`  
seamcarver_serial:	seamcarver_serial.cpp
			g++ seamcarver_serial.cpp -o seamcarver_serial `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0` 
seamcarver_SHMEM:	seamcarver_SHMEM.cpp
			g++ seamcarver_SHMEM.cpp -o seamcarver_SHMEM `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`
seamcarver_int64:	seamcarver_int64.cpp
			g++ seamcarver_int64.cpp -o seamcarver_int64 `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`
		

