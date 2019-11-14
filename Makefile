all :	seamcarver_pthreads seamcarver_serial seamcarver_SHMEM seamcarver_int64serial seamcarver_int64pthreads seamcarver_int64SHMEM
seamcarver_pthreads:	seamcarver_pthreads.cpp
			g++ -std=c++11 seamcarver_pthreads.cpp -o seamcarver_pthreads `pkg-config --cflags lqr-1` `freetype-config --cflags`\
		-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`  
seamcarver_serial:	seamcarver_serial.cpp
			g++ -std=c++11 seamcarver_serial.cpp -o seamcarver_serial `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0` 
seamcarver_SHMEM:	seamcarver_SHMEM.cpp
			oshc++ -std=c++11 seamcarver_SHMEM.cpp -o seamcarver_SHMEM `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`
seamcarver_int64serial:	seamcarver_int64serial.cpp
			g++ -std=c++11 seamcarver_int64serial.cpp -o seamcarver_int64serial `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`
seamcarver_int64pthreads:    seamcarver_int64pthreads.cpp
			g++ -std=c++11 seamcarver_int64pthreads.cpp -o seamcarver_int64pthreads `pkg-config --cflags lqr-1` `freetype-config --cflags`\
                -lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`
seamcarver_int64SHMEM:	seamcarver_int64SHMEM.cpp 
			oshc++ -std=c++11 seamcarver_int64SHMEM.cpp -o seamcarver_int64SHMEM `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`
