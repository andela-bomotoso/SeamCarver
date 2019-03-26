seamcarver:	seamcarver_pthreads.cpp
		g++ seamcarver_pthreads.cpp -o seamcarver_pthreads `pkg-config --cflags lqr-1` `freetype-config --cflags`\
		-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0`  
seamcarver_serial:	seamcarver_serial.cpp
			g++ seamcarver_serial.cpp -o seamcarver_serial `pkg-config --cflags lqr-1` `freetype-config --cflags`\
			-lPNGwriter -lpng -lz -lfreetype  `pkg-config --libs lqr-1` `pkg-config --libs gthread-2.0` 
		

