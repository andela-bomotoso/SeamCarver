#include <pngwriter.h>
#include <lqr.h>
#include <getopt.h>
#include "liquidrescale.h"
#include <math.h>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <mpi.h>


using namespace std;
pngwriter pngwrt(1,1,0,"out_MPI.png");
int BASE_ENERGY = 1000;
int ROWS_PER_THREAD = 64;
int width = 0;
int height = 0;
gfloat rigidity = 0;
gint max_step = 1;
int channels = 3;
int** energyArray;
guchar* seams;
guchar* buffer;
int* verticalSeams;
int** distTo;
int** edgeTo;
int* flattenedEnergyArray;
int* flattenedDistTo;
int* flattenedEdgeTo;
MPI_Win win;
MPI_Win win1;
MPI_Win win2;

/*Copied from the liblqr example
 convert the image in the right format */
guchar * rgb_buffer_from_image(pngwriter *png)
{
    gint x, y, k, channels;
    gint w, h;
    guchar *buffer;

    /* get info from the image */
    w = png->getwidth();
    h = png->getheight();
    channels = 3;                       // we assume an RGB image here 

    /* allocate memory to store w * h * channels unsigned chars */
    buffer = g_try_new(guchar, channels * w * h);
    g_assert(buffer != NULL);

    /* start iteration (always y first, then x, then colours) */
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            for (k = 0; k < channels; k++) {
                /* read the image channel k at position x,y */
                buffer[(y * w + x) * channels + k] = (guchar) (png->dread(x + 1, y + 1, k + 1) * 255);
                /* note : the x+1,y+1,k+1 on the right side are
                 *        specific the pngwriter library */
            }
        }
    }

    return buffer;
}


int computeEnergy(int x, int y, guchar* buffer){
	 if (x == 0 || y == 0 || (x == width - 1) || (y == height- 1))
            return BASE_ENERGY;

	//Declare values for the RGB top, bottom, left and right neighbours
	double top_RB, top_GB, top_BB = 0.0;
	double bottom_RB, bottom_GB, bottom_BB = 0.0;
	double left_RB, left_GB, left_BB = 0.0;
	double right_RB, right_GB, right_BB = 0.0;

	/*Declare values for the differences of the vertical and horizontal 
		RGB neighbours difference*/
	double redDiffH, greenDiffH, blueDiffH = 0.0;
	double redDiffV, greenDiffV, blueDiffV = 0.0;

	/*Declare values for the sum of the vertical and horizontal RGB 	
		neighbours differences*/
	double valueH, valueV,  valueSum;


	//Get the RGB values of the top neighbor of the current pixel

	top_RB = buffer[((y-1)*width+(x))*3+0];
	top_GB = buffer[((y-1)*width+(x))*3+1];
	top_BB = buffer[((y-1)*width+(x))*3+2];

	//Get the RGB values of the bottom neighbor of the current pixel
	bottom_RB = buffer[((y+1)*width+(x))*3+0];
	bottom_GB = buffer[((y+1)*width+(x))*3+1];
	bottom_BB = buffer[((y+1)*width+(x))*3+2];

	//Get the RGB values of the right neigbor of the current pixel
	right_RB = buffer[((y)*width+(x+1))*3+0];
	right_GB = buffer[((y)*width+(x+1))*3+1];
	right_BB = buffer[((y)*width+(x+1))*3+2];

	//Get the RGB values of the left neigbor of the current pixel
	left_RB = buffer[((y)*width+(x-1))*3+0];
	left_GB = buffer[((y)*width+(x-1))*3+1];
	left_BB = buffer[((y)*width+(x-1))*3+2];


	//Get the absolute difference of the red horizontal neighbours
	redDiffV  = abs(top_RB - bottom_RB);

	//Get the absolute difference of the red horizontal neighbours
	redDiffH  = abs(right_RB - left_RB);

	//Get the absolute difference of the green horizontal neighbours
	greenDiffV  = abs(top_GB - bottom_GB);

	//Get the absolute difference of the green horizontal neighbours
	greenDiffH  = abs(right_GB - left_GB);

	//Get the absolute difference of the blue horizontal neighbours
	blueDiffV  = abs(top_BB - bottom_BB);

	//Get the absolute difference of the blue horizontal neighbours
	blueDiffH  = abs(right_BB - left_BB);

	//Get ths sum of the square of differences of vertical neighbours
	valueH = pow(redDiffH,2) + pow(greenDiffH,2) + pow(blueDiffH,2);

	//Get ths sum of the square of differences of horizontal neighbours
	valueV = pow(redDiffV,2) + pow(greenDiffV,2) + pow(blueDiffV,2);
	valueSum = valueV + valueH;

	//Return the squareroot of the sum of differences
	return round(sqrt(valueSum));
}

//Build memory structure for energy matrix
int** initializeEnergyArray(int rows, int columns){
	int** energyArray;
	energyArray = (int**) malloc(rows*sizeof(int*));
	for (int i = 0; i < rows; i++)
		energyArray[i] = (int*) malloc(columns*sizeof(int));
	return energyArray;
}


int* initialize1DArray(int rows, int columns){
        int* array1D = (int*)malloc(rows*columns*sizeof(int));
	 for ( int i = 0; i < rows*columns; i++)
                array1D[i] = 0;
        return array1D;
}


/*int* initializeDistTo1D(int rows, int columns){
        int* array1D = (int*)malloc(rows*columns*sizeof(int));
         for ( int i = 0; i < rows*columns; i++)
		int col = getCol(i, columns);
		int row = getRow(i, columns);
		if (row == 0)
                	array1D[i] = BASE_ENERGY;
      		else 
			array1D[i] = std::numeric_limits<int>::max();

        return array1D;
}*/

int** initializeEdgeTo(int rows, int columns){
        int** edgeTo;
        edgeTo = (int**) malloc(rows*sizeof(int*));
        for (int i = 0; i < rows; i++)
                edgeTo[i] = (int*) malloc(columns*sizeof(int));
        return edgeTo;
}

int** initializeDistTo(int rows, int columns){
        int** distTo;
        distTo = (int**) malloc(rows*sizeof(int*));
        for (int i = 0; i < rows; i++)
                distTo[i] = (int*) malloc(columns*sizeof(int));
	for (int row = 0; row < rows; row++) {
            for (int col = 0; col < columns; col++) {
                if (row == 0)
                    distTo[row][col] = BASE_ENERGY;
                else
                    distTo[row][col] = std::numeric_limits<int>::max();
            }
       }

        return distTo;
}

//Convert a 2D array to 1D

int* flattenArray(int** arrayToFlatten, int num_rows, int num_cols){
	int* flattenedArray = (int*) malloc(num_rows*num_cols*sizeof(int));
	int counter = 0;
	
	for (int row = 0; row < num_rows; row++){
		for(int col = 0; col < num_cols; col++){
			flattenedArray[counter] = arrayToFlatten[row][col];
			counter++;
		}
	}
		return flattenedArray;
}

//Convert a 1D array to 2D
int** unflattenArray(int* flattenedArray, int num_rows, int num_cols){
	int** unflattenedArray =  (int**) malloc(num_rows*sizeof(int*));
        for (int i = 0; i < num_rows; i++)
                unflattenedArray[i] = (int*) malloc(num_cols*sizeof(int));
	int counter = 0;
	for (int row = 0; row < num_rows; row++){
		for (int col = 0; col < num_cols; col++){
		
			unflattenedArray[row][col] = flattenedArray[counter];
			counter++;
		}
	}
	return unflattenedArray;
}

//Given a 1D array, find the row of a given cell if it were to be a 2D array
int getRow(int cell, int num_cols){
		return (cell/num_cols);
}

//Given a 1D array, find the col of a given cell if it were to be a 2D array
int getCol(int cell, int num_cols){
	return cell%num_cols;
}

int* initializeDistTo1D(int rows, int columns){
        int* array1D = (int*)malloc(rows*columns*sizeof(int));
         for ( int cell = 0; cell < rows*columns; cell++){
               // int col = getCol(i, columns);
                int row = getRow(cell, columns);
                if (row == 0)
                        array1D[cell] = BASE_ENERGY;
                else
                        array1D[cell] = std::numeric_limits<int>::max();
	}
        return array1D;
}


void generateEnergyMatrix(int num_rows, int start_col, int stop_col, int me, int num_pes, char* orientation){
    for (int row = 1; row < num_rows; row++)    {
                for(int column = start_col; column < stop_col; column++)    {
			int cell = row*width+column;
                        if (orientation[0] == 'v')
                                 flattenedEnergyArray[cell] = computeEnergy(column, row, buffer);
                        else
                                flattenedEnergyArray[cell] = computeEnergy(row, column, buffer);
		
	/*if (cell == 1900 || cell == 1250){
		cout<< "Printing values inside the energy generation matrix" <<endl;
        	cout<<"Cell "<<cell<< " processed by PE: "<<me<<" has a value of: "<<flattenedEnergyArray[cell]<<endl;
	}*/

		/*Put the new value in PE 0*/
		MPI_Win_lock_all(0, win);
			if(me > 0)
				MPI_Put(&flattenedEnergyArray[cell], 1, MPI_INT, 0, cell, 1, MPI_INT, win);
		MPI_Win_flush_all(win);
		MPI_Win_unlock_all(win);
        }
    }
}
                                                                               


int* backTrack(int** edgeTo, int** distTo, int height, int width){
// Backtrack from the last row to get a shortest path
	int* seams = new int[height];
        int minCol = 0;
        int minDist = std::numeric_limits<int>::max();
        for (int col = 0; col < width; col++) {
            if (distTo[height - 1][col] < minDist) {
                minDist = distTo[height - 1][col];
                minCol = col;
            }
        }
	 for (int row = height - 1; row >= 0; row--) {
            verticalSeams[row] = minCol;
            minCol = minCol - edgeTo[row][minCol];
        }
	return verticalSeams;
	
}

//A Transpose Matrix that makes findVerticalSeams re-usable
 guchar*  transposeRGBuffer(guchar* buffer, int width, int height) {
	guchar* transposedRGBuffer;
	int size = 3 * width * height;
    	transposedRGBuffer = g_try_new(guchar,size);
    	g_assert(transposedRGBuffer != NULL);
	for(int col = 0; col < width; col++){
		for(int row = 0; row < height; row++){
			for (int color = 0; color < 3; color++){
				transposedRGBuffer[(col*height + row)*3 + color] = buffer[(row*width+col)*3+color];
			}			
		}
	}
		return transposedRGBuffer;
}

/*Declare a relax function to optimize the computation of a
 *         shortest path energy values*/

void relax(int row, int col, int width, int start_col, int stop_col, int me, int npes) {
        int relax = 0;
        int nextRow = row + 1;
        int cell = row*width + col;
        for (int i = -1; i <= 1; i++) {
            int nextCol = col + i;
            int nextCell = nextRow*width + nextCol;
            if (nextCol < 0 || nextCol >= width)
                continue;

		if (flattenedDistTo[nextCell] >= flattenedDistTo[cell] + flattenedEnergyArray[nextCell]){
			flattenedDistTo[nextCell] = flattenedDistTo[cell] + flattenedEnergyArray[nextCell];
			flattenedEdgeTo[nextCell] = i;

                /*Put the new values in PE 0, other PEs will fetch these values later*/
                
			MPI_Win_lock_all(0, win1);
			if (me > 0){
                        	MPI_Put(&flattenedDistTo[nextCell], 1, MPI_INT, 0, nextCell, 1, MPI_INT, win1);
			}
			MPI_Win_flush_all(win1);
		        MPI_Win_unlock_all(win1);		
			
			MPI_Win_lock_all(0, win2);
			if (me > 0){
                        	MPI_Put(&flattenedEdgeTo[nextCell], 1, MPI_INT, 0, nextCell, 1, MPI_INT, win2);
			}
			MPI_Win_flush_all(win2);
			MPI_Win_unlock_all(win2);
			
                }
		
        }
    }



//Find the indexes of the vertical seams to be carved out
int * identifySeams( int width, int height, int start_col, int stop_col, int me, int npes){
	 for (int row = 0; row < height - 1; row++) {
             for (int col = start_col; col < stop_col; col++) {
                relax(row, col, width, start_col, stop_col, me, npes);
		
            }
		MPI_Barrier(MPI_COMM_WORLD);
		int nextCell = (row+1)*width;

		/*Other PEs get the new values for the next row from PE 0*/

		if (npes > 1){
                	MPI_Bcast(&flattenedDistTo[nextCell], width, MPI_INT, 0, MPI_COMM_WORLD);
		
                	MPI_Bcast(&flattenedEdgeTo[nextCell], width, MPI_INT, 0, MPI_COMM_WORLD);
        	}	
		MPI_Barrier(MPI_COMM_WORLD);	
	}
	             
}



//Carve out the vertical seams
guchar* carveVertically(int* vertical_seams, guchar* buffer, int width, int height){
	guchar* carved_imageV;
	int size = 3 * width * height;
	
	/* allocate memory to store w * h * channels unsigned chars */
    	carved_imageV = g_try_new(guchar, size);
    	g_assert(carved_imageV != NULL);
	seams = g_try_new(guchar, size);
    	g_assert(seams != NULL);


	for (int row = 0; row < height; row++){
		//Get the RGB value before the seam
		for (int col = 0; col < vertical_seams[row]; col++){
			for (int color = 0; color < 3; color++){
				carved_imageV[(row * width + col) * 3 + color] = buffer[(row * width + col) * 3 + color];
				seams[(row * width + col) * 3 + color] = buffer[(row * width + col) * 3 + color];
			}
		}

		//Get the RGB values after the seams
		for (int col = vertical_seams[row]; col < width-1; col++) {
			for (int color = 0; color < 3; color++){
				carved_imageV[(row * width + col) * 3 + color] = buffer[(row * width + col+1) * 3 + color];
				seams[(row * width + col+1) * 3 + color] = buffer[(row * width + col+1) * 3 + color];
			}
       	 	}
	}

	return carved_imageV;
}

LqrRetVal printSeams(LqrCarver *carver, pngwriter *pngwrt){

    gint x, y;
    guchar *rgb;
    gdouble red, green, blue;

    lqr_carver_scan_reset(carver);

    /* readout (no need to init rgb) */
    while (lqr_carver_scan(carver, &x, &y, &rgb)) {
        /* convert the output into doubles */
        red = (gdouble) rgb[0] / 255;
        green = (gdouble) rgb[1] / 255;
        blue = (gdouble) rgb[2] / 255;

	if(!rgb[0] && !rgb[1] && !rgb[2])
		 pngwrt->plot(x + 1, y + 1, 1.0, 0.0, 0.0);
	else 
		 pngwrt->plot(x + 1, y + 1, red, green, blue);
			
}
	    return LQR_OK;
}

LqrRetVal write_carver_to_image(LqrCarver *carver, pngwriter *pngwrt, char* orientation){

	 gint x, y;
    	 guchar *rgb;
	 gdouble red, green, blue;
    /* resize the image canvas as needed to
     * fit for the new size */
   

	//Resize based on the orientation
	if(orientation[0] == 'v'){
		TRAP(lqr_carver_resize(carver, width-1 , height));
    	 	pngwrt->resize(width-1, height);
	}
	else{
		TRAP(lqr_carver_resize(carver, width , height-1));
		pngwrt->resize(width,height-1);
		}
	 lqr_carver_scan_reset(carver);
	

    /* readout (no need to init rgb) */
    while (lqr_carver_scan(carver, &x, &y, &rgb)) {
        /* convert the output into doubles */
        red = (gdouble) rgb[0] / 255;
        green = (gdouble) rgb[1] / 255;
        blue = (gdouble) rgb[2] / 255;

        /* plot (pngwriter's coordinates start from 1,1) */
        pngwrt->plot(x + 1, y + 1, red, green, blue);
	}
	    return LQR_OK;
}

/* Get current time*/

double timestamp()
{
    struct timeval tval;
    
    gettimeofday( &tval, ( struct timezone * ) 0 );
    return ( tval.tv_sec + (tval.tv_usec / 1000000.0) );
}

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);
	int me, npes;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &me); 

	double total_begin = timestamp();
	
	char * original_img = argv[1]; 
	char* orientation = argv[2];
	
	/*read the image and get the dimension*/

	pngwrt.readfromfile(original_img);
	width = pngwrt.getwidth();
	height = pngwrt.getheight();

	int cells_num = height*width;
	//flattenedEnergyArray = initialize1DArray(height, width);
	
	// Allocate a window for each buffer
	MPI_Win_allocate(cells_num*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &flattenedEnergyArray, &win); 
	MPI_Win_allocate(cells_num*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &flattenedDistTo, &win1);
	//flattenedDistTo = initializeDistTo1D(height, width);
	
	MPI_Win_allocate(cells_num*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &flattenedEdgeTo, &win2);

        if (me == 0)
		cout<<"Width: "<<width<<" Height: "<<height<<endl;

	double begin, end;

	int size = 3 * width * height;
    	buffer = g_try_new(guchar,size);
	buffer = rgb_buffer_from_image(&pngwrt);
	
    	g_assert(buffer != NULL);
	

	LqrCarver *carver;
	LqrCarver *carved_seams;
	//Check the orientation to determine how to carve
	begin = timestamp();
		
	if(orientation[0] == 'v'){
	
		verticalSeams = new int[height];
		
		/*initialize energy 1D and 2D Array*/
		energyArray = initializeEnergyArray(height, width);
		flattenedEnergyArray = initialize1DArray(height, width);
	
		/*initialize distTo 1D and 2D arrays*/
		distTo = initializeDistTo(height, width);
		flattenedDistTo = initializeDistTo1D(height, width);
         	
		/*initialize edgeTo 1D and 2D arrays*/

		edgeTo = initializeEdgeTo(height, width);
		flattenedEdgeTo = initialize1DArray(height, width);

		/*Divide the work among the available PEs*/
		int cell_per_pe = (width*height + npes - 1)/npes;
		int col_per_per = (width + npes - 1)/npes;
		
		int start_cell = me * cell_per_pe;
		int stop_cell = (me + 1) * cell_per_pe;

		int start_col = me*col_per_per;
		int stop_col = (me+1) * col_per_per;

	/*Ensure that the last pe does not exceed the last row*/
	if (me == npes - 1)
		stop_cell  = width*height;
	
	if (me == npes - 1 )
                stop_col  = width;

	int start_col_en =  start_col;
	int stop_col_en = stop_col;
	
	//consider adjacent pixels for energy computation
	if(me > 0 )
			start_col_en = start_col - 1;
		if (me < npes - 1)
			stop_col_en = stop_col + 1;

		
		//Fill the energy matrix with the energy values of each pixel

		generateEnergyMatrix(height, start_col_en, stop_col_en, me, npes, orientation);
	
		MPI_Barrier(MPI_COMM_WORLD);

		/* PE 0 broadcasts the values of the 1D energy Array*/ 
 		if (npes > 1)
			MPI_Bcast(&flattenedEnergyArray[0], width*height, MPI_INT, 0, MPI_COMM_WORLD);
		 
		MPI_Barrier(MPI_COMM_WORLD);

	       /*cout<< "Printing values after PE 0 does a broadcast" <<endl;
                cout<<"Cell 1250 on PE: "<<me<<" has a value of: "<<flattenedEnergyArray[1250]<<endl;
		cout<<"Cell 1900 on PE: "<<me<<" has a value of: "<<flattenedEnergyArray[1900]<<endl;	*/

		identifySeams(width, height, start_col, stop_col, me, npes);
		if (me == 0){
			distTo = unflattenArray(flattenedDistTo, height, width);
			edgeTo = unflattenArray(flattenedEdgeTo, height, width);
                cout<<"Removing vertical seams"<<endl;

		int* v_seams =  backTrack(edgeTo, distTo, height, width);
		
		/*PE 0 will remove the identified seams
		  there is no need to parallelize this process*/
			guchar* carved_imageV = carveVertically(v_seams,buffer, width, height);
			carver = lqr_carver_new(carved_imageV, width, height, 3);
			carved_seams = lqr_carver_new(seams, width, height, 3);

			end = timestamp();
			cout<<"Total Seam Carving Time: "<<(end-begin)<<endl; 
		}
	}
	else{
		verticalSeams = new int[width];

		//initialize  energy array	
		energyArray = initializeEnergyArray(width, height);
		
		/*divide the work among the available PEs*/
		int col_per_pe = (height + npes - 1)/npes;
		int start_col = me * col_per_pe;
		int stop_col = (me + 1) * col_per_pe;
		
		/*Ensure the last PE does not exceed the last row*/
		if (me == npes - 1)
			stop_col = height;
		 
		/*consider adjacent pixels for energy computation*/  
		int start_col_en = start_col;
		int stop_col_en = stop_col;

                if(me > 0 )  
                      start_col_en = start_col - 1;
               
		if (me < npes - 1)
                        stop_col_en = stop_col + 1;
                
		/*Fill the energy matrix with the energy values of each pixel  */                                                       generateEnergyMatrix(width, start_col_en, stop_col_en, me,npes, orientation);
                //shmem_barrier_all();
			
		/*Find the horizontal seams in the image*/
		identifySeams(height, width, start_col, stop_col, me, npes);
		 if (me == 0){
                        cout<<"Removing horizontal seams"<<endl;

		int* h_seams =  backTrack(edgeTo, distTo, width, height);
		
		/*Revome the identified horizontal seams*/
			guchar* transBuffer = transposeRGBuffer(buffer, width, height);
			guchar* carved_imageH = carveVertically(h_seams, transBuffer, height, width);
			carver = lqr_carver_new(transposeRGBuffer(carved_imageH, height,width), width, height, 3);
			carved_seams = lqr_carver_new(transposeRGBuffer(seams, height,width), width, height, 3);
		 
		 end = timestamp();
                 cout<<"Total Seam Carving Time: "<<(end-begin)<<endl;
		}
	}

	//Create a Carver object with the carved image buffer
   	if (me == 0){
		TRAP(lqr_carver_init(carver, max_step, rigidity));
	//write_carver_to_image(carver, &pngwrt, orientation);
		printSeams(carved_seams, &pngwrt);
		lqr_carver_destroy(carver);
		pngwrt.close();

		double total_end = timestamp();
		printf("%s%5.2f\n","Total Processing Time: ", (total_end-total_begin));
	}
	MPI_Win_free(&win);
	MPI_Finalize();
	return 0;
	
}

