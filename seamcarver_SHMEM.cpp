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
#include <shmem.h>


using namespace std;
pngwriter pngwrt(1,1,0,"out_SHMEM.png");
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
static long pSync[SHMEM_BCAST_SYNC_SIZE];


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
	energyArray = (int**) shmem_malloc(rows*sizeof(int*));
	for (int i = 0; i < rows; i++)
		energyArray[i] = (int*) shmem_malloc(columns*sizeof(int));
	return energyArray;
}

int** initializeEdgeTo(int rows, int columns){
        int** edgeTo;
        edgeTo = (int**) shmem_malloc(rows*sizeof(int*));
        for (int i = 0; i < rows; i++)
                edgeTo[i] = (int*) shmem_malloc(columns*sizeof(int));
        return edgeTo;
}

int** initializeDistTo(int rows, int columns){
        int** distTo;
        distTo = (int**) shmem_malloc(rows*sizeof(int*));
        for (int i = 0; i < rows; i++)
                distTo[i] = (int*) shmem_malloc(columns*sizeof(int));
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


void generateEnergyMatrix(int num_rows, int start_col, int stop_col, int num_pes, char* orientation){
    for (int row = 1; row < num_rows; row++)    {
                for(int column = start_col; column < stop_col; column++)    {
                        if (orientation[0] == 'v')
                                 energyArray[row][column] = computeEnergy(column, row, buffer);
                        else
                                energyArray[row][column] = computeEnergy(row, column, buffer);
		//Put the new value in PE 0
		shmem_int_put(&energyArray[row][column], &energyArray[row][column], 1, 0);
		shmem_quiet();
        }
    }
}
                                                                               

/*Declare a relax function to optimize the computation of a 
	shortest path energy values*/

void relax(int row, int col, int** edgeTo, int** distTo, int width) {
	int relax = 0;
        int nextRow = row + 1;
        for (int i = -1; i <= 1; i++) {
            int nextCol = col + i;
            if (nextCol < 0 || nextCol >= width)
                continue;
            if (distTo[nextRow][nextCol] >= distTo[row][col] + energyArray[nextRow][nextCol]) {
                distTo[nextRow][nextCol] = distTo[row][col] + energyArray[nextRow][nextCol];
                edgeTo[nextRow][nextCol] = i;
		/*Put the new values in PE 0, it will be broadcast later*/
		shmem_int_put(&distTo[nextRow][nextCol], &distTo[nextRow][nextCol], 1, 0);
		shmem_int_put(&edgeTo[nextRow][nextCol], &edgeTo[nextRow][nextCol], 1, 0);
		shmem_quiet();
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

//Find the indexes of the vertical seams to be carved out
int * identifySeams( int width, int height, int start_col, int stop_col, int npes){
	distTo = initializeDistTo(height, width);
	edgeTo = initializeEdgeTo(height, width);

	//Initialize distTo to maximum values

	 for (int row = 0; row < height - 1; row++) {
            for (int col = start_col; col < stop_col; col++) {
                relax(row, col, edgeTo, distTo, width);
            }
		//All PEs wait for one another
		shmem_barrier_all();

		/*PE 0 broadcasts its new edges and distance values
			so all PEs have up-to-date values*/
		 if (npes > 1){
                          shmem_broadcast64(&distTo[row][0], &distTo[row][0], width, 0, 0, 0, npes, pSync);
			  shmem_broadcast64(&edgeTo[row][0], &edgeTo[row][0], width, 0, 0, 0, npes, pSync);
		}
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
	//static long pSync[SHMEM_BCAST_SYNC_SIZE];
	 for (int i = 0; i < SHMEM_BCAST_SYNC_SIZE; i++)
                pSync[i] = SHMEM_SYNC_VALUE;

	shmem_init();
	int me = shmem_my_pe();
	int npes = shmem_n_pes();
	
	double total_begin = timestamp();
	
	char * original_img = argv[1]; 
	char* orientation = argv[2];
	
	/*read the image and get the dimension*/

	pngwrt.readfromfile(original_img);
	width = pngwrt.getwidth();
	height = pngwrt.getheight();
	

        if (me == 1)
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

		//initialize energy array
		energyArray = initializeEnergyArray(height, width);
		
		/*Divide the work among the available PEs*/
		int col_per_pe = (width + npes - 1)/npes;
		int start_col = me * col_per_pe;
		int stop_col = (me + 1) * col_per_pe;

		/*Ensure that the last pe does not exceed the last row*/
		if (me == npes - 1)
			stop_col  = width;
		
		//Fill the energy matrix with the energy values of each pixel
		generateEnergyMatrix(height, start_col, stop_col, npes, orientation);
		shmem_barrier_all();

       		 if (npes > 1)
               		  shmem_broadcast64(energyArray, energyArray, width*height, 0, 0, 0, npes, pSync);

	
		if (me  == 0)
			cout<<"Removing vertical seams"<<endl;
		
		//Find the vertical seam in the image
		identifySeams(width, height, start_col, stop_col, npes);
		int* v_seams =  backTrack(edgeTo, distTo, height, width);

		/*PE 0 will remove the identified seams
		  there is no need to parallelize this process
		  no significant speedup is attained with its parallelizations*/

		if (me == 0){
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
		
		/*Compute energy values*/
		generateEnergyMatrix(width, start_col, stop_col, npes, orientation);
		shmem_barrier_all();

                 if (npes > 1)
                          shmem_broadcast64(energyArray, energyArray, width*height, 0, 0, 0, npes, pSync);

		if (me == 0)
			cout<<"Removing horizontal seams"<<endl;
		
		/*Find the horizontal seams in the image*/
		identifySeams(height, width, start_col, stop_col, npes);
		int* h_seams =  backTrack(edgeTo, distTo, width, height);
		
		/*Revome the identified horizontal seams*/
		if (me == 0){
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
	//if (me == 0) {
		double total_end = timestamp();
		printf("%s%5.2f\n","Total Processing Time: ", (total_end-total_begin));
	}

	shmem_finalize();
	return 0;
	
}

