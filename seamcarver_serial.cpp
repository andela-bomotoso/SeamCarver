#include <pngwriter.h>
#include <lqr.h>
#include <getopt.h>
#include "liquidrescale.h"
#include <math.h>
#include <limits>
#include <sys/time.h>
using namespace std;
pngwriter pngwrt(1,1,0,"out.png");
int BASE_ENERGY = 1000;
int width = 0;
int height = 0;
gfloat rigidity = 0;
gint max_step = 1;
double** energyArray;
guchar* seams;

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


double computeEnergy(int x, int y, guchar* buffer){

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
	return sqrt(valueSum);
}

// energy of pixel at column x and row y
double getEnergy(int x, int y, guchar* buffer) {

        if (x == 0 || y == 0 || (x == width - 1) || (y == height- 1))
            return BASE_ENERGY;
        else
            return computeEnergy(x, y, buffer);
    }           


/*Declare a relax function to optimize the computation of a 
	shortest path energy values*/

void relax(int row, int col, int** edgeTo, double** distTo) {
        int nextRow = row + 1;
        for (int i = -1; i <= 1; i++) {
            int nextCol = col + i;
            if (nextCol < 0 || nextCol >= width)
                continue;
            if (distTo[nextRow][nextCol] >= distTo[row][col] + energyArray[nextRow][nextCol]) {
                distTo[nextRow][nextCol] = distTo[row][col] + energyArray[nextRow][nextCol];
                edgeTo[nextRow][nextCol] = i;
            }
        }
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

//TO BE PARALLELIZED
//Find the indexes of the vertical seams to be carved out
int * findVerticalSeams( guchar* buffer, int width, int height){
	//Declare an array to hold the indexes of the vertical seam
	 int * verticalSeams = new int[height];
	//Declare a dynamic 2D array to hold the energy values for all pixels
	energyArray = new double*[height];
	for (int i = 0; i < height; i++)
		energyArray[i] = new double[width];

	/*Declare an array to hold the distance to reach a pixel 
		using the SSSP approach*/
	double** distTo = new double*[height];
	for (int i = 0; i < height; i++)
		distTo[i] = new double[width];

	//Declare an array to hold the paths taken to reach a pixel
	int** edgeTo = new int*[height];
	for (int i = 0; i < height; i++)
		edgeTo[i] = new int[width];


	//Fill up the energyArray with energy values of each image pixel
	for (int i = 1; i < height; i++)	{
		for(int j = 0; j < width; j++)	{
			energyArray[i][j] = getEnergy(j,i, buffer );
		}
	}

	//Initialize distTo to maximum values
	
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                if (row == 0)
                    distTo[row][col] = BASE_ENERGY;
                else
                    distTo[row][col] = std::numeric_limits<double>::infinity();
            }
       }
	 for (int row = 0; row < height - 1; row++) {
            for (int col = 0; col < width; col++) {
                relax(row, col, edgeTo, distTo);
            }
        }
	// Backtrack from the last row to get a shortest path
        int minCol = 0;
        double minDist = std::numeric_limits<double>::infinity();
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

//Find the indexes of the horizontal seams to be carved out
int* findHorizontalSeams(guchar* buffer){
	//Declare an array to hold the indexes of the horizontal seams
	int temp = height;
	int new_height = width;
	int new_width = temp;
        guchar* transBuffer = transposeRGBuffer(buffer,width, height );
        int* horizontalSeams = findVerticalSeams(transBuffer, new_width, new_height);
        return horizontalSeams;
}
//TO BE PARALLELIZED
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
//Carve out the horizontal seams
guchar* carveHorizontally(int* horizontal_seams, guchar* buffer, int width, int height){
       guchar*  carved_imageH = carveVertically(horizontal_seams, buffer, width, height);
	return carved_imageH;
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

        /* plot (pngwriter's coordinates start from 1,1) */
	if (red != 0 && green != 0 && blue !=0)
       		 pngwrt->plot(x + 1, y + 1, red, green, blue);
	else{
		       pngwrt->plot(x + 1, y + 1, 1.0, 0.0, 0.0);
			
}
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
	char * original_img = argv[1]; 
	char* orientation = argv[2];
	double begin, end;
	begin = timestamp();
	pngwrt.readfromfile(original_img);
	width = pngwrt.getwidth();
	height = pngwrt.getheight();
	
	cout<<"Width: "<<width<<" Height: "<<height<<endl;
	guchar* buffer = rgb_buffer_from_image(&pngwrt);
	LqrCarver *carver;
	LqrCarver *carved_seams;
	//Check the orientation to determine how to carve
	if(orientation[0] == 'v'){
		cout<<"Removing vertical seams"<<endl;
		int* v_seams = findVerticalSeams(buffer, width, height);
		guchar* carved_imageV = carveVertically(v_seams,buffer, width, height);
		carver = lqr_carver_new(carved_imageV, width, height, 3);
		carved_seams = lqr_carver_new(seams, width, height, 3);
	}
	else{
		cout<<"Removing horizontal seams"<<endl;
		//Transpose buffer
		guchar* transBuffer = transposeRGBuffer(buffer, width, height);
		int* h_seams = findHorizontalSeams(buffer);
		guchar* carved_imageH = carveVertically(h_seams, transBuffer, height, width);
		carver = lqr_carver_new(transposeRGBuffer(carved_imageH, height,width), width, height, 3);
		carved_seams = lqr_carver_new(transposeRGBuffer(seams, height,width), width, height, 3);
	}

	//Create a Carver object with the carved image buffer
   	TRAP(lqr_carver_init(carver, max_step, rigidity));
	//write_carver_to_image(carver, &pngwrt, orientation);
	printSeams(carved_seams, &pngwrt);
	lqr_carver_destroy(carver);
	pngwrt.close();
	end = timestamp();
	printf("%s%5.2f\n","TOTAL TIME: ", (end-begin));

	return 0;
}

	


