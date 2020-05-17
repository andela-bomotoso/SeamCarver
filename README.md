# Parallelilizing Seam Carving Using OpenSHMEM

We  parallelized Seam Carving with pthreads, MPI and OpenSHMEM in this project

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Dependencies

The Liquid Rescale Library: Follow the installation instructions on the liblqr gitHub page  https://github.com/carlobaldassi/liblqr
PNGwriter: Follow the installation instructions on the PGNwriter gitHub page https://github.com/pngwriter/pngwriter


## Running the test
Use the make file for the execution. 
A sample run takes 2 arguments; the png file and the orientation. Use v for the vertical orientation and h for the horizontal orientation. The out.png file has the processed image with the seams path. 
e.g ./seamcarver_serial burro.png v: will generate a vertical seam paths in out.png for the image burro.png
    ./seamcarver_serial burro.png h: will generate a horizontal seam paths in out.png for the image burro.png



## Acknowledgments

* Dr. Ferrol Aderholdt
* Dr. Suk Seo
* Dr. Chrisila Pettey

