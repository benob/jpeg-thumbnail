#include <stdio.h>
#include <jpeglib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* we will be using this uninitialized pointer later to store raw, uncompressd image */
unsigned char *raw_image = NULL;
unsigned char *raw_image_in = NULL;

/* dimensions of the image we want to write */
int width = 1000;
int height = 1000;
int max_dimension = 256;
int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

//#define USE_NEAREST
#define USE_BILINEAR_INT
//#define USE_BILINEAR_DOUBLE

/**
 * read_jpeg_file Reads from a jpeg file on disk specified by filename and saves into the 
 * raw_image buffer in an uncompressed format.
 * 
 * \returns positive integer if successful, -1 otherwise
 * \param *filename char string specifying the file name to read from
 *
 */

int read_jpeg_file( char *filename )
{
	/* these are standard libjpeg structures for reading(decompression) */
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	/* libjpeg data structure for storing one row, that is, scanline of an image */
	JSAMPROW row_pointer[1];
	
	FILE *infile = fopen( filename, "rb" );
	unsigned long location = 0;
	int i = 0;
	
	if ( !infile )
	{
		printf("Error opening jpeg file %s\n!", filename );
		return -1;
	}
	/* here we set up the standard libjpeg error handler */
	cinfo.err = jpeg_std_error( &jerr );
	/* setup decompression process and source, then read JPEG header */
	jpeg_create_decompress( &cinfo );
	/* this makes the library read from infile */
	jpeg_stdio_src( &cinfo, infile );
	/* reading the image header which contains image information */
	jpeg_read_header( &cinfo, TRUE );
    fprintf(stderr, "size: %d %d\n", cinfo.image_width, cinfo.image_height);
    cinfo.scale_num = max_dimension;
    if(cinfo.image_width < cinfo.image_height) {
        cinfo.scale_denom = cinfo.image_height;
    } else {
        cinfo.scale_denom = cinfo.image_width;
    }
    cinfo.dct_method = JDCT_FASTEST;
	/* Uncomment the following to output image information, if needed. */
	/*--
	printf( "JPEG File Information: \n" );
	printf( "Image width and height: %d pixels and %d pixels.\n", cinfo.image_width, cinfo.image_height );
	printf( "Color components per pixel: %d.\n", cinfo.num_components );
	printf( "Color space: %d.\n", cinfo.jpeg_color_space );
	--*/
	/* Start decompression jpeg here */
	jpeg_start_decompress( &cinfo );

	/* allocate memory to hold the uncompressed image */
	raw_image_in = (unsigned char*)malloc( cinfo.output_width*cinfo.output_height*cinfo.num_components );
	/* now actually read the jpeg into the raw buffer */
	row_pointer[0] = (unsigned char *)malloc( cinfo.output_width*cinfo.num_components );
	/* read one scan line at a time */
    fprintf(stderr, "presize: %d %d\n", cinfo.output_width, cinfo.output_height);
	while( cinfo.output_scanline < cinfo.output_height )
	{
		jpeg_read_scanlines( &cinfo, row_pointer, 1 );
		for( i=0; i<cinfo.output_width*cinfo.num_components;i++) 
			raw_image_in[location++] = row_pointer[0][i];
	}

    width = max_dimension;
    height = max_dimension * cinfo.image_height / cinfo.image_width;
    if(cinfo.image_height > cinfo.image_width) {
        width = max_dimension * cinfo.image_width / cinfo.image_height;
        height = max_dimension;
    }
    raw_image = (unsigned char*)malloc( width * height * cinfo.num_components);
#ifdef USE_BILINEAR_DOUBLE
    double dx = (double) (cinfo.output_width - 1) / (double) width;
    double dy = (double) (cinfo.output_height - 1) / (double) height;
    double line = 0;
    //fprintf(stderr, "resize: %d %d %d %f %f\n", cinfo.num_components, width, height, dx, dy);
    int j = 0, k = 0;
    // bilinear
    for(i = 0; i < height; i++) {
        double y = (dy * i);
        double factor2 = (y - floor(y)) / dy;
        int line = (int) y;
        for(j = 0; j < width; j++) {
            double x = (dx * j);
            double factor1 = (x - floor(x)) / dx;
            int column = (int) x;
#define get_pixel(a, b) (double) raw_image_in[((a) + (b) * cinfo.output_width) * cinfo.num_components + k]
            for(k = 0; k < cinfo.num_components; k++) {
                double result = 
                    (1 - factor1) * (1 - factor2) * get_pixel(column, line) + factor1 * (1 - factor2) * get_pixel(column + 1, line)
                  + (1 - factor1) * factor2 * get_pixel(column + 1, line + 1) + factor1 * factor2 * get_pixel(column, line + 1)
                    ;
                raw_image[(j + i * width) * cinfo.num_components + k] = (int) result;
            }
        }
    }
#endif
#ifdef USE_BILINEAR_INT
    unsigned int dx = ((cinfo.output_width - 1) << 16) / width;
    unsigned int dy = ((cinfo.output_height - 1) << 16) / height;
    unsigned int line = 0;
    unsigned int j = 0, k = 0;
    for(i = 0; i < height; i++) {
        unsigned int y = (dy * i);
        unsigned int factor2 = (y & 0xffff0000) / dy;
        unsigned int line = y >> 16;
        for(j = 0; j < width; j++) {
            unsigned int x = (dx * j);
            unsigned int factor1 = (x & 0xffff0000) / dx;
            unsigned int column = x >> 16;
            unsigned int factor_1_f1_1_f2 = (((0xffff - factor1) * (0xffff - factor2)) >> 16);
            unsigned int factor_f1_1_f2 = (((factor1) * (0xffff - factor2)) >> 16);
            unsigned int factor_1_f1_f2 = (((0xffff - factor1) * (factor2)) >> 16);
            unsigned int factor_f1_f2 = (((factor1) * (factor2)) >> 16);
#define get_pixel(a, b) (unsigned int) raw_image_in[((a) + (b) * cinfo.output_width) * cinfo.num_components + k]
            for(k = 0; k < cinfo.num_components; k++) {
                unsigned int result = (factor_1_f1_1_f2 * get_pixel(column, line)) >> 16;
                result += (factor_f1_1_f2 * get_pixel(column + 1, line)) >> 16;
                result += (factor_1_f1_f2 * get_pixel(column + 1, line + 1)) >> 16;
                result += (factor_f1_f2 * get_pixel(column, line + 1)) >> 16;
                /*if(result < 0) result = 0; // clip
                if(result > 255) result = 255;*/
                raw_image[(j + i * width) * cinfo.num_components + k] = result & 0xff;
                //fprintf(stdout, "%d %d %d %d %d\n", line, column, factor1, factor2, result);
            }
        }
    }
#endif
#ifdef USE_NEAREST
    double dx = (double) (cinfo.output_width - 1) / (double) width;
    double dy = (double) (cinfo.output_height - 1) / (double) height;
    double line = 0;
    int j = 0, k = 0;
    // nearest neighbor
    for(line = 0, i = 0; i < height && line < cinfo.output_height; i++, line += dy) {
        int int_line = (int) line;
        double column;
        for(column = 0, j = 0; j < width && column < cinfo.output_width; j++, column += dx) {
            int int_column = (int) column;
            for(k = 0; k < cinfo.num_components; k++) {
                raw_image[(j + i * width) * cinfo.num_components + k] = raw_image_in[(int_column + int_line * cinfo.output_width) * cinfo.num_components + k];
            }
        }
    }
#endif
	/* wrap up decompression, destroy objects, free pointers and close open files */
	jpeg_finish_decompress( &cinfo );
	jpeg_destroy_decompress( &cinfo );
	free( row_pointer[0] );
    //free(raw_image_in);
	fclose( infile );
	/* yup, we succeeded! */
	return 1;
}

/**
 * write_jpeg_file Writes the raw image data stored in the raw_image buffer
 * to a jpeg image with default compression and smoothing options in the file
 * specified by *filename.
 *
 * \returns positive integer if successful, -1 otherwise
 * \param *filename char string specifying the file name to save to
 *
 */
int write_jpeg_file( char *filename )
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	
	/* this is a pointer to one row of image data */
	JSAMPROW row_pointer[1];
	FILE *outfile = fopen( filename, "wb" );
	
	if ( !outfile )
	{
		printf("Error opening output jpeg file %s\n!", filename );
		return -1;
	}
	cinfo.err = jpeg_std_error( &jerr );
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, outfile);

	/* Setting the parameters of the output file here */
	cinfo.image_width = width;	
	cinfo.image_height = height;
	cinfo.input_components = bytes_per_pixel;
	cinfo.in_color_space = color_space;
    //jpeg_set_quality(&cinfo, 80, 0);
    cinfo.dct_method = JDCT_FASTEST;
    /* default compression parameters, we shouldn't be worried about these */
	jpeg_set_defaults( &cinfo );
	/* Now do the compression .. */
	jpeg_start_compress( &cinfo, TRUE );
	/* like reading a file, this time write one row at a time */
	while( cinfo.next_scanline < cinfo.image_height )
	{
		row_pointer[0] = &raw_image[ cinfo.next_scanline * cinfo.image_width *  cinfo.input_components];
		jpeg_write_scanlines( &cinfo, row_pointer, 1 );
	}
	/* similar to read file, clean up after we're done compressing */
	jpeg_finish_compress( &cinfo );
	jpeg_destroy_compress( &cinfo );
	fclose( outfile );
	/* success code is 1! */
	return 1;
}

int main(int argc, char** argv)
{
    if(argc != 4) {
        fprintf(stdout, "USAGE: %s <infile> <outfile> <maxdim>\n", argv[0]);
        exit(1);
    }

	char *infilename = argv[1], *outfilename = argv[2];
    max_dimension = strtol(argv[3], NULL, 10);

    if(max_dimension <= 0) {
        fprintf(stderr, "ERROR: invalid dimension %s\n", argv[3]);
    }

	/* Try opening a jpeg*/
	if( read_jpeg_file( infilename ) > 0 ) 
	{
		/* then copy it to another file */
		if( write_jpeg_file( outfilename ) < 0 ) return -1;
	}
	else return -1;
	return 0;
}


