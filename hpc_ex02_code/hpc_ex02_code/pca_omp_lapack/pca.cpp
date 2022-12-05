// File       : pca.cpp
// Description: Principal component analysis applied to image compression
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <omp.h>
#include <zlib.h>

// interface for LAPACK routines.
//#include <>


///////////////////////////////////////////////////////////////////////////////
// helpers
double *read_gzfile(char *filename, int frows, int fcols, int rows, int cols)
{
    double *A, *buf;
    gzFile fp;
    int i;

    A = new (std::nothrow) double[rows*cols];
    assert(A != NULL);

    buf = new (std::nothrow) double[fcols];
    assert(buf != NULL);

    fp = gzopen(filename, "rb");
    if (fp == NULL) {
        std::cout << "Input file not available!\n";
        exit(1);
    }

    for (i = 0; i < rows; i++) {
        gzread(fp, buf, fcols*sizeof(double));
        memcpy(&A[i*cols], buf, cols*sizeof(double));
    }
    gzclose(fp);

    delete[] buf;

    return A;
}

void write_ascii(const char* const filename, const double* const data, const int rows, const int cols)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        std::cout << "Failed to create output file\n";
        exit(1);
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(fp, "%.4lf ", data[i*cols+j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
///////////////////////////////////////////////////////////////////////////////

// Caclulate the mean of a feature column.
double calc_column_mean(double *column, int column_length)
{
    // Initialise sum variable.
    double sum = 0;

    // Sum all values.
    for (int i = 0; i < column_length; i++)
    {
        sum += column[i]; 
    }

    // Return mean.
    return sum/column_length;
}

double calc_column_std(double *column, double mean, int column_length)
{
        // Initialise sum variable.
    double sum = 0;

    // Sum all values' deviation from mean squared.
    for (int i = 0; i < column_length; i++)
    {
        sum += pow(column[i] - mean, 2); 
    }

    // Return standard deviation.
    return sum/column_length;
}

///////////////////////////////////////////////////////////////////////////////
//
// elvis.bin.gz:   469x700
// cyclone.bin.gz: 4096x4096
// earth.bin.gz:   9500x9500
//
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    // input parameters (default)
    int image_rows = 469, image_columns = 700;                        // image size (rows, columns)
    int npc = 50;                                // number of principal components
    char *inp_filename = (char *)"../pca_data/elvis.bin.gz"; // input filename (compressed binary)
    char *out_filename = (char *)"elvis.50.bin"; // output filename (text)

    // parse input parameters
    if ((argc != 1) && (argc != 11)) {
        std::cout << "Usage: " << argv[0] << " -m <rows> -n <cols> -npc <number of principal components> -if <input filename> -of <output filename>\n";
        exit(1);
    }

    for(int i = 1; i < argc; i++ ) {
        if( strcmp( argv[i], "-m" ) == 0 ) {
            image_rows = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-n" ) == 0 ) {
            image_columns = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-npc" ) == 0 ) {
            npc = atoi(argv[i+1]);
            i++;
        }
        if( strcmp( argv[i], "-if" ) == 0 ) {
            inp_filename = argv[i+1];
            i++;
        }
        if( strcmp( argv[i], "-of" ) == 0 ) {
            out_filename = argv[i+1];
            i++;
        }
    }

    if (npc > image_columns) npc = image_columns;
    double t_elapsed;

    ///////////////////////////////////////////////////////////////////////////
    // Read image data.  The image dimension is image_rows x image_columns.  The returned pointer
    // points to the data in row-major order.  That is, if (i,j) corresponds to
    // to the row and column index, respectively, you access the data with
    // pixel_{i,j} = I[i*image_columns + j], where 0 <= i < image_rows and 0 <= j < image_columns.
    ///////////////////////////////////////////////////////////////////////////
    double *I = read_gzfile(inp_filename, image_rows, image_columns, image_rows, image_columns);

    // A = transpose(I), so image features (columns) are stored in rows.  More
    // efficient to work with the data in this layout.
    double *A = new (std::nothrow) double[image_columns*image_rows];
    for (int i = 0; i < image_columns; i++) {
        for (int j = 0; j < image_rows; j++) {
            A[i*image_rows + j] = I[j*image_columns+i];
        }
    }
    delete[] I;

    ///////////////////////////////////////////////////////////////////////////
    // TODO: Implement your PCA algorithm here
    // 1. Compute mean and standard deviation of your image features
    // 2. Normalize the data
    // 3. Build the covariance matrix
    // 4. Compute the eigenvalues and eigenvectors of the covariance matrix.
    //    Use LAPACK here.
    // 5. Compute the principal components and report the compression ratio
    // 6. Reconstruct the image from the compressed data and dump the image in
    //    ascii.
    ///////////////////////////////////////////////////////////////////////////
    double start_t = omp_get_wtime();

    ///////////////////////////////////////////////////////////////////////////
    // TODO: 1.
    t_elapsed = -omp_get_wtime();

    double *AMean = new (std::nothrow) double[image_columns];
    double *AStd  = new (std::nothrow) double[image_columns];
    assert(AMean != NULL);
    assert(AStd  != NULL);

    for (int i = 0; i < image_columns; i++)
    {
        // TODO: Make these two an openMP task.
        // Calculate mean of each column.
        AMean[i] = calc_column_mean(&A[i*image_rows], image_rows); 

        AStd[i] = calc_column_std(&A[i*image_rows], AMean[i], image_rows);
    }
    t_elapsed += omp_get_wtime();
    std::cout << "MEAN/STD TIME=" << t_elapsed << " seconds\n";
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // TODO: 2.
    t_elapsed = -omp_get_wtime();

    for (int i = 0; i < image_columns; i++)
    {
        for (int j = 0; j < image_rows; j++)
        {
            // TODO: normalize data here
            // subtract mean and divide by standard deviation.
        }
    }
    t_elapsed += omp_get_wtime();
    std::cout << "NORMAL. TIME=" << t_elapsed << " seconds\n";
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // TODO: 3.
    t_elapsed = -omp_get_wtime();
    double *C = new (std::nothrow) double[image_columns*image_columns];
    assert(C!=NULL);

    // TODO: Compute covariance matrix here

    t_elapsed += omp_get_wtime();
    std::cout << "C-MATRIX TIME=" << t_elapsed << " seconds\n";
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // TODO 4. LAPACK
    t_elapsed = -omp_get_wtime();

    // see also for the interface to dsyev_():
    // http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
    char jobz = '?'; // TODO: compute both, eigenvalues and orthonormal eigenvectors
    char uplo = '?'; // TODO: how did you compute the (symmetric) covariance matrix?
    int info, lwork;

    double *W = new (std::nothrow) double[image_columns]; // eigenvalues
    assert(W != NULL);

    double *work = new (std::nothrow) double[2];
    assert(work != NULL);

    // first call to dsyev_() with lwork = -1 to determine the optimal
    // workspace (cheap call)
    lwork = -1;

    // TODO: call dsyev here

    lwork = (int)work[0];
    delete[] work;

    // allocate optimal workspace
    work = new (std::nothrow) double[lwork];
    assert(work != NULL);

    // second call to dsyev_(), eigenvalues and eigenvectors are computed here
    // TODO: call dsyev here

    t_elapsed += omp_get_wtime();
    std::cout << "DSYEV TIME=" << t_elapsed << " seconds\n";
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // TODO: 5.
    t_elapsed = -omp_get_wtime();
    double *PCReduced = new (std::nothrow) double[image_rows*npc];
    assert(PCReduced != NULL);

    for (int i = 0; i < image_rows; i++)
    {
        for (int j = 0; j < npc; j++)
        {
            // TODO: compute the principal components
        }
    }

    // TODO: Report the compression ratio

    t_elapsed += omp_get_wtime();
    std::cout << "PCREDUCED TIME=" << t_elapsed << " seconds\n";
    ///////////////////////////////////////////////////////////////////////////

    double end_t = omp_get_wtime();
    std::cout << "OVERALL TIME=" << end_t - start_t << " seconds\n";

    ///////////////////////////////////////////////////////////////////////////
    // TODO: 6
    double *Z = new (std::nothrow) double[image_rows*image_columns]; // memory for reconstructed image
    assert(Z != NULL);

    for (int i = 0; i < image_rows; i++)
    {
        for (int j = 0; j < image_columns; j++)
        {
            // TODO: Reconstruct image here.  Don't forget to denormalize.  The
            // dimension of the reconstructed image is m x n (rows x columns).
            // Z[i*image_columns + j] = ...
        }
    }


    // Write the reconstructed image in ascii format.  You can view the image
    // in Matlab with the show_image.m script.
    write_ascii(out_filename, Z, image_rows, image_columns);
    ///////////////////////////////////////////////////////////////////////////

    // cleanup
    delete[] work;
    delete[] W;
    delete[] C;
    delete[] Z;
    delete[] PCReduced;
    delete[] A;
    delete[] AMean;
    delete[] AStd;

    return 0;
}
