#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int** initialPadding(int** Image, int height, int width);
int** array2mat(int *array, int row, int col);
int* mat2Array(int** mat, int row, int col);
void print_matrix(int** mat, int row, int col);
void print_array(int* array,int size);
int** duplicateRow(int** Matrix, int p, int numRows, int numCols);
int convolve(int filter[3][3], int window[3][3]);
int** applyFilter(int** subImage, int filter[3][3], int numRows, int numCols);
int** loadImage(const char* filename);
int saveImage(int** matrix, int row, int col);

int SOBEL_H[3][3] = {
	{1,2,1},
	{0,0,0},
	{-1,-2,-1}
};

int SOBEL_V[3][3] = {
	{-1,0,1},
	{-2,0,2},
	{-1,0,1}
};


int PREWITT_H[3][3] = {
	{1,1,1},
	{0,0,0},
	{-1,-1,-1}
};

int PREWITT_V[3][3] = {
	{-1,0,1},
	{-1,0,1},
	{-1,0,1}
};

int LAPLACIAN[3][3] = {
	{-1,-1,-1},
	{-1,8,-1},
	{-1,-1,-1}
};


int main(int argc, char* argv[]){
	
	int my_rank, p, i, j;
	int imgRows = 512;
	int paddedImgRows;
	int imgCols = 512;
	int paddedImgCols;
	int** Image;
	int** paddedImage;
	int** subImage;
	int** filtImage;
	int** finalImage;
	int* Image1D;
	int* subImage1D;
	int* paddedImage1D;
	int* filtImage1D;
	int* finalImage1D;
	
	int subImageSize;
	
	int filterCode;
	int filter[3][3];
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	if (my_rank == 0){
	Image = loadImage("InputMatrix.txt");	
	finalImage1D = malloc(imgRows * imgCols * sizeof(int));
	/* 
		filter types: 1 -> 6 (the same order above)
	*/
	filterCode = 1; // SOBEL_H

	//duplicate the rows the pad the image (to perform overlapped blocking):
	paddedImage = duplicateRow(Image, p, imgRows, imgCols);
	paddedImgRows = imgRows + ((p-1) * 2);
	paddedImgCols = imgCols;
	paddedImage = initialPadding(paddedImage, paddedImgRows, paddedImgCols);
	paddedImgRows += 2;
	paddedImgCols += 2;
	
	// Transforming 2D Image to 1D array;
	paddedImage1D = mat2Array(paddedImage, paddedImgRows, paddedImgCols);
	}
	
	MPI_Bcast(&imgRows, 1 , MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&imgCols, 1 , MPI_INT, 0, MPI_COMM_WORLD);
	paddedImgRows = imgRows + (2*p);
	paddedImgCols = imgCols + 2;
	
	MPI_Bcast(&filterCode, 1 , MPI_INT, 0, MPI_COMM_WORLD);
	
	subImageSize = paddedImgRows * paddedImgCols / p;
	subImage1D = malloc(subImageSize * sizeof(int));
	MPI_Scatter (&(paddedImage1D[0]), subImageSize, MPI_INT, subImage1D, subImageSize, MPI_INT, 0, MPI_COMM_WORLD);
	//printf("%d \n",subImage1D[1]);
	subImage = array2mat(subImage1D, paddedImgRows/p, paddedImgCols);
	filtImage = applyFilter(subImage, LAPLACIAN,  paddedImgRows/p, paddedImgCols);
	filtImage1D = mat2Array(filtImage, imgRows/p, imgCols);
	MPI_Gather(&(filtImage1D[0]), imgRows * imgCols / p, MPI_INT, finalImage1D, imgRows * imgCols / p, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (my_rank == 0){
		finalImage = array2mat(finalImage1D, imgRows, imgCols);
		if(!saveImage(finalImage, imgRows, imgCols))
		{
			printf("Oh how does it feel the success : GREAT\n");	
		}else
		{
			printf("Oh Shit\n");
		}
	}
	MPI_Finalize();
}


 /**
 	The initial padding function
 **/
int** initialPadding(int** Image, int height, int width)
{
	int** paddMatrix = malloc(sizeof(int*) * (height + 2));
	int i, j;
	for(i = 0; i<(height + 2); i++)
	{
		paddMatrix[i] = malloc(sizeof(int) * (width + 2));
	}
	for(i = 0; i<(height + 2); i++)
	{
		for (j = 0; j < (width + 2); j++)
		{
			if((i == 0 || j==0) || (i==(height+1)))
			{
				paddMatrix[i][j] = 0;		
			}else{
				paddMatrix[i][j] = Image[i-1][j-1];
			}
		}
	}
	return paddMatrix;
}

/*
	The function to convert array to matrix
*/

int** array2mat(int *array, int row, int col)
{
	int** mat = malloc(sizeof(int*) * row);
	int i, j = 0, k = 0;
	for (i = 0; i < row; i++)
	{
		mat[i] = malloc(sizeof(int) * col);
	}
	for (i = 0; i < row*col; i++)
	{
		if (i % col == 0 && i!=0)
		{
			j++;
		}
		mat[j][i%col] = array[i];
	}
	return mat;
}

/*
	The function to convert matrix to array
*/
int* mat2Array(int** mat, int row, int col)
{
	int* array = malloc(sizeof(int) * row * col);
	int i, j, k=0;
	for(i=0; i<row; i++)
	{
		for(j=0; j<col; j++)
		{
			array[k++] = mat[i][j];
		}
	}
	return array;
}

/*
	The function to print matrix in a nice form
*/
void print_matrix(int** mat, int row, int col)
{
	int i, j;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			printf("%d ", mat[i][j]);
		}
		printf("\n");
	}

}

/*
	The function to print arrayr in a nice form
*/
void print_array(int* array,int size)
{
	int i;
	for (i = 0; i < size; i++)
	{
		printf("%d ", array[i]);
	}
}


int** duplicateRow(int** Matrix, int p, int numRows, int numCols){
    int newNumRows = numRows + (2 * (p-1));

    int i, j;
    int** newMatrix = malloc(newNumRows * sizeof(int*));
    for (i=0; i<newNumRows; i++)
        newMatrix[i] = malloc(numCols * sizeof(int));

    int step = numRows / p;
    int startIndex = step - 1;
    int nextRep = startIndex;
    int numDuplication = p - 1;
    int destIndex = 0;
    int sourceIndex = 0;
    while(sourceIndex < numRows){
        if (sourceIndex == nextRep && sourceIndex != numRows-1){
            for (j=0; j<numCols; j++){
                newMatrix[destIndex][j] = Matrix[sourceIndex][j];
                newMatrix[destIndex+1][j] = Matrix[sourceIndex+1][j];
                newMatrix[destIndex+2][j] = Matrix[sourceIndex][j];
                newMatrix[destIndex+3][j] = Matrix[sourceIndex+1][j];
            }
            destIndex += 4;
            sourceIndex += 2;
            nextRep = nextRep + step;
        }
        else {

            for (j=0; j<numCols; j++){
                newMatrix[destIndex][j] = Matrix[sourceIndex][j];
            }
            destIndex += 1;
            sourceIndex +=1;
        }
    }

    return newMatrix;

}

int** applyFilter(int** subImage, int filter[3][3], int numRows, int numCols){

    int i, j, k, l;
    int** newMatrix = malloc((numRows-2) * sizeof(int*));
    for (i=0; i<(numRows-2); i++)
        newMatrix[i] = malloc((numCols-2) * sizeof(int));
    int window[3][3];

    for (i=1; i<(numRows-1); i++){
        for (j=1; j<(numCols-1); j++){
            for (k=-1; k<2; k++){
                for (l=-1; l<2; l++){
                    window[k+1][l+1] = subImage[i+k][j+l];
                }
            }
            newMatrix[i-1][j-1] = convolve(filter, window);
        }
    }
    return newMatrix;


}

int convolve(int filter[3][3], int window[3][3]){
    int i, j, result = 0;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            result += filter[i][j] * window[i][j];
        }
    }
    return abs(result);
}
/*
	@Description: The function to read an image into a 2D-matrix
	@param : filename (char *)
	@return : 2D- matrix of the grayscale data
*/
int** loadImage(const char* filename)
{
	FILE* fp;
	fp = fopen(filename, "r");
	int** matrix = malloc(sizeof(int*)*512);
	int i, j;
	if(fp != NULL)
	{
		for (i = 0; i < 512; i++)
		{
			for(j = 0; j < 512; j++){
				matrix[i] = malloc(sizeof(int) * 512);
			}
		}
		for (i = 0; i < 512; i++)
		{
			for(j = 0; j < 512; j++){
				fscanf(fp, "%d", &matrix[i][j]);
			}
		}
		//print_matrix(matrix, 512, 512);
	}else
	{
		printf("Failed to read file\n");
		return NULL;
	}
	fclose(fp);	
	return matrix;
}


int saveImage(int** matrix, int row, int col)
{
	FILE* fp;
	fp = fopen("OutputMatrix.txt", "w");
	if(fp != NULL)
	{
		int i,j;
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				fprintf(fp, "%d ", matrix[i][j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		return 0;
	}
	return -1;
}

/*int main(int argc, char *argv[])
{
	int** mt = loadImage("imageMat.txt");
	//print_matrix(mt, 512, 512);
	if(!saveImage(mt, 512, 512))
	{
		printf("OOHHHH Love success\n");
	}
	return 0;
}*/