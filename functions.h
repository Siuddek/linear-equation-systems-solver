#pragma once
#define SIZE 996 //9cd - 996

//A----------------------------------------------------

inline double **initialization() {
	//initializing matrix as 2d array in O(n)
	double **matrix = new double*[SIZE];
	for (int i = 0; i < SIZE; i++) {
		matrix[i] = new double[SIZE];
	}

	//filling matrix with zeros O(n^2)
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			matrix[i][j] = 0;
		}
	}
	return matrix;
}

inline void fillMatrix(double **matrix, double a1, double a2, double a3) {
	//creating band matrix in O(n^2)
	for (int i = 0; i < SIZE - 2; i++) {
		for (int j = 0; j < SIZE - 2; j++) {
			if (i == j) {//we are on diagonal
				matrix[i][j] = a1;
				//filing parallel diagonals
				matrix[i + 1][j] = matrix[i][j + 1] = a2;
				matrix[i + 2][j] = matrix[i][j + 2] = a3;
			}
		}
	}

	matrix[SIZE - 2][SIZE - 2] = matrix[SIZE - 1][SIZE - 1] = a1;
	matrix[SIZE - 1][SIZE - 2] = matrix[SIZE - 2][SIZE - 1] = a2;
}

inline void fillVectorB(double *vector) {
	//initializing vector according to pattern vec[i] = sin(n*(f+1))
	for (int i = 0; i < SIZE; i++){
		vector[i] = sin(2 * i);
	}
}

//B-------------------------------------------------\

inline void fillRevNegDiagonal(double **diag, double a1) {
	double tmp = -(1 / a1);
	for (int i = 0; i < SIZE; i++){
		diag[i][i] = tmp;
	}
}

inline void fillRevDiagonal(double **diag, double a1) {
	double tmp = 1 / a1;
	for (int i = 0; i < SIZE; i++) {
		diag[i][i] = tmp;
	}
}

inline void fillLowerUpper(double **matrix, double a1, double a2) {
	//creating band matrix in O(n^2)

	for (int i = 0; i < SIZE - 2; i++) {
		for (int j = 0; j < SIZE - 2; j++) {
			if (i == j) {//we are on diagonal
				matrix[i][j] = 0;
				//filing parallel diagonals
				matrix[i + 1][j] = matrix[i][j + 1] = a1;
				matrix[i + 2][j] = matrix[i][j + 2] = a2;
			}
		}
	}

	matrix[SIZE - 2][SIZE - 2] = matrix[SIZE - 1][SIZE - 1] = 0;
	matrix[SIZE - 1][SIZE - 2] = matrix[SIZE - 2][SIZE - 1] = a1;
}

inline void initVector(double *vector, double number) {
	for (int i = 0; i < SIZE; i++){
		vector[i] = number;
	}
}

inline void multipleMatrixByConstant(double **matrix, double number) {

}

inline void copyMatrix(double **InMat, double **OutMat) {
	for (int i = 0; i < SIZE; i++) {
		memcpy(OutMat[i], InMat[i], SIZE * sizeof(double));
	}
}

inline void multipleMatrixByDiagonalMatrix(double **matrix, double **diag, double ** tmpMatrix)
{
	const double multiplier = diag[0][0];//number we will be multiplying by
	//matrix we will return 
	copyMatrix(matrix, tmpMatrix);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			if(tmpMatrix[i][j] != 0)
				tmpMatrix[i][j] *= multiplier;
		}
	}
}

inline void multipleMatrixByVector(double *outVector, double ** tmpMatrix, double *InVector)
{
	 for (int i = 0; i < SIZE; i++){
		 for (int j = 0; j < SIZE; j++){
			 outVector[i] += (tmpMatrix[i][j] * InVector[j]);
		 }
	 }
}

inline void multipleVectorByDiagonalMatrix(double **diag, double *vector, double *res)
{
	double multiplier = diag[0][0];
	for(int i=0; i<SIZE; i++){
		res[i] = vector[i] * multiplier;
	}
}

inline void addVectors(double *vec1, double *vec2, double *result)
{
	for (int i = 0; i < SIZE; i++) {
		result[i] = vec1[i] + vec2[i];
	}
}

inline double calculateNorm(double *vector)
{
	double norm = 0;
	for(int i=0; i<SIZE; i++){
		norm += vector[i] * vector[i];
	}
	return sqrt(norm);
}

inline void print2D(double **matrix)
{
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			std::cout << matrix[i][j]<<" ";
		}
		std::cout << std::endl;
	}
}
inline void print1D(double *matrix)
{
	for(int i=0; i<SIZE; i++)
	{
		std::cout << matrix[i] << std::endl;
	}
}

//----------GAUSS-SEIDEL----------
inline void fillLowerMatrix(double **matrix, double **lower)
{
	for(int i=0; i<SIZE; i++){
		for(int j=0; j<SIZE; j++){
			if (i >= j)
				lower[i][j] = matrix[i][j];
		}
	}
}

inline void fillUpperMatrix(double **matrix, double **upper)
{
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			if (i < j)
				upper[i][j] = matrix[i][j];
		}
	}
}

inline void forwardSubstitution(double **matrix, double *vector, double *result)
{
	double sum;

	for (int i = 0; i < SIZE; i++) {
		sum = 0;
		for (int j = 0; j < i; j++) {
			sum = sum + matrix[i][j] * result[j];
		}
		result[i] = (vector[i] - sum) / matrix[i][i];
	}
}

inline void backSubstitution(double **matrix, double *vector, double *result)
{
	for(int i = SIZE -1; i>= 0; i--){
		double tmp = vector[i];
		for(int j=i; j<SIZE; j++){
			if (j != i)
				tmp -= matrix[i][j] * result[j];
		}
		result[i] = tmp / matrix[i][i];
	}
}

