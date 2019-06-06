#include <iostream>
#include "functions.h"
#include <fstream>
#include <ctime>
using namespace std;


int main()
{
	const double a1 = 13.0;//13
	const double a2 = -1.0;
	const double a3 = -1.0;

	ofstream jacobiIter;
	ofstream gaussIter;


	double **matrix = initialization();//inicjalization of matrix with given parameters
	fillMatrix(matrix, a1, a2, a3);

	//vecotr b (right side of bacis equation)
	double vectorB[SIZE];
	fillVectorB(vectorB);
	
	//our system of equations, so lets start with Jacobian.
	//the base we need are 3 matrix, lower, upper and diagonal.
	//if we study the equation:
	//r = [-D^(-1)]*(L+U)*r + [-D^(-1)]*b 
	//we are able to se that computing lower and upper matrix separately is redundant.
	//Much more efficient way is just to create matrix without the diagonal from orginal matrix

	//in equation we need 2 reversed diagonal matrix (one negated one normal)
	double **revDiagonal = initialization();
	fillRevDiagonal(revDiagonal, a1);

	//next sum of lower and upper matrix
	double **lowerUpper = initialization();
	fillLowerUpper(lowerUpper, -a2, -a3);

	//result vector r
	double result[SIZE];
	initVector(result, 1.0/SIZE);

	double residuum[SIZE];
	initVector(residuum, 0);
	const double prec = 1.0e-9;
	double norm = 1.0;

	double **tmpMatrix = initialization();
	double tmpLeftSide[SIZE];
	initVector(tmpLeftSide, 0);
	double tmpRightSide[SIZE];

	//tmpMatrix = -D^(-1)*(L+U) <- matrix multiplayed by negated reversed diagonal 
	multipleMatrixByDiagonalMatrix(lowerUpper, revDiagonal, tmpMatrix);
	//tmpRightSide = D^(-1)*b
	multipleVectorByDiagonalMatrix(revDiagonal, vectorB, tmpRightSide);

	cout << "Metoda Jacobiego" << endl;
	jacobiIter.open("jacobiIter.txt", std::ofstream::out | std::ofstream::trunc);
	clock_t begin = clock();
	while(norm > prec)
	{
		//tmpLeftSide = -D^(-1)*(L+U)*r <- vector
		initVector(tmpLeftSide, 0);
		multipleMatrixByVector(tmpLeftSide, tmpMatrix, result);
		//result = tmpLeftSide + tmpRightSide
		addVectors(tmpLeftSide, tmpRightSide, result);

		initVector(residuum, 0);
		//residuum = M*r
		multipleMatrixByVector(residuum, matrix, result);
		//residuum = M*r - b <- vector
		for(int i=0; i<SIZE; i++){
			residuum[i] -= vectorB[i];
		}
		norm = calculateNorm(residuum);
		jacobiIter << norm << ", ";
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << elapsed_secs << endl;
	jacobiIter.close();

	//----------METODA GAUSSA-SEIDLA----------
	//we will try to solve our set of equations using Gauss-Seidl method, the pattern is:
	//r = -(D+L)^(-1)*(U*r) + (D+L)^(-1)*b
	//we see that we should revers 2 matrix, but it's very naive due to very high RAM and performence costs, so
	//we will use forward substitution twice (actually once per iteration, coz is (D+L)^(-1)*b is const) and add results per iterations 

	//lower matrix contains diagonal as is it D + L
	double **lowerD = initialization();
	fillLowerMatrix(matrix, lowerD);

	//rightSide = (D+L)^(-1)*b <- vector
	double rightSide[SIZE];
	forwardSubstitution(lowerD, vectorB, rightSide);

	//U matrix
	double **upper = initialization();
	fillUpperMatrix(matrix, upper);

	double leftSidePart[SIZE];
	initVector(tmpLeftSide, 0);

	//give previous values 
	initVector(residuum, 0);
	norm = 1;
	initVector(result, 1.0 / SIZE);

	cout << "Metoda Gaussa-Seidla" << endl;
	gaussIter.open("gaussIter.txt", std::ofstream::out | std::ofstream::trunc);
	clock_t begin1 = clock();
	while (norm > prec)
	{
		//leftSideParr = U*r <- vector
		initVector(leftSidePart, 0);
		multipleMatrixByVector(leftSidePart, upper, result);

		//tmpLeftSide = (D+L)^(-1)*(U*r) <- vector
		forwardSubstitution(lowerD, leftSidePart, tmpLeftSide);
		//negate all signs so tmpLeftSide = -(D+L)^(-1)*(U*r)
		for (double& i : tmpLeftSide){
			i = -i;
		}
		//result = -(D+L)^(-1)*(U*r) + (D+L)^(-1)*b
		addVectors(tmpLeftSide, rightSide, result);
		initVector(residuum, 0);
		//residuum = M*r
		multipleMatrixByVector(residuum, matrix, result);
		//residuum = M*r - b <- vector
		for (int i = 0; i < SIZE; i++) {
			residuum[i] -= vectorB[i];
		}
		norm = calculateNorm(residuum);
		cout << norm << endl;
		gaussIter << norm << ", ";
	}
	clock_t end1 = clock();
	double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;
	cout << elapsed_secs1 << endl;
	gaussIter.close();

	//ZADANIE D - FAKTORYZACJA LU
	cout << "Faktoryzacja LU" << endl;

	//fill low with 0 and its diag with 1
	//up only 0, copping main matrix now is wasting time 
	double **low = initialization();
	for (int i = 0; i < SIZE; i++)
		low[i][i] = 1.0;
	double **up = initialization();

	double y[SIZE];
	initVector(y, 1);

	double x[SIZE];
	initVector(x, 1);

	clock_t begin2 = clock();
	for(int j=0; j<SIZE; j++){
		for(int i=0; i<SIZE; i++){
			up[i][j] += matrix[i][j];
			for (int k = 0; k <= i - 1; k++)
				up[i][j] -= low[i][k] * up[k][j];
		}

		for(int i= j+1; i<SIZE; i++){
			for (int k = 0; k <= j - 1; k++)
				low[i][j] -= low[i][k] * up[k][j];

			low[i][j] += matrix[i][j];
			low[i][j] /= up[j][j];
		}
	}

	forwardSubstitution(low, vectorB, y);
	backSubstitution(up, y, x);
	clock_t end2 = clock();
	double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
	cout << elapsed_secs2 << endl;

	initVector(residuum, 0);
	multipleMatrixByVector(residuum, matrix, x);
	for (int i = 0; i < SIZE; i++) {
		residuum[i] -= vectorB[i];
	}
	norm = calculateNorm(residuum);
	cout << norm << endl;

	//delete allocated memory
	for (int i = 0; i < SIZE; ++i) {
		delete[] matrix[i];
		delete[] revDiagonal[i];
		delete[] lowerUpper[i];
		delete[] lowerD[i];
		delete[] upper[i];
		delete[] low[i];
		delete[] up[i];
	}
	delete[] matrix;
	delete[] revDiagonal;
	delete[] lowerUpper;
	delete[] lowerD;
	delete[] upper;
	delete[] low;
	delete[] up;
	return 0;
}