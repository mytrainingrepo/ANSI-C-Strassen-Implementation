
#include "MatrixOps.h"

int main(int argc, char *argv[])
{
	int size;
	double optimised_time, naive_time;
	struct timeval system_start, system_end, user_start, user_end; // Structs used by rusage
	struct rusage usage;

	// 1st argument is 'n', determining the matrice's 'n*n' size
	if (argc > 1) size=atoi(argv[1]);
	else size = 1024;

	matrix_t *matrix_A = GenerateSimpleMatrix(size);     // generates a matrix containing consecutive values
	matrix_t *matrix_B = GenerateIdentityMatrix(size);   // generates an identity matrix
	matrix_t *matrix_C = GenerateZeroMatrix(size);       // empty matrix where the 1st result will be stored
	matrix_t *matrix_D = GenerateZeroMatrix(size);       // empty matrix where the 2nd result will be stored


	StartTiming(&usage, &user_start, &system_start);
	// calculate A + A*B using the Strassen algorithm
	ComputeResultOptimised(matrix_A, matrix_B, matrix_C);
	StopTiming(&usage, &user_end, &system_end);
	ComputeTotalTime(&optimised_time, &user_start, &user_end, &system_start, &system_end);

	StartTiming(&usage, &user_start, &system_start);
	// calculate A + A*B using the basic iterative algorithm
	ComputeResultNaive(matrix_A, matrix_B, matrix_D);
	StopTiming(&usage, &user_end, &system_end);
	ComputeTotalTime(&naive_time, &user_start, &user_end, &system_start, &system_end);

	if ((argc > 2) && (*argv[2]=='v'))
	{
		printf("\n A + (A x B) \n A = \n");
		PrintMatrix(matrix_A);
		printf("\n B = \n");
		PrintMatrix(matrix_B);
		printf("\n C (optimised) = \n");
		PrintMatrix(matrix_C);
		printf("\n D (naive) = \n");
		PrintMatrix(matrix_D);
		printf("\n");
	}

	// check results of the two equivalent operations
	if (CompareMatrix(matrix_C, 0, 0, 0, 0, size, matrix_D)) printf("\n Results match : algorithm is correct. \n");
	else printf("\n Results don't match : algorithm is incorrect. \n");

	// show both timings so their efficiency can be compared
	printf("\n Optimised algorithm took %f milliseconds \n Naive algorithm took %f milliseconds \n", optimised_time, naive_time);

	return 0;
}
