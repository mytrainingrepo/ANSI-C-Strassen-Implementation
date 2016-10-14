

#ifndef MATRIXOPS_H_
#define MATRIXOPS_H_

/*********************************************************************************
 * Includes
 *********************************************************************************/
#include "Compiler_Cfg.h"
#include "Platform_Types.h"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

/*********************************************************************************
 * Type declarations
 *********************************************************************************/

typedef sint32 DATA_TYPE;

typedef struct matrix
{
	sint32 size;
	DATA_TYPE *elem;
} matrix_t;

/*********************************************************************************
 * Macros
 *********************************************************************************/

#define min(a,b) ((a<b) ? a : b)

#define ACCESS(mat, i, j)  (mat->elem[i*(mat->size)+j])

/*********************************************************************************
 * Function declarations
 *********************************************************************************/

void ComputeResultOptimised(matrix_t *src_A, matrix_t *src_B, matrix_t *dest_C);
void ComputeResultNaive();
void PrintMatrix(matrix_t *A);
void MultiplyMatrixNaive(matrix_t *A, matrix_t *B, matrix_t *C);
void MultiplyMatrixStrassen(matrix_t *matrix_A, matrix_t *matrix_B, sint32 size,
                            matrix_t *matrix_C, matrix_t *matrix_resultA, matrix_t *matrix_resultB);
sint32 CompareMatrix(matrix_t *A, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, matrix_t *B);

void GetElapsedTime(double* timing, struct timeval* time_start, struct timeval* time_end);
void StartTiming(struct rusage* usage, struct timeval* user_start, struct timeval* system_start);
void StopTiming(struct rusage* usage, struct timeval* user_end, struct timeval* system_end);
void ComputeTotalTime(double* timing, struct timeval* time_start1, struct timeval* time_end1,
										   struct timeval* time_start2, struct timeval* time_end2);
matrix_t * GenerateIdentityMatrix(sint32 size);
matrix_t * GenerateSimpleMatrix(sint32 size);
matrix_t * GenerateZeroMatrix(sint32 size);


#endif /* MATRIXOPS_H_ */
