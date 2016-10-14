

/*********************************************************************************
 * Includes
 *********************************************************************************/
#include "MatrixOps.h"

/*********************************************************************************
 * Global variables
 *********************************************************************************/

sint32 g_truncateSize = 32;    // Size of a matrix side when it is small enough to run on standard multiplication.
sint32 g_cacheBlockSize = 16;  // Block size for loop tilling/blocking.

void *MemoryCurrent;
void *MemoryStart;

/*********************************************************************************
 * Function declarations
 *********************************************************************************/

LOCAL_INLINE sint32 FindNextPowerOf2(sint32 number);
LOCAL_INLINE sint32 FindPrevPowerOf2(sint32 number);
LOCAL_INLINE size_t CalculateHeapSize(sint32 lvl, sint32 max);

LOCAL_INLINE void CopyMatrix(matrix_t *src_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, matrix_t *dest_mat);
LOCAL_INLINE void SubMatrix(matrix_t *a_mat, matrix_t *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  matrix_t *res_mat);
LOCAL_INLINE void AddMatrix(matrix_t *a_mat, matrix_t *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  matrix_t *res_mat);

LOCAL_INLINE void MultiplyMatrixOptimized(matrix_t *A, matrix_t *B, matrix_t *C, sint32 size);
LOCAL_INLINE void MultiplyMatrixTiny(matrix_t *matrix_A, matrix_t *matrix_B, matrix_t *matrix_C);
LOCAL_INLINE void AllocateAll (size_t size);
LOCAL_INLINE void FreeAll (void);
LOCAL_INLINE matrix_t *GetSquareMatrix(sint32 n);
LOCAL_INLINE void ZeroMatrix (matrix_t *A);

/*********************************************************************************
 * Function definitions
 *********************************************************************************/

void GetElapsedTime(double* timing, struct timeval* time_start, struct timeval* time_end)
{
	*timing = 0;

	*timing = (double)(time_end->tv_sec- time_start->tv_sec)*1000;
	*timing += (double)(time_end->tv_usec - time_start->tv_usec)/1000;

	return;
}

/*  */
void StartTiming(struct rusage* usage, struct timeval* user_start, struct timeval* system_start)
{
	getrusage(RUSAGE_SELF, usage);
    *user_start = usage->ru_utime;
    *system_start = usage->ru_stime;
}

void StopTiming(struct rusage* usage, struct timeval* user_end, struct timeval* system_end)
{
	getrusage(RUSAGE_SELF, usage);
    *user_end = usage->ru_utime;
    *system_end = usage->ru_stime;
}

void ComputeTotalTime(double* timing, struct timeval* time_start1, struct timeval* time_end1,
										   struct timeval* time_start2, struct timeval* time_end2)
{
	double aux;
	GetElapsedTime(&aux, time_start1, time_end1);
	GetElapsedTime(timing, time_start2, time_end2);
	*timing += aux;
}

LOCAL_INLINE sint32 FindNextPowerOf2(sint32 number)
{
	int count=0;
	int comp = number;
	while(number!=0)
	{
		number=number>>1;
		count++;
	}

	if ((comp & -comp ) == comp)
	{
		return (count-1);
	}
	else
	{
		return count;
	}
}

LOCAL_INLINE sint32 FindPrevPowerOf2(sint32 number)
{
	int count=0;
	int comp = number;
	while(number!=0)
	{
		number=number>>1;
		count++;
	}
	return (count-1);
}

LOCAL_INLINE size_t CalculateHeapSize(sint32 lvl, sint32 max)
{
	sint32 i,j;
	size_t heap_req = 0;
	sint32 power7 = 0;

	for (i = 0; i< lvl; ++i)
	{
		power7 = 1;

		for (j = 0; j < i; ++j)
		{
			power7 *=7;
		}

		heap_req += 10 * power7 * \
					(sizeof(matrix_t) + \
							sizeof(DATA_TYPE) * \
							((1<<(max-i-1)) * (1<<(max-i-1)))  );
	}

	// temporary result storage matrices
	heap_req += 2 * \
			(sizeof(matrix_t) + \
					sizeof(DATA_TYPE) * \
					((1<<(max-1)) * (1<<(max-1)))  );

	// top level operand and result matrices
	heap_req += 3 * \
			(sizeof(matrix_t) + \
					sizeof(DATA_TYPE) * \
					((1<<(max)) * (1<<(max)))  );
	return heap_req;
}

sint32 CompareMatrix(matrix_t *A, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, matrix_t *B)
{
	sint32 size_A = A->size;
	sint32 size_B = B->size;

	sint32 i,j;
	for (i = 0; i < size; ++i)
	{
		for (j = 0; j < size; ++j)
	    {
			if (B->elem[size_B*(i+ini_i1) + (j+ini_j1)] != A->elem[size_A*(i+ini_i0) + j+ini_j0]) return 0;
		}
	}
	return 1;
}

LOCAL_INLINE void CopyMatrix(matrix_t *src_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, matrix_t *dest_mat)
{
	sint32 size_src = src_mat->size;
	sint32 size_dest = dest_mat->size;

	sint32 i,j;
	for (i = 0; i < size; ++i)
	{
		for (j = 0; j < size; ++j)
	    {
			dest_mat->elem[size_dest*(i+ini_i1) + (j+ini_j1)] = src_mat->elem[size_src*(i+ini_i0) + j+ini_j0];
		}
	}
	return;
}


LOCAL_INLINE void PrintMatrix(matrix_t *A)
{
   sint32 i,j;
   printf("\n");
   for (i=0;i < A->size;i++)
   {
      for (j=0;j  < A->size;j++)
      {
         printf(" %3ld", (A->elem[(A->size)*(i) + (j)]) );
      }
      printf("\n");
   }
}

LOCAL_INLINE void SubMatrix(matrix_t *a_mat, matrix_t *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  matrix_t *res_mat)
{
	sint32 size_a = a_mat->size;
	sint32 size_b = b_mat->size;
	sint32 size_res = res_mat->size;
	sint32 i,j;
	for (i = 0; i < size; ++i)
   {
    	for (j = 0; j < size; ++j)
      {
    		res_mat->elem[size_res*(i+ini_i2) + (j+ini_j2)] = a_mat->elem[size_a*(i+ini_i2) + (j+ini_j2)] - b_mat->elem[size_b*(i+ini_i1) + (j+ini_j1)];
    	}
	}	   		   
	return;		   
}


LOCAL_INLINE void AddMatrix(matrix_t *a_mat, matrix_t *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  matrix_t *res_mat)
{
	sint32 size_a = a_mat->size;
	sint32 size_b = b_mat->size;
	sint32 size_res = res_mat->size;
	sint32 i,j;
	for (i = 0; i < size; ++i)
    {
    	for (j = 0; j < size; ++j)
        {
    		res_mat->elem[size_res*(i+ini_i2) + (j+ini_j2)] = a_mat->elem[size_a*(i+ini_i0) + (j+ini_j0)] + b_mat->elem[size_b*(i+ini_i1) + (j+ini_j1)];
    	}
	}		
	return;		   
}	


void MultiplyMatrixNaive(matrix_t *A, matrix_t *B, matrix_t *C)
{
	int i, j, k;
	int size_A = A->size, size_B = B->size, size_C = C->size;

	for (i = 0; i < size_A; i++)
	{
		for (j = 0; j < size_B; j++)
		{
			C->elem[i*size_C + j] = 0;
	    	for (k = 0; k < size_A; k++)
	    	{
	        	C->elem[i*size_C + j] += A->elem[i*size_A + k] * B->elem[k*size_B + j];
	      	}
	  	}
   	}
   	return;
}


// Calculates matrix multiplication between matrix_A and matrix_B. Store the result in matrix_C.
// Loop unroll and loop blocking applied for performance purposes
LOCAL_INLINE void MultiplyMatrixOptimized(matrix_t *A, matrix_t *B, matrix_t *C, sint32 size)
{
	sint32 x = 0;
    sint32 ini_i = 0, ini_j = 0;
	sint32 i, j, k;
	sint32 limit0 = size + ini_i; // Index i limit
	sint32 limit1 = size + ini_j; // Index j limit
	sint32 limit2 = size;         // Index k limit
    sint32 aux_i, aux_j, aux_k;
    sint32 aux_limit_i; 		   // Block index limit i
    sint32 aux_limit_j; 		   // Block index limit j
    sint32 aux_limit_k; 		   // Block index limit k
    sint32 unroll_factor = 5;
    sint32 unroll_limit; 		        // Loop unroll index limit
    DATA_TYPE acc0, acc1, acc2, acc3, acc4; // Accumulators, eliminate data dependencies

    for (i = ini_i; i < limit0; i += g_cacheBlockSize) 
    {
        // Blocking index i limit
        aux_limit_i = min((i+g_cacheBlockSize), limit0);
        
    	for (j = ini_j; j < limit1; j += g_cacheBlockSize) 
      {
    	    // Blocking index j limit
    	    aux_limit_j = min((j+g_cacheBlockSize), limit1);
    	    
    		for (k = 0; k < limit2; k += g_cacheBlockSize) 
         {
    		    // Blocking index k limit
    		    aux_limit_k = min((k+g_cacheBlockSize), limit2);
    		    
                unroll_limit = aux_limit_k-4; // Unrolling by factor of 5
                
              	for(aux_i = i; aux_i < aux_limit_i; ++aux_i) 
               {
                	for(aux_j = j; aux_j < aux_limit_j; ++aux_j) 
                  {

                    	acc0 = 0; acc1 = 0; acc2 = 0; acc3 = 0; acc4 = 0;
						
						sint32 endA = aux_i*(A->size) + ini_j;
						
						// Unrolling for k loop
                    	for(aux_k = k; aux_k < unroll_limit; aux_k+=unroll_factor) 
                     {
                        	acc0 += A->elem [endA+aux_k]   * B->elem [(B->size)*(aux_k+ini_i)   + (aux_j)];
                        	acc1 += A->elem [endA+aux_k+1] * B->elem [(B->size)*(aux_k+ini_i+1) + (aux_j)];
                        	acc2 += A->elem [endA+aux_k+2] * B->elem [(B->size)*(aux_k+ini_i+2) + (aux_j)];
                        	acc3 += A->elem [endA+aux_k+3] * B->elem [(B->size)*(aux_k+ini_i+3) + (aux_j)];
                        	acc4 += A->elem [endA+aux_k+4] * B->elem [(B->size)*(aux_k+ini_i+4) + (aux_j)];
                     } 
                        
                        sint32 endC = aux_i*(C->size) + aux_j;
                        
                        // Gather possible uncounted elements
                        for (; aux_k < aux_limit_k; ++aux_k)
                        	C->elem[endC] += A->elem[endA+aux_k] * B->elem [(B->size)*(aux_k+ini_i) + (aux_j)];
                       
                        // Sum up everything
                        C->elem[endC] += acc0 + acc1 + acc2 + acc3 + acc4;

                	}  
            	}   
			}
		}
	}
    return;
}



LOCAL_INLINE void MultiplyMatrixTiny(matrix_t *matrix_A, matrix_t *matrix_B, matrix_t *matrix_C)
{
	matrix_C->elem[0] = matrix_A->elem[0] * matrix_B->elem[0];
}

LOCAL_INLINE void AllocateAll (size_t size)
{
	// allocate a chunk of memory that can fit all objects allocated
	//  during a a call of the ComputeResultOptimized() function
	MemoryStart = (void *) calloc(size, 1);
	MemoryCurrent = MemoryStart;
	if (MemoryStart == NULL) perror("Memory Allocation failed:");
	return;
}

LOCAL_INLINE void FreeAll (void)
{
	free(MemoryStart);
}

LOCAL_INLINE matrix_t *GetSquareMatrix(sint32 n)
{
   matrix_t *ret;
   sint32 i;

   ret = (matrix_t *) MemoryCurrent;
   MemoryCurrent += n*n*sizeof(DATA_TYPE) + sizeof(matrix_t);
   
   ret->size = n;
   ret->elem = (DATA_TYPE*) (ret + 1);

   return ret;
}

LOCAL_INLINE void ZeroMatrix (matrix_t *A)
{
   sint32 size = A->size;
   memset(A->elem,0,size*size*sizeof(*(A->elem)));
}


void MultiplyMatrixStrassen(matrix_t *matrix_A, matrix_t *matrix_B, sint32 size,
                            matrix_t *matrix_C, matrix_t *matrix_resultA, matrix_t *matrix_resultB) {
    if (size <= g_truncateSize){ // Small enough to run standard multiplication
    	MultiplyMatrixOptimized(matrix_A, matrix_B, matrix_C, size);
    } else {
        sint32 half_size = size/2;		   // Size for each four quadrant
        sint32 mid_i = half_size; // Index middle of side1
        sint32 mid_j = half_size; // Index middle of side2
        sint32 end_i = size;	   // Index end limit of side1
        sint32 end_j = size; 	   // Index end limit of side2
        sint32 i,j;

        // Matrix A quadrants
        matrix_t *submatrix_A11 = GetSquareMatrix(half_size);
        matrix_t *submatrix_A12 = GetSquareMatrix(half_size);
        matrix_t *submatrix_A21 = GetSquareMatrix(half_size);
        matrix_t *submatrix_A22 = GetSquareMatrix(half_size);
      
        // Matrix B quadrants
        matrix_t *submatrix_B11 = GetSquareMatrix(half_size);
        matrix_t *submatrix_B12 = GetSquareMatrix(half_size);
        matrix_t *submatrix_B21 = GetSquareMatrix(half_size);
        matrix_t *submatrix_B22 = GetSquareMatrix(half_size);

        // Matrices M_i are calculated and stored directly into matrix_C.
		// We save up some memory by using the same matrix_M1 as a buffer for each
		// computation and then filling it with zeros for the next calculation
		// Two auxiliar matrices are needed to hold M4 and M5 values
        matrix_t *matrix_M1 = GetSquareMatrix(half_size);
        matrix_t *matrix_M2 = GetSquareMatrix(half_size);
        ZeroMatrix(matrix_M1);
        ZeroMatrix(matrix_M2);
        
        // Getting submatrices from A and B
        for (i = 0; i < half_size; ++i)
        {
        	sint32 endM = half_size*i;
 			for (j = 0; j < half_size; ++j)
 			{ //[n*(i) + (j)]
 				submatrix_A11->elem[endM + j] = matrix_A->elem[(matrix_A->size)*(i) + (j)];
 				submatrix_A12->elem[endM + j] = matrix_A->elem[(matrix_A->size)*(i) + (j+half_size)];
 				submatrix_A21->elem[endM + j] = matrix_A->elem[(matrix_A->size)*(i+half_size) + (j)];
 				submatrix_A22->elem[endM + j] = matrix_A->elem[(matrix_A->size)*(i+half_size) + (j+half_size)];

 				submatrix_B11->elem[endM + j] = matrix_B->elem[(matrix_B->size)*(i) + (j)];
 				submatrix_B12->elem[endM + j] = matrix_B->elem[(matrix_B->size)*(i) + (j+half_size)];
 				submatrix_B21->elem[endM + j] = matrix_B->elem[(matrix_B->size)*(i+half_size) + (j)];
 				submatrix_B22->elem[endM + j] = matrix_B->elem[(matrix_B->size)*(i+half_size) + (j+half_size)];
 			}
		}

		// Calculating M1 = (A11 + A22) * (B11 + B22)
        AddMatrix(submatrix_A11, submatrix_A22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        AddMatrix(submatrix_B11, submatrix_B22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        // Copying directly to matrix_C (answer matrix)
        CopyMatrix(matrix_M1, 0, 0, 0, 0, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculating M2 = (A21 + A22) * (B11)
        AddMatrix(submatrix_A21, submatrix_A22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        CopyMatrix(submatrix_B11, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, mid_i, 0, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculating M3 = (A11) * (B12 - B22)
        CopyMatrix(submatrix_A11, 0, 0, 0, 0, half_size, matrix_resultA);
        SubMatrix(submatrix_B12, submatrix_B22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, 0, mid_j, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculating M6 = (A21 - A11) * (B11 + B12)
        SubMatrix(submatrix_A21, submatrix_A11, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        AddMatrix(submatrix_B11, submatrix_B12, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, mid_i, mid_j, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculate C22 = M1 - M2 + M3 + M6
        AddMatrix(matrix_C, matrix_C, mid_i, mid_j, 0, 0, mid_i, mid_j, half_size, matrix_C);     // C22 = M1 + M6
        AddMatrix(matrix_C, matrix_C, 0, mid_j, mid_i, mid_j, mid_i, mid_j, half_size, matrix_C); // C22 = M3 + C22
        SubMatrix(matrix_C, matrix_C, mid_i, mid_j, mid_i, 0, mid_i, mid_j, half_size, matrix_C); // C22 = C22 - M2
          
        // Calculating M4 = (A22) * (B21 - B11)
        CopyMatrix(submatrix_A22, 0, 0, 0, 0, half_size, matrix_resultA);
        SubMatrix(submatrix_B21, submatrix_B11, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        // Calculate C21 = M2 + M4
        AddMatrix(matrix_C, matrix_M1, mid_i, 0, 0, 0, mid_i, 0, half_size, matrix_C);
        
		// Not resetting auxiliar matrix_M1 ! Holding M4 value...

        // Calculating M5 = (A11 + A12) * (B22)
        AddMatrix(submatrix_A11, submatrix_A12, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        CopyMatrix(submatrix_B22, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M2, matrix_resultA, matrix_resultB);
        // Calculate C12 = M3 + M5
        AddMatrix(matrix_C, matrix_M2, 0, mid_j, 0, 0, 0, mid_j, half_size, matrix_C);
        
        // auxiliar matrix_M1 = M4, auxiliar matrix_M2 = M5 
        // M4 = M4 - M5
        SubMatrix(matrix_M1, matrix_M2, 0, 0, 0, 0, 0, 0, half_size, matrix_M1);
        
        // Calculating partially C11 = M1 + M4 - M5
        AddMatrix(matrix_C, matrix_M1, 0, 0, 0, 0, 0, 0, half_size, matrix_C);
        
        ZeroMatrix(matrix_M1);
        
        // Calculating M7 = (A12 - A22) * (B21 + B22)
        SubMatrix(submatrix_A12, submatrix_A22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        AddMatrix(submatrix_B21, submatrix_B22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        MultiplyMatrixStrassen(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        // Calculate C11 = M1 + M4 - M5 + M7
        AddMatrix(matrix_C, matrix_M1, 0, 0, 0, 0, 0, 0, half_size, matrix_C);
    }
    return;
}

matrix_t * GenerateIdentityMatrix(sint32 size)
{
	matrix_t *ret;
	sint32 i,j;

	ret = (matrix_t *) calloc(size*size*sizeof(DATA_TYPE) + sizeof(matrix_t),1);
	ret->size = size;
	ret->elem = (DATA_TYPE*) (ret + 1);

	for (i = 0; i < size; ++i)
	{
	   for (j = 0; j < size; ++j)
	   {
		   if (i==j) ret->elem[i*size + j] = 1;
	   }
	}
	return ret;
}

matrix_t * GenerateSimpleMatrix(sint32 size)
{
	matrix_t *ret;
	sint32 i,j;
	ret = (matrix_t *) calloc(size*size*sizeof(DATA_TYPE) + sizeof(matrix_t),1);
	ret->size = size;
	ret->elem = (DATA_TYPE*) (ret + 1);

	for (i = 0; i < size; ++i)
	{
	   for (j = 0; j < size; ++j)
	   {
		   ret->elem[i*size + j] = i*size + j;
	   }
	}
	return ret;
}

matrix_t * GenerateZeroMatrix(sint32 size)
{
	matrix_t *ret;
	ret = (matrix_t *) calloc(size*size*sizeof(DATA_TYPE) + sizeof(matrix_t),1);
	ret->size = size;
	ret->elem = (DATA_TYPE*) (ret + 1);
	return ret;
}

// compute A + (A*B) using an optimised algorithm
void ComputeResultOptimised(matrix_t *src_A, matrix_t *src_B, matrix_t *dest_C)
{
	sint32 init_size = dest_C->size;
	sint32 high2power, low2power, recursion_lvl, new_size;
	size_t mem_req;

	if (init_size == 1)
	{
		MultiplyMatrixTiny(src_A, src_B, dest_C);
		return;
	}

	// determine the top power of 2 that fits the init_size
	high2power = FindNextPowerOf2(init_size);
	low2power  = FindPrevPowerOf2(g_truncateSize);
	// based on the maximum size of the matrix to be computed
	//  and the minimum size of a matrix that doesn't a recursive
	//  call, determine max level of recursion of the
	//  MultiplyMatrixStrassen() function
	recursion_lvl = high2power - low2power;

	// new matrix size is higher or equal to the initial size
	//  and must be a power of two
	new_size = 1 << high2power;

	// based on the number of dynamic memory allocations,
	//  the level of recursion and the size of the objects
	//  determine how much heap memory the program needs
	mem_req = CalculateHeapSize(recursion_lvl, high2power);
	// allocate at once all memory required
	AllocateAll(mem_req);

	// allocate the matrices that will fit and pad the original inputs and output
	matrix_t *matrix_A = GetSquareMatrix(new_size);
	matrix_t *matrix_B = GetSquareMatrix(new_size);
	matrix_t *matrix_C = GetSquareMatrix(new_size);
	matrix_t *matrix_resultA, *matrix_resultB;

	CopyMatrix(src_A, 0, 0, 0, 0, init_size, matrix_A); // padding A until size is power of two
	CopyMatrix(src_B, 0, 0, 0, 0, init_size, matrix_B); // padding B until size is power of two

	matrix_resultA = GetSquareMatrix(new_size/2); // Buffer for addition, subtraction and copy results
	matrix_resultB = GetSquareMatrix(new_size/2); // Buffer for addition, subtraction and copy results

	// calculate the result of the matrix multiplication recursively
	MultiplyMatrixStrassen(matrix_A, matrix_B, new_size, matrix_C, matrix_resultA, matrix_resultB);
	// perform last addition operation
	AddMatrix(matrix_A, matrix_C, 0, 0, 0, 0, 0, 0, new_size, matrix_C);

	// copying results back to the original size matrix
	CopyMatrix(matrix_C, 0, 0, 0, 0, init_size, dest_C);

	// free all heap memory alocated previously
    FreeAll();
    return;
}

// compute A + (A*B) in a naive fashion
void ComputeResultNaive(matrix_t *src_A, matrix_t *src_B, matrix_t *dest_C)
{
	MultiplyMatrixNaive(src_A, src_B, dest_C);
	AddMatrix(src_A, dest_C, 0, 0, 0, 0, 0, 0, dest_C->size, dest_C);
	return;
}
