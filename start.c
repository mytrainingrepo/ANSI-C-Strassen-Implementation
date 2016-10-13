

/*********************************************************************************
 * Includes
 *********************************************************************************/
#include "Compiler_Cfg.h"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>


typedef int sint32;
typedef sint32 BBQ;

typedef struct matrix
{
	sint32 size;
	BBQ *elem;
} matrix_t;

// __attribute__((always_inline))
#define WTF matrix_t
#define CONST
#define min(a,b) ((a<b) ? a : b)

/*********************************************************************************
 * Global variables
 *********************************************************************************/


CONST sint32 g_maxThreads = 7; 	   // Maximum number of threads this program do use ( dont need more ).
CONST sint32 g_truncateSize = 32;    // Size of a matrix side when it is small enough to run on standard multiplication.
CONST sint32 g_outputWidth = 9; 	   // Set the output width for elements.
CONST sint32 g_outputPrecision = 15;  // Matrix C elements output precision.
sint32 g_cacheBlockSize = 16;   	   // Block size for loop tilling/blocking.

void *MemoryCurrent;
void *MemoryStart;

/*********************************************************************************
 * Function declarations
 *********************************************************************************/


LOCAL_INLINE void PrintMatrix(WTF *A);
LOCAL_INLINE void GetElapsedTime(double* timing, struct timeval* time_start, struct timeval* time_end);
LOCAL_INLINE void StartTiming(struct rusage* usage, struct timeval* user_start, struct timeval* system_start);
LOCAL_INLINE void StopTiming(struct rusage* usage, struct timeval* user_end, struct timeval* system_end);
LOCAL_INLINE void ComputeTotalTime(double* timing, struct timeval* time_start1, struct timeval* time_end1,
										   struct timeval* time_start2, struct timeval* time_end2);
LOCAL_INLINE int FindNextPowerOf2(int number);
LOCAL_INLINE int FindPrevPowerOf2(int number);
LOCAL_INLINE size_t CalculateHeapSize(sint32 lvl, sint32 max);
LOCAL_INLINE void CopyMatrix(CONST WTF *src_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, WTF *dest_mat);
LOCAL_INLINE void SubMatrix(CONST WTF *a_mat, CONST WTF *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  WTF *res_mat);
LOCAL_INLINE void AddMatrix(CONST WTF *a_mat, CONST WTF *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  WTF *res_mat);
LOCAL_INLINE void MultiplyMatrixNaive(WTF *A, WTF *B, WTF *C);
LOCAL_INLINE void MultiplyMatrixOptimized(CONST WTF *A, CONST WTF *B, WTF *C, sint32 size);
LOCAL_INLINE void MultiplyMatrixTiny(CONST WTF *matrix_A, CONST WTF *matrix_B, WTF *matrix_C, sint32 size);
LOCAL_INLINE void AllocateAll (size_t size);
LOCAL_INLINE void FreeAll (void);
LOCAL_INLINE WTF *GetSquareMatrix(sint32 n);
LOCAL_INLINE void ZeroMatrix (WTF *A);

void strassenMultiplication(CONST WTF *matrix_A, CONST WTF *matrix_B, sint32 size,
                            WTF *matrix_C, WTF *matrix_resultA, WTF *matrix_resultB);

/*********************************************************************************
 * Function definitions
 *********************************************************************************/

LOCAL_INLINE void GetElapsedTime(double* timing, struct timeval* time_start, struct timeval* time_end)
{
	*timing = 0;

	*timing = (double)(time_end->tv_sec- time_start->tv_sec)*1000;
	*timing += (double)(time_end->tv_usec - time_start->tv_usec)/1000;

	return;
}

LOCAL_INLINE void StartTiming(struct rusage* usage, struct timeval* user_start, struct timeval* system_start)
{
	getrusage(RUSAGE_SELF, usage);
    *user_start = usage->ru_utime;
    *system_start = usage->ru_stime;
}

LOCAL_INLINE void StopTiming(struct rusage* usage, struct timeval* user_end, struct timeval* system_end)
{
	getrusage(RUSAGE_SELF, usage);
    *user_end = usage->ru_utime;
    *system_end = usage->ru_stime;
}

LOCAL_INLINE void ComputeTotalTime(double* timing, struct timeval* time_start1, struct timeval* time_end1,
										   struct timeval* time_start2, struct timeval* time_end2)
{
	double aux;
	GetElapsedTime(&aux, time_start1, time_end1);//printf("\n\n     Time spent in user mode   = %f \n", aux);
	GetElapsedTime(timing, time_start2, time_end2);//printf("\n     Time spent in kernel mode = %f \n", *timing);
	*timing += aux;
	//printf("\n          Total time spent     = %f \n", *timing);
}

LOCAL_INLINE int FindNextPowerOf2(int number)
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

LOCAL_INLINE int FindPrevPowerOf2(int number)
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
	size_t heap_req = 0;
	int power7 = 0;

	for (int i = 0; i< lvl; ++i)
	{
		power7 = 1;

		for (int j = 0; j < i; ++j)
		{
			power7 *=7;
		}

		heap_req += 10 * power7 * \
					(sizeof(WTF) + \
							sizeof(BBQ) * \
							((1<<(max-i-1)) * (1<<(max-i-1)))  );
	}


	heap_req += 2 * \
			(sizeof(WTF) + \
					sizeof(BBQ) * \
					((1<<(max-1)) * (1<<(max-1)))  );

	heap_req += 3 * \
			(sizeof(WTF) + \
					sizeof(BBQ) * \
					((1<<(max)) * (1<<(max)))  );
	return heap_req;
}

LOCAL_INLINE void CopyMatrix(CONST WTF *src_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, WTF *dest_mat)
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


LOCAL_INLINE void PrintMatrix(WTF *A)
{
   sint32 i,j;
   printf("\n");
   for (i=0;i < A->size;i++)
   {
      for (j=0;j  < A->size;j++)
      {
         printf(" %3d",A->elem[(A->size)*(i) + (j)]);
      }
      printf("\n");
   }
}

LOCAL_INLINE void SubMatrix(CONST WTF *a_mat, CONST WTF *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  WTF *res_mat)
{
	sint32 size_a = a_mat->size;
	sint32 size_b = b_mat->size;
	sint32 size_res = res_mat->size;
	for (sint32 i = 0; i < size; ++i)
   {
    	for (sint32 j = 0; j < size; ++j)
      {
    		res_mat->elem[size_res*(i+ini_i2) + (j+ini_j2)] = a_mat->elem[size_a*(i+ini_i2) + (j+ini_j2)] - b_mat->elem[size_b*(i+ini_i1) + (j+ini_j1)];
    	}
	}	   		   
	return;		   
}


LOCAL_INLINE void AddMatrix(CONST WTF *a_mat, CONST WTF *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
               sint32 ini_i2, sint32 ini_j2, sint32 size,  WTF *res_mat)
{
	sint32 size_a = a_mat->size;
	sint32 size_b = b_mat->size;
	sint32 size_res = res_mat->size;

	for (sint32 i = 0; i < size; ++i)
    {
    	for (sint32 j = 0; j < size; ++j)
        {
    		res_mat->elem[size_res*(i+ini_i2) + (j+ini_j2)] = a_mat->elem[size_a*(i+ini_i0) + (j+ini_j0)] + b_mat->elem[size_b*(i+ini_i1) + (j+ini_j1)];
    	}
	}		
	return;		   
}	


LOCAL_INLINE void MultiplyMatrixNaive(WTF *A, WTF *B, WTF *C)
{
	int i, j, k;// i*(A->size) + j
	int size_A = A->size, size_B = B->size, size_C = C->size;
	for (i = 0; i < size_A; i++)
	{
		for (j = 0; j < size_B; j++)
		{
	    	for (k = 0; k < size_A; k++)
	    	{
	        	C->elem[i*size_C + j] += A->elem[i+size_A + k]* B->elem[k*size_B + j];
	      	}
	  	}
   	}
   	return;
}


// Calculates matrix multiplication between matrix_A and matrix_B. Store the result in matrix_C.
// ini_i means starting index i and ini_j means starting index j.
// Looping through lines ini_i <= i <= ini_i+size and columns ini_j <= j <= ini_j+size from A and B
// Loop unroll and loop blocking applied for performance purposes
LOCAL_INLINE void MultiplyMatrixOptimized(CONST WTF *A, CONST WTF *B, WTF *C, sint32 size)
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
    BBQ acc0, acc1, acc2, acc3, acc4; // Accumulators, eliminate data dependencies

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






// n = 0 or n = 1
LOCAL_INLINE void MultiplyMatrixTiny(CONST WTF *matrix_A, CONST WTF *matrix_B, WTF *matrix_C, sint32 size)
{
	if ( size == 1 )
	{
	   matrix_C->elem[0] = matrix_A->elem[0] * matrix_B->elem[0];
	}
	else if ( size == 2 )
	{
	   MultiplyMatrixOptimized(matrix_A, matrix_B, matrix_C, size);
	}
}

LOCAL_INLINE void AllocateAll (size_t size)
{
	MemoryStart = (void *) malloc(size);
	MemoryCurrent = MemoryStart;
	memset(MemoryStart, 0, size); //??
	return;
}

LOCAL_INLINE void FreeAll (void)
{
	free(MemoryStart);
}

LOCAL_INLINE WTF *GetSquareMatrix(sint32 n)
{
   WTF *ret;
   sint32 i;
   
   // individual allocation
   //ret = (WTF *) malloc(n*n*sizeof(BBQ) + sizeof(WTF));

   ret = (WTF *) MemoryCurrent;
   MemoryCurrent += n*n*sizeof(BBQ) + sizeof(WTF);
   
   ret->size = n;
   ret->elem = (BBQ*) (ret + 1);

   return ret;
}

LOCAL_INLINE void ZeroMatrix (WTF *A)
{
   sint32 size = A->size;
   memset(A->elem,0,size*size*sizeof(*(A->elem)));
}

// Calculates Strassen multiplication within the range: lines ini_i <= i <= ini_i+size and columns ini_j <= j <= ini_j+size
void strassenMultiplication(CONST WTF *matrix_A, CONST WTF *matrix_B, sint32 size,
                            WTF *matrix_C, WTF *matrix_resultA, WTF *matrix_resultB) {
    if (size <= g_truncateSize){ // Small enough to run standard multiplication
    	MultiplyMatrixOptimized(matrix_A, matrix_B, matrix_C, size);
    } else {
        sint32 half_size = size/2;		   // Size for each four quadrant
        sint32 mid_i = half_size; // Index middle of side1
        sint32 mid_j = half_size; // Index middle of side2
        sint32 end_i = size;	   // Index end limit of side1
        sint32 end_j = size; 	   // Index end limit of side2

        // Matrix A quadrants
        WTF *submatrix_A11 = GetSquareMatrix(half_size);
        WTF *submatrix_A12 = GetSquareMatrix(half_size);
        WTF *submatrix_A21 = GetSquareMatrix(half_size);
        WTF *submatrix_A22 = GetSquareMatrix(half_size);
      
        // Matrix B quadrants
        WTF *submatrix_B11 = GetSquareMatrix(half_size);
        WTF *submatrix_B12 = GetSquareMatrix(half_size);
        WTF *submatrix_B21 = GetSquareMatrix(half_size);
        WTF *submatrix_B22 = GetSquareMatrix(half_size);

        // Matrices M_i are calculated and stored directly into matrix_C.
		// We save up some memory by using the same matrix_M1 as a buffer for each
		// computation and then filling it with zeros for the next calculation
		// Two auxiliar matrices are needed to hold M4 and M5 values
        WTF *matrix_M1 = GetSquareMatrix(half_size);
        WTF *matrix_M2 = GetSquareMatrix(half_size);
        ZeroMatrix(matrix_M1);
        ZeroMatrix(matrix_M2);
        
        // Getting submatrices from A and B
        for (sint32 i = 0; i < half_size; ++i)
        {
        	sint32 endM = half_size*i;
 			for (sint32 j = 0; j < half_size; ++j)
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
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        // Copying directly to matrix_C (answer matrix)
        CopyMatrix(matrix_M1, 0, 0, 0, 0, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculating M2 = (A21 + A22) * (B11)
        AddMatrix(submatrix_A21, submatrix_A22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        CopyMatrix(submatrix_B11, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, mid_i, 0, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculating M3 = (A11) * (B12 - B22)
        CopyMatrix(submatrix_A11, 0, 0, 0, 0, half_size, matrix_resultA);
        SubMatrix(submatrix_B12, submatrix_B22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, 0, mid_j, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculating M6 = (A21 - A11) * (B11 + B12)
        SubMatrix(submatrix_A21, submatrix_A11, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        AddMatrix(submatrix_B11, submatrix_B12, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, mid_i, mid_j, half_size, matrix_C);
        ZeroMatrix(matrix_M1);
        
        // Calculate C22 = M1 - M2 + M3 + M6
        AddMatrix(matrix_C, matrix_C, mid_i, mid_j, 0, 0, mid_i, mid_j, half_size, matrix_C);     // C22 = M1 + M6
        AddMatrix(matrix_C, matrix_C, 0, mid_j, mid_i, mid_j, mid_i, mid_j, half_size, matrix_C); // C22 = M3 + C22
        SubMatrix(matrix_C, matrix_C, mid_i, mid_j, mid_i, 0, mid_i, mid_j, half_size, matrix_C); // C22 = C22 - M2
          
        // Calculating M4 = (A22) * (B21 - B11)
        CopyMatrix(submatrix_A22, 0, 0, 0, 0, half_size, matrix_resultA);
        SubMatrix(submatrix_B21, submatrix_B11, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        // Calculate C21 = M2 + M4
        AddMatrix(matrix_C, matrix_M1, mid_i, 0, 0, 0, mid_i, 0, half_size, matrix_C);
        
		// Not resetting auxiliar matrix_M1 ! Holding M4 value...

        // Calculating M5 = (A11 + A12) * (B22)
        AddMatrix(submatrix_A11, submatrix_A12, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        CopyMatrix(submatrix_B22, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M2, matrix_resultA, matrix_resultB);   
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
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        // Calculate C11 = M1 + M4 - M5 + M7
        AddMatrix(matrix_C, matrix_M1, 0, 0, 0, 0, 0, 0, half_size, matrix_C);
    }
    return;
}



sint32 main(sint32 argc, char *argv[])
{	
    double read_time, process_time, write_time;
	struct timeval system_start, system_end, user_start, user_end; // Structs used by rusage
	struct rusage usage;
	struct timeval time_start, time_end; // Structs gettimeofday

	sint32 array_size = 8;
	sint32 direct_comp_size = 4;

	sint32 high2power, low2power, recursion_lvl;

	if (argc > 1) array_size = atoi(argv[1]);
	if (argc > 2) g_truncateSize = atoi(argv[2]);

	StartTiming(&usage, &user_start, &system_start);

	high2power = FindNextPowerOf2(array_size);
	low2power  = FindPrevPowerOf2(g_truncateSize);
	recursion_lvl = high2power - low2power;

	array_size = 1 << high2power;
	g_truncateSize = 1 << low2power;

	size_t mem_req;
	mem_req = CalculateHeapSize(recursion_lvl, high2power);
	AllocateAll(mem_req);

    sint32 m, k, n, cacheBlockSize, dataType;

	
	// if (!readDimensions(&m, &k, &n, &cacheBlockSize, &dataType)) { 
		// printf("Couldn't read in.txt file! \n");
		// return 0;	
	// }	
	g_cacheBlockSize = 16;//cacheBlockSize;
	
	// Timing start.
	//StartTiming(&usage, &user_start, &system_start);
		
	// Get matrices dimension.
	sint32 size = array_size;
	
   WTF *matrix_A = GetSquareMatrix(array_size);
   WTF *matrix_B = GetSquareMatrix(array_size);
   WTF *matrix_C = GetSquareMatrix(array_size);
   WTF *matrix_resultA, *matrix_resultB;

   matrix_A->size = array_size;
   for (sint32 i = 0; i < array_size; ++i)
   {
	   for (sint32 j = 0; j < array_size; ++j)
	   {
		   matrix_A->elem[i*matrix_A->size + j] = 1;//i*matrix_A->size + j;   // boring content
		   //if (i==j) matrix_A->elem[i*matrix_A->size + j] = 1;          // unit matrix
	   }
   }

   matrix_B->size = array_size;
   //matrix_B->elem = MU;
   for (sint32 i = 0; i < array_size; ++i)
   {
	   for (sint32 j = 0; j < array_size; ++j)
	   {
		   //if (i==j) matrix_B->elem[i*matrix_B->size + j] = 1;          // unit matrix
		   matrix_B->elem[i*matrix_B->size + j] = 1;
	   }
   }

    // Read Matrices A and B from input
    // if (readInput(matrix_A, matrix_B, size, size, size) != 1)  // If return is not 1, something is wrong
		// return 0;

   if ( size <= 2 ) 
   {
	   MultiplyMatrixTiny(matrix_A, matrix_B, matrix_C, size); // Corner cases, n = 0 or n = 1
	} 
   else 
   {
	   matrix_resultA = GetSquareMatrix(size/2); // Buffer for addition, subtraction and copy results
	   matrix_resultB = GetSquareMatrix(size/2); // Buffer for addition, subtraction and copy results
	   strassenMultiplication(matrix_A, matrix_B, size, matrix_C, matrix_resultA, matrix_resultB);
   }

   StopTiming(&usage, &user_end, &system_end);
   ComputeTotalTime(&process_time, &user_start, &user_end, &system_start, &system_end);

   printf("\n          Total time spent     = %f \n", process_time);
    
//    PrintMatrix(matrix_A);
//    printf("\n     x     \n");
//    PrintMatrix(matrix_B);
//    printf("\n     =     \n");
//    PrintMatrix(matrix_C);
    


    FreeAll();
    return 0;
}
