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


#define WTF matrix_t
#define CONST
#define min(a,b) ((a<b) ? a : b)

CONST sint32 g_maxThreads = 7; 	   // Maximum number of threads this program do use ( dont need more ).
CONST sint32 g_truncateSize = 32;    // Size of a matrix side when it is small enough to run on standard multiplication.
CONST sint32 g_outputWidth = 9; 	   // Set the output width for elements.
CONST sint32 g_outputPrecision = 15;  // Matrix C elements output precision.
sint32 g_cacheBlockSize = 16;   	   // Block size for loop tilling/blocking.

void print_mat(WTF *x);
void print_2mats(WTF *x, WTF *y);

int nextHighest2Power(int number)
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

int prevHighest2Power(int number)
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

size_t calculate_heap_usage(sint32 lvl, sint32 max)
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

		heap_req += 10 * power7 * (sizeof(WTF) + sizeof(BBQ) * (1<<(max-i)));
	}

	return heap_req;
}

void CopyMatrix(CONST WTF *src_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, sint32 size, WTF *dest_mat)
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


void print_mat(WTF *x)
{
   sint32 i,j;
   printf("\n");
   for (i=0;i < x->size;i++)
   {
      for (j=0;j  < x->size;j++)
      {
         printf(" %3d",x->elem[(x->size)*(i) + (j)]);
      }
      printf("\n");
   }
}

void print_2mats(WTF *x, WTF *y)
{
	   sint32 i,j;
	   printf("\n");
	   for (i=0;i < x->size;i++)
	   {
	      for (j=0;j  < x->size;j++)
	      {
	         printf(" %3d",x->elem[(x->size)*(i) + (j)]);
	      }

	      printf("          ");

	      for (j=0;j  < y->size;j++)
	      {
	         printf(" %3d",y->elem[(y->size)*(i) + (j)]);
	      }
	      printf("\n");
	   }
}

void SubMatrix(CONST WTF *a_mat, CONST WTF *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
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


void AddMatrix(CONST WTF *a_mat, CONST WTF *b_mat, sint32 ini_i0, sint32 ini_j0, sint32 ini_i1, sint32 ini_j1, \
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





// Calculates matrix multiplication between matrix_A and matrix_B. Store the result in matrix_C.
// ini_i means starting index i and ini_j means starting index j.
// Looping through lines ini_i <= i <= ini_i+size and columns ini_j <= j <= ini_j+size from A and B
// Loop unroll and loop blocking applied for performance purposes
void MultiplyMatrix(CONST WTF *A, CONST WTF *B, WTF *C, sint32 size)
{
	sint32 x = 0;
   sint32 ini_i = 0, ini_j = 0;
	sint32 i, j, k, l;
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






// Corner case when n = 0 or n = 1
void handleCornerCase(CONST WTF *matrix_A, CONST WTF *matrix_B, WTF *matrix_C, sint32 size)
{
	if ( size == 1 ) {
	   matrix_C->elem[0] = matrix_A->elem[0] * matrix_B->elem[0];
	} else if ( size == 2 )
	{
	   // This size would crash multi core forced solution ???
	   MultiplyMatrix(matrix_A, matrix_B, matrix_C, size);
	}
}


WTF *GetSquareMatrix(sint32 n)
{
   WTF *ret;
   sint32 i;
   
   ret = (WTF *) malloc(n*n*sizeof(BBQ) + sizeof(WTF));
   
   ret->size = n;
   ret->elem = (BBQ*) (ret + 1);
   memset(ret->elem, 0, n*n*sizeof(BBQ));

   return ret;
}

void clear (WTF *X)
{
   sint32 i,j;
   sint32 size = X->size;
   for (i = 0; i < size; i++)
   {
      for (j = 0; j < size; j++)
      {
         X->elem[size*(i) + (j)] = 0;
      }
   }
}

// Calculates Strassen multiplication within the range: lines ini_i <= i <= ini_i+size and columns ini_j <= j <= ini_j+size
void strassenMultiplication(CONST WTF *matrix_A, CONST WTF *matrix_B, sint32 size,
                            WTF *matrix_C, WTF *matrix_resultA, WTF *matrix_resultB) {
    if (size <= g_truncateSize){ // Small enough to run standard multiplication
    	MultiplyMatrix(matrix_A, matrix_B, matrix_C, size);
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
        clear(matrix_M1);
        clear(matrix_M2);
        
//        static unsigned int cnt = 0;
//        cnt++;
//        printf("\n %5d .    10 x (overhead [%2lu] + data(%3d^2)[%4lu] )%4lu", cnt, sizeof(WTF),half_size, half_size*half_size*sizeof(BBQ), half_size*half_size*sizeof(BBQ) + sizeof(WTF));

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
        clear(matrix_M1);
        
        // Calculating M2 = (A21 + A22) * (B11)
        AddMatrix(submatrix_A21, submatrix_A22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        CopyMatrix(submatrix_B11, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, mid_i, 0, half_size, matrix_C);
        clear(matrix_M1);
        
        // Calculating M3 = (A11) * (B12 - B22)
        CopyMatrix(submatrix_A11, 0, 0, 0, 0, half_size, matrix_resultA);
        SubMatrix(submatrix_B12, submatrix_B22, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, 0, mid_j, half_size, matrix_C);
        clear(matrix_M1);
        
        // Calculating M6 = (A21 - A11) * (B11 + B12)
        SubMatrix(submatrix_A21, submatrix_A11, 0, 0, 0, 0, 0, 0, half_size, matrix_resultA);
        AddMatrix(submatrix_B11, submatrix_B12, 0, 0, 0, 0, 0, 0, half_size, matrix_resultB);
        strassenMultiplication(matrix_resultA, matrix_resultB, half_size, matrix_M1, matrix_resultA, matrix_resultB);
        
        CopyMatrix(matrix_M1, 0, 0, mid_i, mid_j, half_size, matrix_C);
        clear(matrix_M1);
        
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
        
        clear(matrix_M1);
        
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
	sint32 array_size = 8;
	sint32 direct_comp_size = 4;

	sint32 high2power, low2power, recursion_lvl;

//	struct rlimit hehe;
//	hehe.rlim_cur = 30000000;
//	hehe.rlim_max = 30000000;
//	setrlimit(RLIMIT_DATA, &hehe);


	if (argc > 1) array_size = atoi(argv[1]);
	if (argc > 2) g_truncateSize = atoi(argv[2]);

	high2power = nextHighest2Power(array_size);
	low2power  = prevHighest2Power(g_truncateSize);
	recursion_lvl = high2power - low2power;printf("\n this cfg will recurse %d times", recursion_lvl);

	array_size = 1 << high2power;
	printf("\n N   = %d \n", array_size);
	g_truncateSize = 1 << low2power;
	printf("\n min = %d \n", g_truncateSize);

	unsigned long mem_req;
	mem_req = calculate_heap_usage(recursion_lvl, high2power);
	printf ("\n Glorious heap usage = %lu \n", mem_req + 3*1040 + 2*272);

    sint32 m, k, n, cacheBlockSize, dataType;

	
	// if (!readDimensions(&m, &k, &n, &cacheBlockSize, &dataType)) { 
		// printf("Couldn't read in.txt file! \n");
		// return 0;	
	// }	
	g_cacheBlockSize = 16;//cacheBlockSize;
	
	// Timing start.
	//start_timing(&usage, &user_start, &system_start);
		
	// Get matrices dimension.
	sint32 size = array_size;
	
   WTF *matrix_A = GetSquareMatrix(array_size);
   WTF *matrix_B = GetSquareMatrix(array_size);
   WTF *matrix_C = GetSquareMatrix(array_size);
   WTF *matrix_resultA, *matrix_resultB;

   matrix_A->size = array_size;
   //matrix_A->elem = MA;
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
    	handleCornerCase(matrix_A, matrix_B, matrix_C, size); // Corner cases, n = 0 or n = 1
	} 
   else 
   {
      //MultiplyMatrix(matrix_A, matrix_B, matrix_C, size);
	   matrix_resultA = GetSquareMatrix(size/2); // Buffer for addition, subtraction and copy results
	   matrix_resultB = GetSquareMatrix(size/2); // Buffer for addition, subtraction and copy results
	   strassenMultiplication(matrix_A, matrix_B, size, matrix_C, matrix_resultA, matrix_resultB);
   }
    
//    print_mat(matrix_A);
//    printf("\n     x     \n");
//    print_mat(matrix_B);
//    printf("\n     =     \n");
//    print_mat(matrix_C);
    
    // // Processing time end.
	// #ifdef FORCE_SINGLE_CORE
    	// end_timing(&usage, &user_end, &system_end);
    	// getTiming_User_System(&process_time, &user_start, &user_end, &system_start, &system_end);
    // #else
    	// gettimeofday(&time_end, NULL);
    	// getTiming(&process_time, &time_start, &time_end);
    // #endif
	
	// // Write processing time.
	// writeProcessTime(process_time);
	
	// // Timing start.
	// start_timing(&usage, &user_start, &system_start);
    
    // writeOutput(matrix_C, size, size);
    
    // // Timing end.
	// end_timing(&usage, &user_end, &system_end);
	
	// // Get Time.
	// getTiming_User_System(&write_time, &user_start, &user_end, &system_start, &system_end);
	
	// // Write writing time.
	// writeWriteTime(write_time);

    return 0;
}
