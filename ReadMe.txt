
Test Task. Fast Matrix Multiplication.
======================================


Description
======================================

The program seeks to solve the A+(A*B) matrix operation as fast as possible, using the Strassen algorithm for multiplication.

Platforms tested:
======================================
PC (Intel I7 / x86_64) - standard GCC 5.4.0 compiler
BeagleBone Black (TI AM335 / ARM Cortex A8) - GCC compiler (arm-linux-muslgnueabihf-gcc)


Implementation:
======================================
The program is based on an efficient method of multiplying two matrices using the Strassen algorithm presented by Eduardo Ferreira in a research project ("Code Optimization Techniques") completed at Federal University of Rio de Janeiro.

Test program walkthrough (main() from Test.c):
	1. The program's 1st represents the matrix' dimensions

	2. (optional) if 2nd argument is 'v' the input test matrices and output
	    of the operation are printed on standard output

	3. Two simple matrices are generated:
		identity matrix, having 1 on the main diagonal and 0 everywhere else
		a simple matrix storing all the numbers from 0 to its total size - 1
	   The multiplication and addition can be easily observed from the output matrix

	4. Calculate the A + S(A,B), store the result in matrix_C and determine the total 
	    time (spent in user and kernel space) it took for the CPU to finish the oepration.

	5. Calculate A + B*C using simplest iterative method and determine the time

	6. (optional) print the two input matrices and the two output matrices, 
	    the latter must be identical if the algorithms functioned correctly.

	7. Display the status message and the timings of the two equivalent operations.

Computation function walkthrough (ComputeResultOptimised from MatrixOps.c):
	1. if matrix size = 1, basically the matrix is a scalar, perform the operation directly

	2. fit the input matrices inside a provisional set of matrices, which can be divided
	    recursively into 4 equal sub-matrices. The difference is padded with 0s.

	3. By knowing the matrix size and minimal size of a matrix couple we can multiply directly,
	    we can find out how many times we must split the input matrixes recursively, and thus
	    the recursion level.

	4. Calculate the exact amount of memory needed and allocate it all.

	5. Initialize the temporary matrices

	6. Run the Strassen algorithm and add the A matrix to the resulting C matrix

	7. free all memory

Known shortcomings:
======================================
- SIMD not supported for any architecture (x86 SSE or ARM NEON)
- only float data type supported
- data alignment not taken into consideration
- if the matrix size not a power of two, the matrix is temporarily padded 
  with zeros so that Strassen's recursion works seamlessly. Depending on the
  number, this can slow the program down dramatically.
- potentially a very large amount of data is allocated dynamically at one time,
  depending on system resources, for very large matrices the allocation might fail.
