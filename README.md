# DataSample-and-Test-Suite
This project is for generating a sequence that passs all the NIST test suite,  TestU01 batteries and Hamano et al.'s new test,
but will be rejected by our proposed DM-1 and DM-2 test.

##1. Generate the Sample Sequence 

   The code for generating the sample data is in the file "GenerateSample.c", in which the data file "rdseed.data" is used as input.  Compile this single c source file and run it, then you will get the output file "biggap10_4_1", which is the same as "sample.data".

##2. Test the Sample Data by NIST and TestU01 Test Suite

   Unzip these two compression packs，there are README or INSTALL files in both of them. Follow the README or INSTALL files and you will get the corresponding executable files after compiling. Then  you may test the generated file "biggap10_4_1". Notice that the size of this data file is 10^9 bits. For the NIST test suite, the input sequence size should be 10^6, and the number of sequences should be 1000.  While for TestU01 batteries, the input size is 10^9.

##3. Test the Generated Sample Data by NIST Linear Complexity Test, HSY Test and our proposed Deviation Measure Tests

   These source files are in the directory ./LinearComplexityRandomTest exclude CRandTest.c. Enter this directory and compile these c source files exclude CRandTest.c. You can firstly modify the suffix name of CRandTest.c, for example by using the command
            
            mv CRandTest.c CRandTest.cpp

And then compile all the C files,

            gcc -o SampleTest *.c -lm -fopenmp 
   
   Notice that the input parameter "-fopenmp" in this command is for parallelization,  since there are some codes for parallelization in  "LC_Random_Test.c" by the instruction of omp. Then, there will be an executable file named "FileTest", and you can test the sample data by
           
           ./SampleTest ../biggap10_4_1 TestResult.txt
   
   Experimental results of these 4 test methods for "biggap_10_4_1"  are in the file "TestResult.txt", including both P and U results, and the number of rejected sequences according to the criteria of U.
   
##4. Test the traditional rand function in the C library by NIST Linear Complexity Test, HSY Test and our proposed Deviation Measure Tests

   These source files are in the directory ./LinearComplexityRandomTest exclude SampleTest.c. Enter this directory and compile these c source files exclude SampleTest.c. You can firstly modify the suffix name of SampleTest.c, for example by using the command
            
            mv SampleTest.c SampleTest.cpp

And then compile all the C files,

            gcc -o CRandTest *.c -lm -fopenmp 
   
   Then, there will be an executable file named "CRandTest", and you can test the Rand function by
           
           ./CRandTest  1  1RandSeed100.txt
The first input parameter of the program CRandTest is the position where we choose the bit from the return value of rand() each time. The legal value of first parameter is from 1 to 15, and the second input parameter is the seed file for srand function, which can be uncompressed from file 10GroupsRandomSeeds.rar in the upper folder.   
  
   Experimental results of CRandTest for position x are in the generated file "C_Rand_posx_result.txt"(for instance, if x is 1, then the corresponding file is "C_Rand_pos1_result.txt"), including both P and U results, and the number of rejected sequences according to the criteria of U.
   