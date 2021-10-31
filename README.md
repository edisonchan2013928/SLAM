# SLAM
This is the implementation of our paper "SLAM: Efficient Sweep Line Algorithms for Kernel Density Visualization", which is submitted to SIGMOD 2022.

**Compiling and running our code:**

We use the shell script file "compile_and_run_SLAM.sh" to compile our code with Cygwin (i.e., Windows OS). In this shell script file, we describe different input parameters and provide an example for running our code. Please check this shell script file for details. 

Moreover, the zip file "Z_order_data_sampling.7z" includes the implementation of the Z-order method. We first use this file to generate the reduced dataset from each dataset (e.g., Seattle). Then, we use the RQS_{ball} method (the fastest exact method) to generate KDV for each reduced dataset to obtain the approximate KDV result. As a remark, the shell script file "sampling.sh" (in the "Z_order_data_sampling.7z" file) compiles the code for this Z-order method, includes the explanation of each input variable, and provides an example for running our code. Please also check this shell script file for details.

Our code is written in C++ and no specific library has been adopted in our code (only the C++ standard template library (STL) is used).
