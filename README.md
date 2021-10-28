# SLAM
This is the implementation of our paper "SLAM: Efficient Sweep Line Algorithms for Kernel Density Visualization", which is submitted to SIGMOD 2022.

Compiling and running our code:
We use the shell script file "compile_and_run_SLAM.sh" to compile our code in Cygwin. In this shell script file, we also describe different input parameters and provide one example for running our code. Please check this file for details. Moreover, the zip file "Z_order_data_sampling.zip" includes the implementation of the Z-order method. We first use this file to generate the reduced dataset from each dataset (e.g., Seattle). Then, we use the RQS_{ball} method (the fastest exact method) to generate KDV for each reduced dataset to obtain the KDV result.
