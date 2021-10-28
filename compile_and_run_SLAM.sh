#Compile our code
g++ -c init_visual.cpp -w -o init_visual.o
g++ -c Euclid_Bound.cpp -w -o Euclid_Bound.o
g++ -c Validation.cpp -w -o Validation.o
g++ -c SS_visual.cpp -w -o SS_visual.o
g++ -c kd_tree.cpp -w -o kd_tree.o
g++ -c ball_tree.cpp -w -o ball_tree.o
g++ -c alg_visual.cpp -w -o alg_visual.o
g++ -c SLAM.cpp -w -o SLAM.o
g++ main.cpp -O3 -o main init_visual.o Euclid_Bound.o Validation.o SS_visual.o kd_tree.o ball_tree.o alg_visual.o SLAM.o
exit

#Explanation of the input parameters in this code
#char*dataFileName = (char*)argv[1]; #Dataset file name

#stat.outMatrixFileName = (char*)argv[2]; #Output file name of the KDV result

#stat.b_para = atof(argv[3]); #This parameter is used to vary the bandwidth value by multiplying the default bandwidth with this value (default value is 1)

#stat.method = atoi(argv[4]);  #SCAN (Method = 0), aKDE (Method = 1), QUAD (Method = 2), RQS_{kd} (Method = 3), RAQ_{ball} (Method = 4),  
							   #SLAM_{SORT} (Method = 5), SLAM_{BUCKET} (Method = 6), SLAM_{SORT}^{(RAO)} (Method = 7), SLAM_{BUCKET}^{(RAO)} (Method = 8)

#stat.n_row = atoi(argv[5]); #Number of pixels in each row

#stat.n_col = atoi(argv[6]); #Number of pixels in each column

#stat.kernel_type = atoi(argv[7]); #Uniform kernel (kernel_type = 0), Epanechnikov kernel (kernel_type = 1), quartic kernel (kernel_type = 2)

#stat.epsilon = atof(argv[8]); #This parameter is only used in the aKDE and QUAD methods (By default: epsilon = 0.05)


#Example of running our code (with the SLAM_{BUCKET}^{(RAO)} method in the Seattle dataset)
dir="../Datasets/"
out_dir="./Results/"
n_row=1280
n_col=960
kernel_type=1
epsilon=0.05
b_para=1
select_region=1

method=8
dataset="Seattle"
timeout 14400 ./main $dir$dataset"/"$dataset $out_dir$dataset"_M"$method"_r"$n_row"_c"$n_col $b_para $method $n_row $n_col $kernel_type $epsilon