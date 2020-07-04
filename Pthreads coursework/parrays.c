#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>


//Argument structure used to pass arguments to the function called by pthread_create
struct args{
	int thread_num; 	//Current thread number/id
	int total_threads;	//Total number of threads
	int dim;			//array dimension
	double precision;	//precision to work towards
	int *in_precision;	//variable acting as a flag to know when to stop loop
	double **num_array;	//pointer to array
	double **new_array;	//pointer to array
	pthread_barrier_t *barr; //pointer to barrier
};



//calc_column function is the function called by the pthread_create which calculates
//the average value for each column for the thread calling the function
void *calc_column(void *arg_ptr)
{
	struct args *arg = (struct args *)arg_ptr;
	int y,x;
	
	//This while loop will loop until all threads reach 
	//the while and the in_precision variable is 0
	do{
		//synchronise all barriers before setting the in_precision variable to 0
		pthread_barrier_wait(arg->barr);
		*arg->in_precision = 0;
		
		//Nested for loop which iterates the columns of the array
		for (x = arg->thread_num; x < arg->dim - 1;  x=x+arg->total_threads)
		{
			for (y =1 ; y < arg->dim - 1; y++)
			{
				//Calculate average for current slot
				arg->new_array[x][y] = (arg->num_array[x + 1][y] + arg->num_array[x - 1][y] + arg->num_array[x][y + 1] + arg->num_array[x][y - 1]) / 4;
				
				//If in_precision is still 0 and the previous number is not 
				//within the precision given set variable to 1
				if (  fabs(arg->num_array[x][y] - arg->new_array[x][y]) > arg->precision && (*arg->in_precision == 0) )
				{
					//Set variable to false
					*arg->in_precision = 1;
				}
				
			}
		}
		
		//Swap pointers of the old and new array
		double **temp_address = arg->new_array;
		arg->new_array = arg->num_array;
		arg->num_array = temp_address;
		
		//Synchronise threads before checking the end of the loop
		pthread_barrier_wait(arg->barr);
	}while(*arg->in_precision >= 1 );
	
	return(NULL);
}

//This function creates the pthreads and barriers and 
//calls the calculate column function thorugh pthread_create
void array_relax(double **arr, double **arr_copy, int dimension, int thrdcount, double precision)
{
	int i;
	
	//Allocate memory for the variable
	int *precision_flag = malloc(sizeof(int));
	*precision_flag = 0;
	
	//Declare threads
	pthread_t thread_ids[thrdcount];
	pthread_barrier_t barrier;
	
	//Initialise barrier
	pthread_barrier_init(&barrier, NULL, thrdcount);
	
	//Allocate memory for the argument structure
	struct args *arg = malloc(sizeof(struct args) * thrdcount);
	
	//Create array of arguments one for each thread num.
	//thread num is used to know what columns to process once the calc_column fucntion is called
	for (i = 0; i < thrdcount; i++)
	{
		arg[i].num_array = arr;
		arg[i].new_array = arr_copy;
		arg[i].dim = dimension;
		arg[i].in_precision = precision_flag;
		arg[i].precision = precision;
		arg[i].total_threads = thrdcount;
		arg[i].thread_num = i+1;
		arg[i].barr = &barrier;

		
		pthread_create(&thread_ids[i], NULL, calc_column, &arg[i]);
	}
	
	//Join threads
	for (i=0; i < thrdcount; i++) {
    	pthread_join(thread_ids[i], NULL);
	}
	//Destroy barrier
	pthread_barrier_destroy(&barrier);

	//Free memory
	free(precision_flag);
	free(arg);
}


int main()
{
	
	//Number of threads to start with
	int thrdcount = 1;
	
	//Up to how many threads to test the program
	int threads_to_test = 8;
	
	//Precision of the array average to work towards
	double precision = 0.00001;
	
	//Dimension of the array
	int dimension = 150;
	
	//Initialise counter variables
	int y,x,i;
	
	//Allocate memory for the two arrays used in the program and one for
	//reseting the array's values between tests
	double **arr = (double**)calloc(sizeof(*arr) * dimension, 1);
	double **arr_copy = (double**)calloc(sizeof(*arr_copy) * dimension, 1);
	double **test = (double**)calloc(sizeof(*arr) * dimension, 1);
	for ( i = 0; i < dimension; ++i)
	{
		arr[i] = (double*)calloc(sizeof(**arr) * dimension, 1);
		arr_copy[i] = (double*)calloc(sizeof(**arr_copy) * dimension, 1);
		test[i] = (double*)calloc(sizeof(**test) * dimension, 1);
	}

	
	//Assign random double value to the arrays
	for (y = 0; y < dimension; y++)
	{
		for (x = 0; x < dimension; x++)
		{
			int new_num = rand() % 10;
			arr[x][y] = new_num;
			arr_copy[x][y] = new_num;
			test[x][y] = new_num;
		}
	}


	int k ;
	
	//Iterate through the program 16 times, one for each thread, 
	//each time testing the runtime of the program and outputting the time taken
	for (k = thrdcount; k < threads_to_test+1; k++){
	
		//time measurement
		struct timespec start, finish;
		clock_gettime(CLOCK_MONOTONIC, &start);		
		
		//Call the fucntion
		array_relax(arr, arr_copy, dimension, k, precision);
		
		clock_gettime(CLOCK_MONOTONIC, &finish);
		if(start.tv_nsec > finish.tv_nsec){
			finish.tv_nsec += 1000000000;
			finish.tv_sec--;
		}
		//Print time taken
		printf("Array size: %d Threads used: %d TIME TAKEN: %ld.%09ld  \n",dimension,k, (long)(finish.tv_sec-start.tv_sec),finish.tv_nsec - start.tv_nsec);
	
		//Recopy arrays for next execution
		for (y = 0; y < dimension; y++)
		{
			for (x = 0; x < dimension; x++)
			{
				arr[x][y]=test[x][y];
				arr_copy[x][y]=test[x][y];
			}
		}
	}

	//Free arrays
	for (i = 0; i < dimension; i++)
	{
		free(arr[i]);
		free(arr_copy[i]);
		free(test[i]);
	}
	free(arr);
	free(arr_copy);
	free(test);


	return 0;
};







