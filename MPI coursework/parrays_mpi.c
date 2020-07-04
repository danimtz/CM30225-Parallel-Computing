
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//Relax array function that implements the basic math of averaging an array slot and calculating in_precision flag
int relax_array(double *row_array, double *row_array_copy, int dimension, int rows_per_proc, double precision, int in_precision){
    
    int i,j;
    for (i = 1; i < rows_per_proc - 1; i++) { 
		for (j = 1; j < dimension - 1; j++) {
		    
			// Store relaxed number into new array
			row_array_copy[i*dimension+j] = (row_array[(i-1)*dimension+j] + row_array[(i+1)*dimension+j] + row_array[i*dimension+(j-1)] + row_array[i*dimension+(j+1)]) / 4;
					  
						  
			//If not within precision set in_precision flag to false (0)
			if (in_precision == 1 && fabs(row_array[i*dimension+j] - row_array_copy[i*dimension+j]) > precision) {
				in_precision = 0;
			}
		}
	}
	
	//Return in_precision flag
	return in_precision;
}

int main(int argc, char** argv) {
    
    
    
    // Initialize the MPI environment
    int rc, world_size, proc_rank;
    rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}


    // Get the number of processesw
   
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    //Initialise MPI_status and MPI_request
    MPI_Status status;
    MPI_Request request;
    
    //Initialize variables
    int root,dimension, core_dim, proc_id, topmost_row, botmost_row; 
    double precision;
    
    //For timing
    double time_dif;
    
    dimension = 10;
    precision = 0.0001;
    
    int tag_1 = 1;
    int tag_2 = 2;
    //Core dim refers to the core of the array: array without the edge numbers
    core_dim = dimension-2; 
    
    root = 0;
    
    //How many rows each process will have access to.
	int rows_per_proc = (floor((core_dim) / world_size)) + 2;
		
	
		
    //Run by root process
    if(proc_rank == root){
        //Allocate memory for array with values
        double *array = (double *)malloc(dimension * dimension * sizeof(double));
        
		int i, j;
		//Assign values to the array
		for (i = 0; i < dimension; i++) {
			for (j = 0; j < dimension; j++) {
                
            	int new_num = rand() % 10;
                array[i*dimension+j] = new_num;
            	
			}
		}
		
		// Allocate memory for the rows stored in the processor before and after the averaging
		double *root_rows = (double *)malloc(dimension * (rows_per_proc+remainder_rows) * sizeof(double));
        double *root_rows_copy = (double *)malloc(dimension * (rows_per_proc+remainder_rows) * sizeof(double));
        
		//How many rows arent accounted for in rows_per_proc if (core_dim/world_size is not an integer)
        int remainder_rows = (core_dim) % world_size; 
        
        //Intialiase variables
        int first_time_setup = 0; 
        int in_precision = 0;
        
        //Start timer
        double tic = MPI_Wtime();
        
        //Loop until array is relaxed within precision
        while (!in_precision) {
			in_precision = 1;
			
			//First time set up assigns rows to each processor
			if (first_time_setup == 0) {
                first_time_setup = 1;
                
                
                
                //First assign rows+remainder rows to the root process
				for (i = 0; i < rows_per_proc+remainder_rows; i++) {
					for (j = 0; j < dimension; j++) {
					    
						root_rows[i*dimension+j] = array[i*dimension+j];
						root_rows_copy[i*dimension+j] = root_rows[i*dimension+j];
					}
				}
				
				botmost_row = (rows_per_proc+remainder_rows)-1;
	
                
                //Iterate through the rest of the processors and assign rows to each
				for (proc_id = 1; proc_id < world_size; proc_id++) {
					
					
                    //Calculate start and end rows for each processor. 
                    //Note: The very first and very last row per processor is shared with the previous and next processor but isnt altered
					topmost_row = botmost_row - 1; 
					botmost_row = topmost_row + (rows_per_proc - 1); 


					//Non-blocking send to the process
					int send_ammount;
					send_ammount = rows_per_proc*dimension;
					
					MPI_Isend( &array[topmost_row*dimension], send_ammount, MPI_DOUBLE, proc_id, tag_1, MPI_COMM_WORLD, &request );
					
				}
				
            }else{
                
                //Send slave process 1 the root process bottommost row and receive process 1 top most row.
				MPI_Isend( &root_rows[(rows_per_proc-2)*dimension], dimension , MPI_DOUBLE, root + 1, tag_1, MPI_COMM_WORLD, &request);
				
				MPI_Recv( &root_rows[(rows_per_proc-1)*dimension], dimension , MPI_DOUBLE, root + 1, tag_1, MPI_COMM_WORLD, &status);
            }
				
            //Perform array relaxation on the root rows
            in_precision = relax_array(root_rows, root_rows_copy, dimension, (rows_per_proc+remainder_rows), precision, in_precision);
            
            
    		// Ensure that processes have received the messages before moving on to the next section
            MPI_Wait(&request, &status);
            
            /* CHECK FOR IN PRECISION OF THE OTHER PROCESSES */
            int process_in_precision = 0;
    
    		//Get the in_precision flag from every process to checl if to keep going or not
    		for (proc_id = 1; proc_id < world_size; proc_id++) {
    			
    			MPI_Recv( &process_in_precision, 1, MPI_INT, proc_id, tag_2, MPI_COMM_WORLD, &status);
    
    			// If root process is in precision but any of the other processes are not, set flag to false (0)
    			if (in_precision && !process_in_precision) {
    				in_precision = 0;
    			}
            }
            
            // Message other processes whether to stop or not (If not in precision, keep going)
    		for (proc_id = 1; proc_id < world_size; proc_id++) {
    			MPI_Isend( &in_precision, 1 , MPI_INT, proc_id, tag_1, MPI_COMM_WORLD, &request);
            }
            
            // Swap the old and new arrays for the next iteration. 
            //(Similar to how this was implemented in the p_threads coursework, however here the array swapped is only for each process's rows and not global)
            
    		double *temp = root_rows;
    		root_rows = root_rows_copy;
    		root_rows_copy = temp;
    
            //Synchronise
            MPI_Wait(&request, &status);
        }
        
        //Precision is done
        //RECLAIM ALL ROWS FROM OTHER PROCESSES INTO ROOT PROCESS AND PRINT

        botmost_row = (rows_per_proc+remainder_rows)-2;
        
        for (proc_id = 1; proc_id < world_size; proc_id++) {
            
			topmost_row = botmost_row - 1; 
		    botmost_row = topmost_row + (rows_per_proc - 1); 

			// Receive the values from process number an_id, into tempVals
			int recv_ammount = (rows_per_proc-2)*dimension;
			
			MPI_Recv( &array[(topmost_row+2)*dimension], recv_ammount, MPI_DOUBLE, proc_id, tag_2, MPI_COMM_WORLD, &status);
			
		}
        
        
        //Copy updated values from the root_rows to the main array
        for (i = 1; i < (rows_per_proc+remainder_rows)-1; i++) { 
			for (j = 1; j < dimension - 1; j++) { 
				array[i*dimension+j] = root_rows[i*dimension+j];
			}
        }
        
        
        double toc = MPI_Wtime();
		time_dif = toc-tic;
		
		printf(" Array size: %d\n\n", dimension);
		printf(" Number of processors: %d\n\n", world_size);
		printf(" Time taken = %lf seconds\n\n", time_dif);

        
    }else{//Run by all other processes
    
        //Variable declaration for each process
        int i,j;
        int first_time_setup;
        int in_precision;
        
        first_time_setup = 1;
        in_precision = 0;

		//Allocate memory for the two arrays
		double *rows = (double*)malloc(rows_per_proc * dimension * sizeof(double));
        double *rows_copy = (double*)malloc(rows_per_proc * dimension * sizeof(double));
        
        
        while (!in_precision) {
            
            in_precision = 1;
            
			//Similar to the root process, on first time setup receive the rows from root. 
			if (first_time_setup) {
			    
				//Receive rows from root section and then copy the rows to the second copy array
				int recv_ammount = (rows_per_proc * dimension);
		        MPI_Recv(rows, recv_ammount, MPI_DOUBLE, root, tag_1, MPI_COMM_WORLD, &status);
		        
		        
		        
		        // Copy this into the new array as well
		        for (i = 0; i < rows_per_proc; i++) {
					for (j = 0; j < dimension; j++) {
						rows_copy[i*dimension+j] = rows[i*dimension+j];
						
					}
				}
				first_time_setup = 0;
                
            } else {
                //Send and receive the top and bottom rows for each process needed
                
				//Send topmost row
				MPI_Isend( &rows[dimension], dimension , MPI_DOUBLE, proc_rank-1, tag_1, MPI_COMM_WORLD, &request);
				
								
				//Send bottommost row unless process is the last one
				if (proc_rank != world_size - 1) {
					MPI_Isend( &rows[(rows_per_proc-2)*dimension], dimension , MPI_DOUBLE, proc_rank+1, tag_1, MPI_COMM_WORLD, &request);
					
				}
				
				
				//Receive topmost row
				MPI_Recv( &rows[0], dimension , MPI_DOUBLE, proc_rank-1, tag_1, MPI_COMM_WORLD, &status);


				// Receive bottommost row unless process is the last one
				if (proc_rank != world_size - 1) {
					MPI_Recv( &rows[(rows_per_proc-1)*dimension], dimension, MPI_DOUBLE, proc_rank+1, tag_1, MPI_COMM_WORLD, &status);
					
				}

                
            }
            
            //Array relaxation
    		in_precision = relax_array(rows, rows_copy, dimension, rows_per_proc, precision, in_precision);
    		
    		
    		
    		//Swap array pointers for next run
			double *temp = rows;
			rows = rows_copy;
            rows_copy = temp;
            
            //Send in_precision flag to root using blocking send since process must wait until it has sent this data to continue
			MPI_Send( &in_precision, 1, MPI_INT, root, tag_2, MPI_COMM_WORLD);

			//Wait to receive in_precision flag from root
			MPI_Recv( &in_precision, 1, MPI_INT, root, tag_1, MPI_COMM_WORLD, &status);

            
        }
        //If program gets here then all processes are within precision so each process must send 
        //back their rows to the root process so it can put them all back together again
        
        //Send only the critical rows
        int send_ammount = (rows_per_proc-2)*dimension;
        MPI_Send( &rows[dimension], send_ammount, MPI_DOUBLE, root, tag_2, MPI_COMM_WORLD);
		
	    //Free arrays
		free(rows);
        free(rows_copy);
    }
    
    //Finalize the MPI environment.
    MPI_Finalize();
}

//mpicc -Wall parrays_mpi.c -o parrays_mpi
//mpicc -Wall parrays_mpi.c -lm -o parrays_mpi
//mpirun -n 4 ./parrays_mpi