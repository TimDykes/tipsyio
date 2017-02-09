#include <string>
#include <fstream>
#include "tipsy_file.h"
#include <mpi.h>
#include <iostream>
int main(int argc, char** argv)
{

    // // Check arguments
    // if(argc != 2)
    // {
    //     printf("Usage: ./tipsy_header infilename\nExiting...\n");
    //     exit(0);
    // }
    // int rank, size;
    // MPI_Init(&argc, &argv);
    
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);

    // std::cout << "I am rank: " << rank << " of " << size << " processes" << std::endl;

    // MPI_Finalize();
    
    // Open tipsy file, read and print header, close it again
    std::string infile(argv[1]);
    TipsyFile file(infile.c_str());
    file.read_header();
    file.report_header();
    file.close();

    return 0;
}

