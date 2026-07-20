#include <mpi.h>

#include <chrono>
#include <iostream>
#include <thread>

#include "MPI_Tail_Work.hpp"

using namespace DaMaSCUS_SUN;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank = 0;
	int processes = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
	if(processes != 2)
	{
		if(rank == 0)
			std::cerr << "test_MPI_Tail_Work requires exactly two ranks" << std::endl;
		MPI_Finalize();
		return 1;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
		std::this_thread::sleep_for(std::chrono::milliseconds(200));

	unsigned long int callback_count = 0;
	const unsigned long int completed = Perform_MPI_Tail_Work(
		MPI_COMM_WORLD,
		1000UL,
		[&]()
		{
			callback_count++;
			std::this_thread::sleep_for(std::chrono::milliseconds(2));
		});

	int local_failure = (completed != callback_count) ? 1 : 0;
	if(rank == 1 && completed == 0)
	{
		std::cerr << "early rank performed no useful work while waiting for the delayed rank" << std::endl;
		local_failure = 1;
	}

	int global_failure = 0;
	MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Finalize();
	return global_failure;
}
