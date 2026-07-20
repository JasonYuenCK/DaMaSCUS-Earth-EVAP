#include "MPI_Tail_Work.hpp"

#include <stdexcept>

namespace DaMaSCUS_SUN
{

unsigned long int Perform_MPI_Tail_Work(
	MPI_Comm communicator,
	unsigned long int maximum_extra_work,
	const std::function<void()>& perform_one_work_item)
{
	if(!perform_one_work_item)
		throw std::invalid_argument("Perform_MPI_Tail_Work requires a work callback");

	MPI_Request rendezvous = MPI_REQUEST_NULL;
	MPI_Ibarrier(communicator, &rendezvous);

	unsigned long int completed_extra_work = 0;
	int all_ranks_ready = 0;
	MPI_Test(&rendezvous, &all_ranks_ready, MPI_STATUS_IGNORE);
	while(!all_ranks_ready && completed_extra_work < maximum_extra_work)
	{
		perform_one_work_item();
		completed_extra_work++;
		MPI_Test(&rendezvous, &all_ranks_ready, MPI_STATUS_IGNORE);
	}

	if(!all_ranks_ready)
		MPI_Wait(&rendezvous, MPI_STATUS_IGNORE);
	return completed_extra_work;
}

} // namespace DaMaSCUS_SUN
