#ifndef __MPI_Tail_Work_hpp_
#define __MPI_Tail_Work_hpp_

#include <functional>
#include <mpi.h>

namespace DaMaSCUS_SUN
{

// Rendezvous at the end of an MPI work round without idling early ranks.
// Every rank must call this function in the same collective order. While the
// nonblocking barrier is waiting for a straggler, ranks that still have local
// capacity execute independent work one item at a time.
unsigned long int Perform_MPI_Tail_Work(
	MPI_Comm communicator,
	unsigned long int maximum_extra_work,
	const std::function<void()>& perform_one_work_item);

} // namespace DaMaSCUS_SUN

#endif
