/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * See RandomAccess.h for a comprehensive description of this test and
 * its goals.
 *
 * This file contains the interface for the star single cpu RandomAccess test.
 * The test runs on all cpus in the MPI universe, but there is no communication
 * between cpus during the process (each cpu runs its own version of the
 * single_cpu test).  The final result is the average of the entire system.
 */
/* This file contains the interface for the single cpu RandomAccess test.  The
 * test is only run on a single (random) node in the MPI universe, with all
 * other CPUs stuck (in theory, idle) in an MPI_Bcast waiting for the selected
 * CPU to finish the RandomAccess test. */

#include <hpcc.h>
#include "RandomAccess.h"

int
StarRandomAccess(HPCC_Params *params)
{
  int myRank, commSize;
  int rv, errCount, failure = 0, failureAll = 0;
  double minGUPs, avgGUPs, maxGUPs, localGUPs;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  minGUPs = avgGUPs = maxGUPs = localGUPs = 0.0;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  rv = RandomAccess( params, 0 == myRank, &localGUPs, &failure );
  MPI_Reduce( &rv, &errCount, 1, MPI_INT, MPI_SUM, 0, comm );
  MPI_Allreduce( &failure, &failureAll, 1, MPI_INT, MPI_MAX, comm );
  if (failureAll) params->Failure = 1;

  MPI_Reduce( &localGUPs, &minGUPs, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &localGUPs, &avgGUPs, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &localGUPs, &maxGUPs, 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  avgGUPs /= commSize;

  MPI_Bcast( &avgGUPs, 1, MPI_DOUBLE, 0, comm );
  params->StarGUPs = avgGUPs;

  BEGIN_IO( myRank, params->outFname, outputFile);
  FPRINTF(  myRank, outputFile, "Node(s) with error %d", errCount );
  FPRINTF(  myRank, outputFile, "Minimum GUP/s %.6f", minGUPs );
  FPRINTF(  myRank, outputFile, "Average GUP/s %.6f", avgGUPs );
  FPRINTF(  myRank, outputFile, "Maximum GUP/s %.6f", maxGUPs );
  END_IO(   myRank, outputFile );

  return 0;
}

int
SingleRandomAccess(HPCC_Params *params)
{
  int myRank, commSize;
  int rv, errCount, rank, failure = 0;
  double localGUPs;
  double scl = 1.0 / RAND_MAX;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  localGUPs = 0.0;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  srand(time(NULL));
  scl *= commSize;

  /* select a node at random, but not node 0 (unless there is just one node) */
  if (1 == commSize)
    rank = 0;
  else
    for (rank = 0; ; rank = (int)(scl * rand())) {
      if (rank > 0 && rank < commSize) break;
    }

  MPI_Bcast( &rank, 1, MPI_INT, 0, comm ); /* broadcast the rank selected on node 0 */

  if (myRank == rank) /* if this node has been selected */
    rv = RandomAccess( params, 0 == myRank, &localGUPs, &failure );

  MPI_Bcast( &rv, 1, MPI_INT, rank, comm ); /* broadcast error code */
  MPI_Bcast( &localGUPs, 1, MPI_DOUBLE, rank, comm ); /* broadcast GUPs */
  MPI_Bcast( &failure, 1, MPI_INT, rank, comm ); /* broadcast failure indication */
  errCount = rv;
  params->SingleGUPs = localGUPs;
  if (failure) params->Failure = 1;

  BEGIN_IO( myRank, params->outFname, outputFile);
  FPRINTF(  myRank, outputFile, "Node(s) with error %d", errCount );
  FPRINTF(  myRank, outputFile, "Node selected %d", rank );
  FPRINTF(  myRank, outputFile, "Single GUP/s %.6f", localGUPs );
  END_IO(   myRank, outputFile );

  return 0;
}
