/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */


#include <hpcc.h>

int
StarStream(HPCC_Params *params) {
  int myRank, commSize;
  int rv, errCount, failure = 0, failureAll = 0;
  double copyLocalGBs, copyMinGBs, copyMaxGBs, copyAvgGBs;
  double scaleLocalGBs, scaleMinGBs, scaleMaxGBs, scaleAvgGBs;
  double addLocalGBs, addMinGBs, addMaxGBs, addAvgGBs;
  double triadLocalGBs, triadMinGBs, triadMaxGBs, triadAvgGBs;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  copyLocalGBs = copyMinGBs = copyMaxGBs = copyAvgGBs =
  scaleLocalGBs = scaleMinGBs = scaleMaxGBs = scaleAvgGBs =
  addLocalGBs = addMinGBs = addMaxGBs = addAvgGBs =
  triadLocalGBs = triadMinGBs = triadMaxGBs = triadAvgGBs = 0.0;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  rv = Stream( params, 0 == myRank, &copyLocalGBs, &scaleLocalGBs, &addLocalGBs, &triadLocalGBs,
               &failure );
  MPI_Reduce( &rv, &errCount, 1, MPI_INT, MPI_SUM, 0, comm );
  MPI_Allreduce( &failure, &failureAll, 1, MPI_INT, MPI_MAX, comm );
  if (failureAll) params->Failure = 1;

  MPI_Reduce( &copyLocalGBs, &copyMinGBs, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &copyLocalGBs, &copyAvgGBs, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &copyLocalGBs, &copyMaxGBs, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
  copyAvgGBs /= commSize;
  MPI_Reduce( &scaleLocalGBs, &scaleMinGBs, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &scaleLocalGBs, &scaleAvgGBs, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &scaleLocalGBs, &scaleMaxGBs, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
  scaleAvgGBs /= commSize;
  MPI_Reduce( &addLocalGBs, &addMinGBs, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &addLocalGBs, &addAvgGBs, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &addLocalGBs, &addMaxGBs, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
  addAvgGBs /= commSize;
  MPI_Reduce( &triadLocalGBs, &triadMinGBs, 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &triadLocalGBs, &triadAvgGBs, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
  MPI_Reduce( &triadLocalGBs, &triadMaxGBs, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
  triadAvgGBs /= commSize;

  MPI_Bcast( &copyAvgGBs, 1, MPI_DOUBLE, 0, comm ); params->StarStreamCopyGBs = copyAvgGBs;
  MPI_Bcast( &scaleAvgGBs, 1, MPI_DOUBLE, 0, comm ); params->StarStreamScaleGBs = scaleAvgGBs;
  MPI_Bcast( &addAvgGBs, 1, MPI_DOUBLE, 0, comm ); params->StarStreamAddGBs = addAvgGBs;
  MPI_Bcast( &triadAvgGBs, 1, MPI_DOUBLE, 0, comm ); params->StarStreamTriadGBs = triadAvgGBs;

  BEGIN_IO( myRank, params->outFname, outputFile);
  FPRINTF(  myRank, outputFile, "Node(s) with error %d", errCount );
  FPRINTF(  myRank, outputFile, "Minimum Copy GB/s %.6f", copyMinGBs );
  FPRINTF(  myRank, outputFile, "Average Copy GB/s %.6f", copyAvgGBs );
  FPRINTF(  myRank, outputFile, "Maximum Copy GB/s %.6f", copyMaxGBs );
  FPRINTF(  myRank, outputFile, "Minimum Scale GB/s %.6f", scaleMinGBs );
  FPRINTF(  myRank, outputFile, "Average Scale GB/s %.6f", scaleAvgGBs );
  FPRINTF(  myRank, outputFile, "Maximum Scale GB/s %.6f", scaleMaxGBs );
  FPRINTF(  myRank, outputFile, "Minimum Add GB/s %.6f", addMinGBs );
  FPRINTF(  myRank, outputFile, "Average Add GB/s %.6f", addAvgGBs );
  FPRINTF(  myRank, outputFile, "Maximum Add GB/s %.6f", addMaxGBs );
  FPRINTF(  myRank, outputFile, "Minimum Triad GB/s %.6f", triadMinGBs );
  FPRINTF(  myRank, outputFile, "Average Triad GB/s %.6f", triadAvgGBs );
  FPRINTF(  myRank, outputFile, "Maximum Triad GB/s %.6f", triadMaxGBs );
  END_IO(   myRank, outputFile );

  return 0;
}

int
SingleStream(HPCC_Params *params) {
  int myRank, commSize;
  int rv, errCount, rank, failure = 0;
  double copyLocalGBs, scaleLocalGBs, addLocalGBs, triadLocalGBs;
  double scl = 1.0 / RAND_MAX;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  copyLocalGBs = scaleLocalGBs = addLocalGBs = triadLocalGBs = 0.0;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  srand(time(NULL));
  scl *= commSize;

  /* select a node at random, but not node 0 */
  for (rank = 0; ; rank = (int)(scl * rand())) {
    if (rank > 0 && rank < commSize) break;
  }

  MPI_Bcast( &rank, 1, MPI_INT, 0, comm ); /* broadcast the rank selected on node 0 */

  if (myRank == rank) /* if this node has been selected */
    rv = Stream( params, 0 == myRank, &copyLocalGBs, &scaleLocalGBs, &addLocalGBs,
                 &triadLocalGBs, &failure );

  MPI_Bcast( &rv, 1, MPI_INT, rank, comm ); /* broadcast error code */
  MPI_Bcast( &failure, 1, MPI_INT, rank, comm ); /* broadcast failure indication */
  if (failure) params->Failure = 1;

  /* broadcast results */
  MPI_Bcast( &copyLocalGBs, 1, MPI_DOUBLE, rank, comm );
  MPI_Bcast( &scaleLocalGBs, 1, MPI_DOUBLE, rank, comm );
  MPI_Bcast( &addLocalGBs, 1, MPI_DOUBLE, rank, comm );
  MPI_Bcast( &triadLocalGBs, 1, MPI_DOUBLE, rank, comm );
  errCount = rv;
  params->SingleStreamCopyGBs = copyLocalGBs;
  params->SingleStreamScaleGBs = scaleLocalGBs;
  params->SingleStreamAddGBs = addLocalGBs;
  params->SingleStreamTriadGBs = triadLocalGBs;

  BEGIN_IO( myRank, params->outFname, outputFile);
  FPRINTF(  myRank, outputFile, "Node(s) with error %d", errCount );
  FPRINTF(  myRank, outputFile, "Node selected %d", rank );
  FPRINTF(  myRank, outputFile, "Single STREAM Copy GB/s %.6f", copyLocalGBs );
  FPRINTF(  myRank, outputFile, "Single STREAM Scale GB/s %.6f", scaleLocalGBs );
  FPRINTF(  myRank, outputFile, "Single STREAM Add GB/s %.6f", addLocalGBs );
  FPRINTF(  myRank, outputFile, "Single STREAM Triad GB/s %.6f", triadLocalGBs );
  END_IO(   myRank, outputFile );

  return 0;
}
