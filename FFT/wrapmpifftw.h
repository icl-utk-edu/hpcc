

#ifdef USING_FFTW
#include <fftw_mpi.h>
#else
#include <mpi.h>
typedef struct hpcc_fftw_mpi_plan_struct *fftw_mpi_plan;
#define fftw_mpi_create_plan  hpcc_fftw_mpi_create_plan
#define fftw_mpi_destroy_plan hpcc_fftw_mpi_destroy_plan
#define fftw_mpi hpcc_fftw_mpi
#define fftw_mpi_local_sizes hpcc_fftw_mpi_local_sizes
#endif

struct hpcc_fftw_mpi_plan_struct {
  MPI_Comm comm;
  fftw_complex *wx, *wy, *wz, *c, *d, *ww, *www;
  int n0, n;
  int flags;
  fftw_direction dir;
};
typedef struct hpcc_fftw_mpi_plan_struct *hpcc_fftw_mpi_plan;

extern hpcc_fftw_mpi_plan
hpcc_fftw_mpi_create_plan(MPI_Comm comm, int n, fftw_direction dir, int flags);
extern void hpcc_fftw_mpi_destroy_plan(hpcc_fftw_mpi_plan plan);
extern void hpcc_fftw_mpi(hpcc_fftw_mpi_plan p, int n_fields, fftw_complex *local_data,
		     fftw_complex *work);
extern void hpcc_fftw_mpi_local_sizes(hpcc_fftw_mpi_plan p, int *local_n, int *local_start,
                          int *local_n_after_transform, int *local_start_after_transform,
                          int *total_local_size);
