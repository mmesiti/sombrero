#ifndef INCLUDE_MAX_CG_API_H_
#define INCLUDE_MAX_CG_API_H_
#include <stddef.h>

/*
 * Struct representing the 3x3 SU(3) matrices that make up a D-slash operator
 */
typedef struct {
   _Complex float c00, c01, c10, c11, c20, c21;
} su3;

/*
 *  Struct representing the spinors that live on the lattice nodes
 */
typedef struct {
   _Complex float s0c0, s1c0, s2c0, s3c0;
   _Complex float s0c1, s1c1, s2c1, s3c1;
   _Complex float s0c2, s1c2, s2c2, s3c2;
} cg_spinor;

/*
 * Struct representing the clovers that make up a clover operator
 */
typedef struct {
   _Complex float c0[18];
   _Complex float c1[18];
} cg_clover;

/*
 * Struct for specifying lattice dimensions
 */
typedef struct {
	int LX, LY, LZ, LT;
} volume_spec;



/*
 * Run conjugate gradient to solution, CPU model version
 * The matrix-vector equation A*x = b will be solved for x. The input x is a trial solution that will be overridden with the
 * CG solution. The input b is a known vector, and A is the matrix formed from the input gauges, the clovers, and the
 * parameter gamma (see documentation for more details on the matrix structure).
 */
 void max_cg_cpu_model(
		cg_spinor *in_x, cg_spinor *in_b, su3 *gauge_u01, cg_clover *clover0, cg_clover *clover1,
		double gamma, double bb,
		int *niter, int *no_convergence, int max_iter, double res,
		int LX, int LY, int LZ, int LT
);

/**
 * Compute A*in_vec, where A is the matrix formed from the input gauges, the clovers, and the parameter gamma (see
 * documentation for more details on the matrix structure).
 */
void max_cg_mat_mult_cpu(
		cg_spinor *result, cg_spinor *in_vec, su3 *gauge_u01, cg_clover *clover0, cg_clover *clover1,
		double gamma,
		int LX, int LY, int LZ, int LT
);


/**
 * Interleave gauges, such that two input vectors of gauge quartets become a single vector of gauge octets.
 */
void max_cg_interleave_gauges_float(su3 *out, su3 *u0, su3 *u1, size_t volumeH);

/**
 * De-interleave gauges, such that one input vector of gauge octets becomes two vectors of gauge quartets.
 */
void max_cg_deinterleave_gauges_float(su3 *out0, su3 *out1, su3 *in, size_t volumeH);

/**
 * Reorder data from one subdomain size to another. The input data should contain data of element_size_bytes bytes per lattice
 * point, subdomain-by-subdomain (i.e. first all data for subdomain 0, then all data for subdomain 1, etc). The output data
 * will contain the same data but reordered according to the new subdomain size.
 */
void max_cg_reorder_subdomains(void *result, void *in_data, size_t element_size_bytes,
	size_t subdomain_LX_in, size_t subdomain_LY_in, size_t subdomain_LZ_in, size_t subdomain_LT_in,
	size_t subdomain_LX_out, size_t subdomain_LY_out, size_t subdomain_LZ_out, size_t subdomain_LT_out,
	size_t domain_LX, size_t domain_LY, size_t domain_LZ, size_t domain_LT
);



/**
 * Read a vector of spinors from disk, in a specified format. Each spinor is additionally expected to have a one-line header.
 */
void max_cg_read_spinor(cg_spinor *result, int *maximum_exponent, char *filename, char *format, volume_spec vol);

/**
 * Read a vector of gauges from disk, in the specified format. Each gauge is additionally expected to have a one-line
 * header. The number of gauges per lattice point can be specified. The maximum exponent of the gauge is determined simultaneously
 * with reading it.
 */
void max_cg_read_gauge(
		su3 *result, int *max_exp, char *filename, char *format,
		size_t num_gauges_per_lattice_point, volume_spec vol
);

/**
 * Read a vector of clovers from disk, in the specified format. Each clover is additionally expected to have a one-line
 * header. The maximum exponents of the clover are also determined simulataneously with reading it.
 */
void max_cg_read_clover(
		cg_clover *result, int *max_exp_diagonal, int *max_exp_offdiagonal,
		char *filename, char *format, volume_spec vol
);

/**
 * Write a vector of spinors to disk using a specified format
 */
void max_cg_write_spinor(char *filename, char *format, cg_spinor *to_write, volume_spec vol);

/**
 * Write a vector of gauges to disk using the specified format
 */
void max_cg_write_gauge(
	char *filename, char *format, su3 *s, size_t num_gauges_per_lattice_point, volume_spec vol
);

/**
 * Write a vector of clovers to disk using a specified format
 */
void max_cg_write_clover(char *filename, char *format, cg_clover *to_write, volume_spec vol);


// DEBUG: inferred from examples
typedef enum {SOLVE, RESIDUE, NEG_MAT_MULT} cg_mode;


cg_spinor* max_cg_mpi(cg_spinor* in_x, cg_spinor* in_b,
                      su3* gauge_u01, cg_clover *clover0,  cg_clover *clover1,
                      int global_max_exp_gauge, int global_max_exp_clover_diag,
                      int global_max_exp_clover_offdiag, double gamma, double bb,
                      int *niter, int *no_convergence, int max_iter, double res,
                      int LX_domain, int LY_domain,
                      int LZ_domain, int LT_domain,
                      int LX_mpi_subdomain, int LY_mpi_subdomain,
                      int LZ_mpi_subdomain, int LT_mpi_subdomain,
                      int my_pe, int n_pes, cg_mode mode);


#endif /* INCLUDE_MAX_CG_API_H_ */
