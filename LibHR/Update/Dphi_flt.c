/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File Dphi_flt.c
*
* Action of the Wilson-Dirac operator D on a given single-precision spinor field
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "error.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "spinor_field.h"
#include "geometry.h"
#include "communications.h"
#include "memory.h"

#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */


/*
 * Init of Dphi_flt
 */

static int init_dirac=1;
static spinor_field_flt *gtmp=NULL;
static spinor_field_flt *etmp=NULL;
static spinor_field_flt *otmp=NULL;

static void free_mem() {
    if (gtmp!=NULL) { free_spinor_field_f_flt(gtmp); etmp=NULL; }
    if (etmp!=NULL) { free_spinor_field_f_flt(etmp); etmp=NULL; }
    if (otmp!=NULL) { free_spinor_field_f_flt(otmp); otmp=NULL; }
    init_dirac=1;
}

static void init_Dirac() {
    if (init_dirac) {
        gtmp=alloc_spinor_field_f_flt(1,&glattice);
        etmp=alloc_spinor_field_f_flt(1,&glat_even);
        otmp=alloc_spinor_field_f_flt(1,&glat_odd);
        atexit(&free_mem);
        init_dirac=0;
    }
}


/*
 * the following variable is used to keep trace of
 * matrix-vector multiplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

unsigned long int getMVM_flt() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	MVMcounter=0; /* reset counter */

	return res;
}

/* Theta Boundary conditions 
 * local copy in single precision of global variable
 */
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
static complex_flt eitheta_flt[4]={{1.f,0.f}};
#endif

/* r=t*u*s */
#ifdef BC_T_THETA

#define _suNf_theta_T_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta_flt[0],vtmp)

#define _suNf_theta_T_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta_flt[0],vtmp)

#else

#define _suNf_theta_T_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_T_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_X_THETA

#define _suNf_theta_X_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta_flt[1],vtmp)

#define _suNf_theta_X_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta_flt[1],vtmp)

#else

#define _suNf_theta_X_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_X_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Y_THETA

#define _suNf_theta_Y_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta_flt[2],vtmp)

#define _suNf_theta_Y_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta_flt[2],vtmp)

#else

#define _suNf_theta_Y_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_Y_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Z_THETA

#define _suNf_theta_Z_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta_flt[3],vtmp)

#define _suNf_theta_Z_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta_flt[3],vtmp)

#else

#define _suNf_theta_Z_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_Z_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif




/*
 * NOTE :
 * here we are making the assumption that the geometry is such that
 * all even sites are in the range [0,VOLUME/2[ and all odd sites are
 * in the range [VOLUME/2,VOLUME[
 */
void Dphi_flt_(spinor_field_flt *out, spinor_field_flt *in)
{

   error((in==NULL)||(out==NULL),1,"Dphi_flt_ [Dphi_flt.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"Dphi_flt_ [Dphi_flt.c]",
         "Input and output fields must be different");

#ifndef CHECK_SPINOR_MATCHING
   error(out->type==&glat_even && in->type!=&glat_odd,1,"Dphi_ [Dphi.c]", "Spinors don't match! (1)");
   error(out->type==&glat_odd && in->type!=&glat_even,1,"Dphi_ [Dphi.c]", "Spinors don't match! (2)");
   error(out->type==&glattice && in->type!=&glattice,1,"Dphi_ [Dphi.c]", "Spinors don't match! (3)");
#endif


   ++MVMcounter; /* count matrix calls */
   if(out->type==&glattice) ++MVMcounter;

/************************ loop over all lattice sites *************************/
   /* start communication of input spinor field */
   start_sf_sendrecv_flt(in);

  _PIECE_FOR(out->type,ixp) {
     if(ixp==out->type->inner_master_pieces) {
       _OMP_PRAGMA( master )
       /* wait for spinor to be transfered */
       complete_sf_sendrecv_flt(in);
       _OMP_PRAGMA( barrier )
     }
     _SITE_FOR(out->type,ixp,ix) {

       int iy;
       suNf_flt *up,*um;
       suNf_vector_flt psi,chi;
       suNf_spinor_flt *r=0,*sp,*sm;
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
       suNf_vector_flt vtmp;
#endif

       r=_FIELD_AT(out,ix);

       /******************************* direction +0 *********************************/

       iy=iup(ix,0);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_flt(ix,0);

       _vector_add_f(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_T_multiply(chi,(*up),psi);

       (*r).c[0]=chi;
       (*r).c[2]=chi;

       _vector_add_f(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_T_multiply(chi,(*up),psi);

       (*r).c[1]=chi;
       (*r).c[3]=chi;

       /******************************* direction -0 *********************************/

       iy=idn(ix,0);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_flt(iy,0);

       _vector_sub_f(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_sub_assign_f((*r).c[2],chi);

       _vector_sub_f(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_sub_assign_f((*r).c[3],chi);

       /******************************* direction +1 *********************************/

       iy=iup(ix,1);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_flt(ix,1);

       _vector_i_add_f(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_X_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_sub_assign_f((*r).c[3],chi);

       _vector_i_add_f(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_X_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_sub_assign_f((*r).c[2],chi);

       /******************************* direction -1 *********************************/

       iy=idn(ix,1);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_flt(iy,1);

       _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_add_assign_f((*r).c[3],chi);

       _vector_i_sub_f(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_add_assign_f((*r).c[2],chi);

       /******************************* direction +2 *********************************/

       iy=iup(ix,2);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_flt(ix,2);

       _vector_add_f(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_Y_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_add_assign_f((*r).c[3],chi);

       _vector_sub_f(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_Y_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_sub_assign_f((*r).c[2],chi);

       /******************************* direction -2 *********************************/

       iy=idn(ix,2);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_flt(iy,2);

       _vector_sub_f(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_sub_assign_f((*r).c[3],chi);

       _vector_add_f(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_add_assign_f((*r).c[2],chi);

       /******************************* direction +3 *********************************/

       iy=iup(ix,3);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_flt(ix,3);

       _vector_i_add_f(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_sub_assign_f((*r).c[2],chi);

       _vector_i_sub_f(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_add_assign_f((*r).c[3],chi);

       /******************************* direction -3 *********************************/

       iy=idn(ix,3);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_flt(iy,3);

       _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_add_assign_f((*r).c[2],chi);

       _vector_i_add_f(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_sub_assign_f((*r).c[3],chi);

       /******************************** end of loop *********************************/

       _spinor_mul_f(*r,-0.5f,*r);

     } /* SITE_FOR */
   } /* PIECE FOR */
}


#ifdef WITH_CLOVER

/*************************************************
 * Dirac operators with clover term:             *
 * Cphi = Dphi + clover                          *
 * Cphi_eopre = D_ee - D_eo D_oo^-1 D_oe         *
 * Cphi_diag = D_oo or D_ee                      *
 * Cphi_diag_inv = D_oo^-1 or D_ee^-1            *
 *************************************************/

#if defined(REPR_ADJOINT)
typedef suNfc_flt clover_type_flt;
#define clover_inverse_multiply(a,b,c) _suNfc_inverse_multiply(a,b,c) 
#define clover_multiply(a,b,c) _suNfc_multiply(a,b,c) 
#else
typedef suNffull_flt clover_type_flt;
#define clover_inverse_multiply(a,b,c) _suNffull_inverse_multiply(a,b,c) 
#define clover_multiply(a,b,c) _suNffull_multiply(a,b,c) 
#endif

static void Cphi_flt_(float mass, spinor_field_flt *dptr, spinor_field_flt *sptr)
{
	// Correct mass term
	mass = (4.f+mass);

	// Loop over local sites
	_MASTER_FOR(dptr->type,ix)
	{
		suNf_vector_flt v1, v2;
		suNf_spinor_flt *out, *in, tmp;
		clover_type_flt *s0, *s1, *s2, *s3;

		// Field pointers
		out = _FIELD_AT(dptr,ix);
		in = _FIELD_AT(sptr,ix);
		s0 = _4FIELD_AT(cl_term_flt,ix,0);
		s1 = _4FIELD_AT(cl_term_flt,ix,1);
		s2 = _4FIELD_AT(cl_term_flt,ix,2);
		s3 = _4FIELD_AT(cl_term_flt,ix,3);

		// Component 0
		clover_multiply(v1, *s0, in->c[0]);
		clover_multiply(v2, *s1, in->c[1]);
		_vector_add_f(tmp.c[0], v1, v2);

		// Component 1
		clover_inverse_multiply(v1, *s1, in->c[0]);
		clover_multiply(v2, *s0, in->c[1]);
		_vector_sub_f(tmp.c[1], v1, v2);

		// Component 2
		clover_multiply(v1, *s2, in->c[2]);
		clover_multiply(v2, *s3, in->c[3]);
		_vector_add_f(tmp.c[2], v1, v2);

		// Component 3
		clover_inverse_multiply(v1, *s3, in->c[2]);
		clover_multiply(v2, *s2, in->c[3]);
		_vector_sub_f(tmp.c[3], v1, v2);

		// Add mass
		_spinor_mul_add_assign_f(tmp, mass, *in);

		// Store
        *out = tmp;
	}
}
#endif

#ifdef WITH_CLOVER


void maxeler_fake_eopre_flt(float mass, spinor_field_flt *dptr, spinor_field_flt *sptr)
{
	if(init_dirac)
	{
		init_Dirac();
	}
	apply_BCs_on_spinor_field_flt(sptr);
    // 1 - D_{ee}^{-1} D_{eo} D_{oo}^{-1} D_{oe}
    {
        // D_{ee}^{-1} D_{eo} D_{oo}^{-1} D_{oe}
	    Dphi_flt_(otmp, sptr);
	    Cphi_flt_(mass, otmp, otmp);
	    apply_BCs_on_spinor_field_flt(otmp);
	    Dphi_flt_(dptr, otmp);
	    apply_BCs_on_spinor_field_flt(dptr);
	    Cphi_flt_(mass, dptr, dptr);
    }

    // 1 - (...)
	spinor_field_sub_f_flt(dptr, sptr, dptr);
}

void maxeler_fake_eopre_sq_flt(float mass, spinor_field_flt *dptr, spinor_field_flt *sptr)
{
	if(init_dirac)
	{
		init_Dirac();
	}

	maxeler_fake_eopre_flt(mass, etmp, sptr);
	maxeler_fake_eopre_flt(mass, dptr, etmp);
}


#endif //#ifdef WITH_CLOVER
