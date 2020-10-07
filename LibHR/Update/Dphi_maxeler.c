/*  This file contains the implementation of all the functions needed */

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
	    Dphi_flt(otmp_flt, sptr);
	    Cphi_flt(mass, otmp_flt, otmp_flt, 0);
	    apply_BCs_on_spinor_field_flt(otmp_flt);
	    Dphi_flt(dptr, otmp_flt);
	    apply_BCs_on_spinor_field_flt(dptr);
	    Cphi_flt(mass, dptr, dptr, 0);
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
