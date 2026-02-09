/**
 * reconstruct.cpp
 * 
 * Contains routines for higher order reconstruction of hydrodynamic primitives.
 */
#include "aiolos.h"

double MonotonizedCentralSlope(double ql, double qm, double qr, double cF, double cB, double dxF, double dxB);
double VanLeerSlope(double ql, double qm, double qr, double cF, double cB, double dxF, double dxB);

/**
 *  General slope limiting function based on mid, left and right values.
 */
double MonotonizedCentralSlope(double ql, double qm, double qr, 
                               double cF=2, double cB=2, double dxF=1, double dxB=1) {
    // From Mignone's (2005) reconstruction paper

    double dF = (qr - qm) / dxF ;
    double dB = (qm - ql) / dxB ;

    if (dF*dB < 0)
        return 0 ;

    auto min_mod = [](double l, double r) {
        if (std::abs(l) < std::abs(r))
            return l;
        else 
            return r; 
    } ;

    return min_mod(0.5*(dF + dB), min_mod(cF*dF, cB*dB)) ;
    //return min_mod(dF, dB) ;
}

double VanLeerSlope(double ql, double qm, double qr, 
                                double cF=2, double cB=2, double dxF=1, double dxB=1) {
     double dF = (qr - qm) / dxF ;
     double dB = (qm - ql) / dxB ;
     

     if (dF*dB <= 0)
         return 0 ;

    double v  = dB/dF;

      //return dF*std::max(0., std::min( 0.5*(1.+v),  std::min(cF, v*cB) )) ;
      return dF*v*(cF*v + cB)/(v*v + v*(cF + cB -2) + 1); //Modified Van leer from Mignone+2014
}

/**
 * Reconstruct the primitive hydro states at the cell edges based on the cell centered value and a slope, which is interpreted from the neighbouring cells.
 */
void c_Species::reconstruct_edge_states( std::vector<double>&u_mask, int orderstep) {
    
    if(base->steps > base->debug_steps && this_species_index == base->debug_species  && base->debug_steps < num_cells+1) {
        cout<<"t="<<base->steps<<"IN RECONSTRUCT Pos 1 T[423]_s = [1]= "<<prim[base->debug_cell].temperature;
        cout<<" manual pressure u_in = "<<(u[base->debug_cell].u3 - 0.5*u[base->debug_cell].u2*u[base->debug_cell].u2/u[base->debug_cell].u1)<<endl;
    }

    const bool is_gas = !is_dust_like ;
    int problematic = 0;

    // Step 1:
    //   1st order hydrostatic reconstruction
    for (int i=0; i < num_cells+2; i++) {
        prim_r[i] = prim_l[i] = prim[i] ;

        if (is_gas) {
            if(i > 0)
                prim_l[i].pres = prim[i].pres + prim[i].density * get_dp_hydrostatic(i,  -1);
            if(i <= num_cells)
                prim_r[i].pres = prim[i].pres + prim[i].density * get_dp_hydrostatic(i,  +1);

            // TODO: add warning?
            if (prim_l[i].pres < 0 || prim_r[i].pres < 0) {
                prim_l[i] = prim_r[i] = prim[i] ;
                problematic += 1;
            }
        }
    }

    // Step 2: Add 2nd order slope-limited correction
    IntegrationType order = base->order ;
    //if (order == IntegrationType::second_order && orderstep == 1) {
    if (orderstep == 1) {
        const std::vector<double>& 
            x_i = base->x_i, 
            x_iVC = base->x_iVC,
            dx = base->dx ;
       
        for (int i=1; i <= num_cells; i++) {
            //double maskmul = u_mask[i]>0.5?0.:1.; //i>num_cells/2?0.:1.; // Multiply all slopes with 1 in the nominal case, or 0 in case we get a message from above that this cell is broken
            //maskmul = 1.;
            double maskmul_l = 1.;
            double maskmul_r = 1.;

            double p_ratio_l = std::log10(prim[i-1].pres/prim[i].pres);
            double p_ratio_r = std::log10(prim[i].pres/prim[i+1].pres);

            int limratio = 1;
            if(p_ratio_l < -limratio || p_ratio_l > limratio)
            maskmul_l = 1.;
            if(p_ratio_r < -limratio || p_ratio_r > limratio)
            maskmul_r = 1.;

            double dp_l = 0, dp_r = 0 ;
            if (is_gas) {
                dp_l = 
                    prim[ i ].density * get_dp_hydrostatic( i ,  -1) -
                    prim[i-1].density * get_dp_hydrostatic(i-1,  +1) ;
                dp_r = 
                    prim[ i ].density * get_dp_hydrostatic( i ,  +1) -
                    prim[i+1].density * get_dp_hydrostatic(i+1,  -1) ;
            }
        
            // Non-uniform grid factors
            double cF = (x_iVC[i+1] - x_iVC[i]) / (x_i[i] - x_iVC[i]) ;
            double cB = (x_iVC[i] - x_iVC[i-1]) / (x_iVC[i] - x_i[i-1]) ;

            double dxF = (x_iVC[i+1] - x_iVC[i]) ;
            double dxB = (x_iVC[i] - x_iVC[i-1]) ;

            // Pressure perturbation
            double slope ;
            slope = reconstruct_pointer(
                prim[i-1].pres -  dp_l, prim[i].pres, prim[i+1].pres - dp_r, cF, cB, dxF, dxB) ;

            //prim_l[i].pres += maskmul_l * slope * (x_i[i-1] - x_iVC[i]) ; 
            //prim_r[i].pres += maskmul_r * slope * (x_i[ i ] - x_iVC[i]) ;
            prim_l[i].pres += slope * (x_i[i-1] - x_iVC[i]) ; 
            prim_r[i].pres += slope * (x_i[ i ] - x_iVC[i]) ;
            

            // Density  //Changed the sloping to be in number densities in order to investigate the relative numerical particle drift (see appendix of Helium paper)
                        //Result: that's not what's causing it, maybe its better to reinstate the sloping in mass density?
             slope = reconstruct_pointer(
                 prim[i-1].density, prim[i].density, prim[i+1].density, cF, cB, dxF, dxB) ;
            // 
             prim_l[i].density +=  slope * (x_i[i-1] - x_iVC[i]) ; 
             prim_r[i].density +=  slope * (x_i[ i ] - x_iVC[i]) ;
            //            
            // Speed
            slope = reconstruct_pointer(
                prim[i-1].speed, prim[i].speed, prim[i+1].speed, cF, cB, dxF, dxB) ;

            //prim_l[i].speed +=  maskmul_l *slope * (x_i[i-1] - x_iVC[i]) ; 
            //prim_r[i].speed +=  maskmul_r *slope * (x_i[ i ] - x_iVC[i]) ;
            prim_l[i].speed +=  slope * (x_i[i-1] - x_iVC[i]) ; 
            prim_r[i].speed +=  slope * (x_i[ i ] - x_iVC[i]) ;
            
             if ((prim_l[i].pres < 0) || (prim_r[i].pres < 0) || 
                 (prim_l[i].density < 0) || (prim_r[i].density < 0)) {
                     prim_l[i] = prim_r[i] = prim[i] ;
                     problematic += 1;
                 }
        }
    }

//
// 2025-01-13: Additional corrections for weird pressure gradients
//
	for (int i=1; i <= num_cells; i++) {
		//Detect under/overshoots
		for(int sw=0; sw <= 1; sw++) {
			if( (prim[i+sw].pres > prim[i+(1-sw)].pres) && (prim_r[i+sw].pres < prim_l[i+(1-sw)].pres ) ) {
				double pmean = 0.5*(prim_r[i].pres + prim_l[i+1].pres);
				prim_r[i].pres   = pmean;
				prim_l[i+1].pres = pmean;
			}
		}

		for(int sw=0; sw <= 1; sw++) {
                        if( (prim[i+sw].density > prim[i+(1-sw)].density) && (prim_r[i+sw].density < prim_l[i+(1-sw)].density ) ) {
                                double dmean = 0.5*(prim_r[i].density + prim_l[i+1].density);
                                prim_r[i].density   = dmean;
                                prim_l[i+1].density = dmean;
                        }
                }

		//Detect undershoot if dp/dr > 0
	}


    // Step 3:
    //   Extra thermodynamic variables
    eos->compute_auxillary(&(prim_l[0]), num_cells+2);
    eos->compute_auxillary(&(prim_r[0]), num_cells+2);
    
    if(base->steps > base->debug_steps && this_species_index == base->debug_species) {
        for(int k=-1;k<=1;k++) {
            cout<<"t="<<base->steps<<"IN RECONSTRUCT Pos -1 T["<<(base->debug_cell+k)<<"]_s = [1]= "<<prim[base->debug_cell+k].temperature;
            cout<<" manual pressure u_in = "<<(u[base->debug_cell+k].u3 - 0.5*u[base->debug_cell+k].u2*u[base->debug_cell+k].u2/u[base->debug_cell+k].u1);
            cout<<" problematic = "<<problematic<<" p_lr = "<<prim_l[base->debug_cell+k].pres<<"/"<<prim_r[base->debug_cell+k].pres;
            cout<<" v_lr = "<<prim_l[base->debug_cell+k].speed<<"/"<<prim_r[base->debug_cell+k].speed;
            cout<<" rho_lr = "<<prim_l[base->debug_cell+k].density<<"/"<<prim_r[base->debug_cell+k].density<<endl;
            
        }

    }

} 