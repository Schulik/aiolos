
#include "aiolos.h"

double MonotonizedCentralSlope(double ql, double qm, double qr, 
                                      double cF=2, double cB=2) {
    // From Mignone's reconstruction paper

    double dF = (qr - qm);
    double dB = (qm - ql);

    if (dF*dB < 0)
        return 0 ;

    auto min_mod = [](double l, double r) {
        if (std::abs(l) < std::abs(r))
            return l;
        else 
            return r; 
    } ;

    return min_mod(0.5*(dF + dB), min_mod(cF*dF, cB*dB)) ;
}


void c_Species::reconstruct_edge_states() {

    // Step 1:
    //   1st order hydrostatic reconstruction
    for (int i=0; i <= num_cells+1; i++) {
        prim_r[i] = prim_l[i] = prim[i] ;

        if(i > 0)
            prim_l[i].pres = prim[i].pres + prim[i].density * get_dp_hydrostatic(i,  -1);
        if(i <= num_cells)
            prim_r[i].pres = prim[i].pres + prim[i].density * get_dp_hydrostatic(i,  +1);

        // TODO: add warning?
        if (prim_l[i].pres < 0 || prim_r[i].pres < 0) {
            prim_l[i] = prim_r[i] = prim[i] ;
        }
    }

    // Step 2: Add 2nd order slope-limited correction
    IntegrationType order = base->order ;
    if (order == IntegrationType::second_order) {
        const std::vector<double>& 
            x_i = base->x_i, 
            x_iVC = base->x_iVC,
            dx = base->dx ;
       
        for (int i=1; i <= num_cells; i++) {
            double dp_l = 
                prim[ i ].density * get_dp_hydrostatic( i ,  -1) -
                prim[i-1].density * get_dp_hydrostatic(i-1,  +1) ;
            double dp_r = 
                prim[ i ].density * get_dp_hydrostatic( i ,  +1) -
                prim[i+1].density * get_dp_hydrostatic(i+1,  -1) ;
        
            // Non-uniform grid factors
            double cF = (x_iVC[i+1] - x_iVC[i]) / (x_i[i] - x_iVC[i]) ;
            double cB = (x_iVC[i] - x_iVC[i-1]) / (x_iVC[i] - x_i[i-1]) ;

            // Pressure perturbation
            double slope ;
            slope = MonotonizedCentralSlope(
                prim[i-1].pres - dp_l, prim[i].pres, prim[i+1].pres - dp_r, cF, cB) ;

            prim_l[i].pres += slope * (x_i[i-1] - x_iVC[i]) / dx[i] ; 
            prim_r[i].pres += slope * (x_i[ i ] - x_iVC[i]) / dx[i] ;


            // Density
            slope = MonotonizedCentralSlope(
                prim[i-1].density, prim[i].density, prim[i+1].density, cF, cB) ;

            prim_l[i].density += slope * (x_i[i-1] - x_iVC[i]) / dx[i] ; 
            prim_r[i].density += slope * (x_i[ i ] - x_iVC[i]) / dx[i] ;

            // Speed
            slope = MonotonizedCentralSlope(
                prim[i-1].speed, prim[i].speed, prim[i+1].speed, cF, cB) ;

            prim_l[i].speed += slope * (x_i[i-1] - x_iVC[i]) / dx[i] ; 
            prim_r[i].speed += slope * (x_i[ i ] - x_iVC[i]) / dx[i] ;
            
            if ((prim_l[i].pres < 0) || (prim_r[i].pres < 0) || 
                (prim_l[i].density < 0) || (prim_r[i].density < 0))
                prim_l[i] = prim_r[i] = prim[i] ;

        }
    }

    // Step 3:
    //   Extra thermodynamic variables
    eos->compute_auxillary(&(prim_l[0]), num_cells+2);
    eos->compute_auxillary(&(prim_r[0]), num_cells+2);

} 
