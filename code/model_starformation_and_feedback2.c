//
//  model_starformation_and_feedback2.c
//  
//
//  Created by Alex Rohl on 18/12/18.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "core_allvars.h"
#include "core_proto.h"


void starformation_and_feedback(int p, int centralgal, double dt, int step)
{
    //From original sfr+feedback file
    double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, SFE_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscStarsSum, DiscPre, ColdPre;
    double r_inner, r_outer;
    double reff, tdyn, cold_crit, strdotfull, H2sum; // For SFprescription==3
    
    double NewStars[N_BINS], NewStarsMetals[N_BINS];
    int i;
    
    double StarsPre = Gal[p].StellarMass;
    check_channel_stars(p);
    
    double ejected_sum = 0.0;
    reheated_mass = 0.0; // initialise
    
    // Checks that the deconstructed disc is being treated properly and not generating NaNs
    DiscGasSum = get_disc_gas(p);
    DiscStarsSum = get_disc_stars(p);
    assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
    assert(DiscStarsSum <= 1.01*Gal[p].StellarMass);
    assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
    assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
    
    f_H2_const = 1.38e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*H2FractionExponent);
    SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR; // This says if SfrEfficiency==1.0, then the true efficiency of SF from H2 is 7.75e-4 h Myr^-1
    
    // Initialise variables
    strdot = 0.0;
    strdotfull = 0.0;
    stars_sum = 0.0;
    
    Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
    Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas; // I believe TAO wants these fields.  Otherwise irrelevant
    //TAO Theoretical Astrophysical Observatory, Download catalogues of models.
    
    //recomputes fraction of gas in form of h1 h2 hydrogren
    update_HI_H2(p);
    

    //Now we want a model that can track the interaction of energy between annuli

    //will ignore 'tester' code for now
    
    //Loops through annuli
    for (i=0; i<N_BINS; i++) {
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        
        //area of annulis
        area = M_PI * (r_outer*r_outer - r_inner*r_inner);
        
        // calculation of star formation rate strdot
        if(Gal[p].Vvir>0) {
            strdot = SfrEfficiency * SFE_H2 * Gal[p].DiscH2[i];
        } else {
            strdot = 0.0;
        }
        
        stars = strdot * dt; //dt time step -> total mass of stars formed

        if(stars < MIN_STARFORMATION) //minimum threshold
            stars = 0.0;
        
        if(stars > Gal[p].DiscGas[i])
            stars = Gal[p].DiscGas[i];
