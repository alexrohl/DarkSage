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


//assume SNRecipe == 1 i.e. "scale with local density"
//assume SFprescription set to 0
void starformation_and_feedback2(int p, int centralgal, double dt, int step)
{
    //printf("feedback2 being called");
    //From original sfr+feedback file
    double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, stars_sum, area, areaIn, areaOut, SFE_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscStarsSum, DiscPre, ColdPre;
    double r_Inner, r_inner, r_outer, r_Outer;
    double reff, tdyn, cold_crit, strdotfull, H2sum; // For SFprescription==3
    
    double NewStars[N_BINS], NewStarsMetals[N_BINS];
    int i, prop;
    
    double StarsPre = Gal[p].StellarMass;
    check_channel_stars(p);
    
    double ejected_sum = 0.0;
    reheated_mass = 0.0; // initialise
    
    DiscGasSum = get_disc_gas(p);
    DiscStarsSum = get_disc_stars(p);
    
    f_H2_const = 1.38e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*H2FractionExponent);
    SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR; // This says if SfrEfficiency==1.0, then the true efficiency of SF from H2 is 7.75e-4 h Myr^-1
    
    // Initialise variables
    strdot = 0.0;
    strdotfull = 0.0;
    stars_sum = 0.0;
    
    Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
    Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas; // I believe TAO wants these fields.  Otherwise irrelevant
    //TAO Theoretical Astrophysical Observatory, Download catalogues of models.
    
    update_HI_H2(p);
    //recomputes fraction of gas in form of h1 h2 hydrogren
    
    for(i=0; i<N_BINS; i++) // Loops through annuli
    {
        if (i>0) {
            r_Inner = Gal[p].DiscRadii[i-1];
        } else {
            r_Inner = 0;
        }
        
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        
        if (i<N_BINS) {
            r_Outer = Gal[p].DiscRadii[i+2];
        } else {
            r_Outer = 0;
        }
        
        area = M_PI * (r_outer*r_outer - r_inner*r_inner); //area of annulis
        areaIn = M_PI * (r_inner*r_inner - r_Inner*r_Inner); //area of annulis
        areaOut = M_PI * (r_Outer*r_Outer - r_outer*r_outer); //area of annulis
        
        if(Gal[p].Vvir>0) //test purposes
        {
                // calculation of star formation rate
                strdot = SfrEfficiency * SFE_H2 * Gal[p].DiscH2[i];
        }
        else // These galaxies (which aren't useful for science) won't have H2 to form stars
            strdot = 0.0;
        
        stars = strdot * dt; //dt time step -> total mass of stars formed
        
        if(stars < MIN_STARFORMATION) //minimum threshold
            stars = 0.0;
        
        if(stars > Gal[p].DiscGas[i])
            stars = Gal[p].DiscGas[i];
        
        //supernoveRecipeOn>0 -> stellar feedback
        if(Gal[p].DiscGas[i] > 0.0 && stars>MIN_STARS_FOR_SN) //too low mass to make supernovas
        {
            //!!!uses equation in the paper for supernova feedback
            //finds Sigma_0gas the reference surface density
            Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / sqr(CM_PER_MPC/1e6 / UnitLength_in_cm);
            //???
            //reheated_mass = FeedbackReheatingEpsilon * stars * pow(Sigma_0gas / (Gal[p].DiscGas[i]/area), FeedbackExponent);
            
            if (i>0 && i < N_BINS-1) {
                prop = areaIn/areaOut;
                reheated_mass = FeedbackReheatingEpsilon * (NewStars[i-1] * prop + NewStars[i+1] * (1-prop));
            }
            
            // Can't use more cold gas than is available, so balance SF and feedback
            if((stars + reheated_mass) > Gal[p].DiscGas[i] && (stars + reheated_mass) > 0.0)
            {
                fac = Gal[p].DiscGas[i] / (stars + reheated_mass);
                stars *= fac;
                reheated_mass *= fac;
            }
            
            if(stars<MIN_STARS_FOR_SN)
            {
                stars = MIN_STARS_FOR_SN;
                reheated_mass = Gal[p].DiscGas[i] - stars; // Used to have (1-RecycleFraction)* in front of stars here, but changed philosophy
            }
            
            if(reheated_mass < MIN_STARFORMATION)
                reheated_mass = 0.0; // Limit doesn't have to be the same as MIN_STARFORMATION, but needs to be something reasonable
            
            //equation already formed from energy arguents
            ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) - FeedbackReheatingEpsilon) * stars;
            if(ejected_mass < MIN_STARFORMATION)
                ejected_mass = 0.0;
            
        }
        else
        {
            reheated_mass = 0.0;
            ejected_mass = 0.0;
        }
        
        Gal[p].DiscSFR[i] += stars / dt;
        //^ stores sfr
        stars_sum += stars; //counting total stars formed accross annuli
        ejected_sum += ejected_mass; //gas ejected from halo around galaxy
        
        //checking values before the function
        DiscPre = Gal[p].DiscGas[i];
        ColdPre = Gal[p].ColdGas;
        
        // Update for star formation
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        if(stars>=MIN_STARS_FOR_SN)
        {
            NewStars[i] = (1 - RecycleFraction) * stars; //builds temporary stellar disc
            //recyclefraction -> mass of stars just formed, immediately returned back to the gas
            NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
        }
        else //no feedback in this case (rare)
        {
            NewStars[i] = stars;
            NewStarsMetals[i] = metallicity * stars;
        }
        if(NewStarsMetals[i] > NewStars[i])
        {
            printf("NewStars, metals = %e, %e\n", NewStars[i], NewStarsMetals[i]);
            printf("Gas, metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        }
        update_from_star_formation(p, stars, metallicity, i);
        //^ takes mass out of the gas disc to conserve mass
        
        if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
        
        //catches bug
        if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
            reheated_mass = Gal[p].DiscGas[i];
        
        // These checks ensure numerical uncertainties don't blow up
        
        DiscPre = Gal[p].DiscGas[i];
        ColdPre = Gal[p].ColdGas;
        
        // Update from SN feedback
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        if(reheated_mass > 0.0)
            update_from_feedback(p, centralgal, reheated_mass, metallicity, i);
        //^ taking reheated gas out of disc and into the hot gas halo
        
        // Inject new metals from SN II
        if(stars>=MIN_STARS_FOR_SN)
        {
            Gal[p].DiscGasMetals[i] += Yield * stars * (1.0 - get_metallicity(NewStars[i],NewStarsMetals[i])); // Could just use metallicity variable here, surely
            Gal[p].MetalsColdGas += Yield * stars * (1.0 - get_metallicity(NewStars[i],NewStarsMetals[i]));
        }
        if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i]) printf("DiscGas, Metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
    }
    
    update_from_ejection(p, centralgal, ejected_sum);
    //^ updating masses after ejection
    
    double NewStarSum = 0.0;
    for(i=N_BINS-1; i>=0; i--) NewStarSum += NewStars[i];
    
    // Sum stellar discs together
    if(NewStarSum>0.0)
        combine_stellar_discs(p, NewStars, NewStarsMetals);
    // combines old and new stellar disc, maps stars to this new disc from old and new stellar discs
    
    // Update the star formation rate
    Gal[p].SfrFromH2[step] += stars_sum / dt;
    Gal[p].StarsFromH2 += NewStarSum;
    
    if(Gal[p].StellarMass >= MIN_STARS_FOR_SN)
    {
        check_channel_stars(p);
    }
    
    
    DiscGasSum = get_disc_gas(p);
    
}