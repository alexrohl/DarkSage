#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "core_allvars.h"
#include "core_proto.h"

//new paramters
//type: -1 : disk instability
//       0 : mergers
//       1 : starform + feedback


    

void feedback(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double strdotfull)
{
    //printf("modular function running\n");
    double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, area, SFE_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscStarsSum, DiscPre, ColdPre;
    double r_inner, r_outer;
    double reff, tdyn, cold_crit, H2sum; // For SFprescription==3
    
    int i;
    
    double StarsPre = Gal[p].StellarMass;
    check_channel_stars(p);
    
    double ejected_sum = 0.0;
    reheated_mass = 0.0; // initialise
    strdot = 0.0;
    
    // Loops through annuli
    for(i=0; i<N_BINS; i++)
        
    {
        if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
        assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
        
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        
        area = M_PI * (r_outer*r_outer - r_inner*r_inner); //area of annulis
        
        
        if(Gal[p].Vvir>0) //test purposes
        {
            if(SFprescription==1 && Gal[p].DiscH2[i]<0.5*Gal[p].DiscHI[i])
            {
                double bb = sqrt(Gal[p].DiscStars[i]*Gal[p].DiscStars[0])/area; // quadratic b term
                double cc = -pow(0.5/f_H2_const, 1.0/0.92);
                double Sig_gas_half = 0.5*(-bb + sqrt(bb*bb-4.0*cc));
                double SFE_gas = SFE_H2 * 0.75 * 1.0/(1.0/0.5 + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i])/1.3 / Sig_gas_half;
                strdot = SfrEfficiency * SFE_gas * sqr(Gal[p].DiscGas[i]) / area;
            }
            else if(SFprescription==2)
            {
                if(Gal[p].ColdGas>0.0)
                    strdot = strdotfull * Gal[p].DiscGas[i] / Gal[p].ColdGas;// * Gal[p].DiscH2[i] / H2sum;
                else
                    strdot = 0.0;
            }
            else
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
        if(SupernovaRecipeOn > 0 && Gal[p].DiscGas[i] > 0.0 && stars>MIN_STARS_FOR_SN) //too low mass to make supernovas
        {
            if(SupernovaRecipeOn == 1)
                //!!!uses equation in the paper for supernova feedback
            {
                //finds Sigma_0gas the reference surface density
                Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / sqr(CM_PER_MPC/1e6 / UnitLength_in_cm);
                
                reheated_mass = FeedbackReheatingEpsilon * stars * pow(Sigma_0gas / (Gal[p].DiscGas[i]/area), FeedbackExponent);
                
            }
            else if(SupernovaRecipeOn == 2)
                //uniform reheated fraction
                reheated_mass = FeedbackReheatingEpsilon * stars;
            
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
            
            assert(stars+reheated_mass < 1.01*Gal[p].DiscGas[i]);
        }
        else // I haven't actually dealt with the situation of Supernovae being turned off here.  But do I even want to turn SN off?
        {
            reheated_mass = 0.0;
            ejected_mass = 0.0;
        }
        
        Gal[p].DiscSFR[i] += stars / dt;
        //^ stores sfr
        stars_sum += stars; //counting total stars formed accross annuli
        ejected_sum += ejected_mass; //gas ejected from halo around galaxy
        
        DiscPre = Gal[p].DiscGas[i];
        ColdPre = Gal[p].ColdGas; //checking values before the function
        
        // Update for star formation
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
        if(stars>=MIN_STARS_FOR_SN)
        {
            //builds temporary stellar disc
            NewStars[i] = (1 - RecycleFraction) * stars;
            
            //recyclefraction -> mass of stars just formed, immediately returned back to the gas
            //temporary disc also
            NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
        }
        else //no feedback in this case (rare)
        {
            NewStars[i] = stars;
            NewStarsMetals[i] = metallicity * stars;
        }
        if(!(NewStarsMetals[i] <= NewStars[i]))
        {
            printf("NewStars, metals = %e, %e\n", NewStars[i], NewStarsMetals[i]);
            printf("Gas, metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        }
        assert(NewStarsMetals[i] <= NewStars[i] + 0.00000000001);
        update_from_star_formation(p, stars, metallicity, i);
        //^ takes mass out of the gas disc to conserve mass
        
        if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
        assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
        
        if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
            reheated_mass = Gal[p].DiscGas[i];
        
        // These checks ensure numerical uncertainties don't blow up
        
        assert(fabs(Gal[p].ColdGas-ColdPre) <= 1.01*fabs(Gal[p].DiscGas[i]-DiscPre) && fabs(Gal[p].ColdGas-ColdPre) >= 0.999*fabs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
        
        DiscPre = Gal[p].DiscGas[i];
        ColdPre = Gal[p].ColdGas;
        
        // Update from SN feedback
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
        assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
        if(reheated_mass > 0.0)
            update_from_feedback(p, centralgal, reheated_mass, metallicity, i);
        //^ taking reheated gas out of disc and into the hot gas halo
        assert(fabs(Gal[p].ColdGas-ColdPre) <= 1.01*fabs(Gal[p].DiscGas[i]-DiscPre) && fabs(Gal[p].ColdGas-ColdPre) >= 0.999*fabs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
        
        // Inject new metals from SN II
        if(SupernovaRecipeOn > 0 && stars>=MIN_STARS_FOR_SN)
        {
            Gal[p].DiscGasMetals[i] += Yield * stars*(1.0 - get_metallicity(NewStars[i],NewStarsMetals[i])); // Could just use metallicity variable here, surely
            Gal[p].MetalsColdGas += Yield * stars*(1.0 - get_metallicity(NewStars[i],NewStarsMetals[i]));
        }
        if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i]) printf("DiscGas, Metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
        if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
        assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
    }

    update_from_ejection(p, centralgal, ejected_sum);
    //^ updating masses after ejection
    return;
}
