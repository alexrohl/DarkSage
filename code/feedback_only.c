#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "core_allvars.h"
#include "core_proto.h"

//new paramaters
//type: -1 : disk instability
//       0 : mergers
//       1 : starform + feedback


//recipe essentially deals with updating all the necessay fields.
void recipe(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double metals_stars_sum, double strdotfull, double ejected_mass, double ejected_sum, double reheated_mass, double metallicity, double stars_angmom, int i, double stars, int feedback_type, double gas_sf, double V_rot)
{
    double strdot, fac, SFE_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscStarsSum, DiscPre, ColdPre;
    double r_inner, r_outer, area, Vvir, j_bin;
    double reff, tdyn, cold_crit, H2sum; // For SFprescription==3
    double StarsPre = Gal[p].StellarMass;
    check_channel_stars(p);
    
    //gas_sf will simply read Gal[p].DiscGas[i]
    //for mergers and passive star formation only
    
    //supernoveRecipeOn>0 -> stellar feedback
    if(SupernovaRecipeOn > 0 && Gal[p].DiscGas[i] > 0.0 && stars>MIN_STARS_FOR_SN) //too low mass to make supernovas
    {
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        area = M_PI * (r_outer*r_outer - r_inner*r_inner);
        
        if(SupernovaRecipeOn == 1)
            //!!!uses equation in the paper for supernova feedback
        {
            //finds Sigma_0gas the reference surface density
            Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / sqr(CM_PER_MPC/1e6 / UnitLength_in_cm);
            
            reheated_mass = FeedbackReheatingEpsilon * stars * pow(Sigma_0gas / (Gal[p].DiscGas[i]/area), FeedbackExponent);
            
            if(feedback_type == -1 && reheated_mass < MIN_STARFORMATION) {
                
                reheated_mass = 0.0;
            }
        }
        else if(SupernovaRecipeOn == 2)
            //uniform reheated fraction
            reheated_mass = FeedbackReheatingEpsilon * stars;
        else
            reheated_mass = 0.0;
            assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
        // Can't use more cold gas than is available, so balance SF and feedback
        if((stars + reheated_mass) > gas_sf && (stars + reheated_mass) > 0.0)
        {
            fac = gas_sf / (stars + reheated_mass);
            stars *= fac;
            reheated_mass *= fac;
            assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
        }
        
        if(stars<MIN_STARS_FOR_SN)
        {
            stars = MIN_STARS_FOR_SN;
            reheated_mass = gas_sf - stars; // Used to have (1-RecycleFraction)* in front of stars here, but changed philosophy

            if(feedback_type == -1 && gas_sf < MIN_STARS_FOR_SN) {
                stars = gas_sf;
                reheated_mass = 0.0;
            }
            ejected_mass = 0.0;
        
        if(feedback_type != -1 && reheated_mass < MIN_STARFORMATION)
            reheated_mass = 0.0; // Limit doesn't have to be the same as MIN_STARFORMATION, but needs to be something reasonable
        }
        
        //equation already formed from energy arguents
        else {
            if (feedback_type != -1) {
                ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) - FeedbackReheatingEpsilon) * stars;
            } else {
             ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (V_rot * V_rot) - FeedbackReheatingEpsilon) * stars;
        }
        if(ejected_mass < MIN_STARFORMATION)
            ejected_mass = 0.0;
        
        assert(stars+reheated_mass < 1.01*Gal[p].DiscGas[i]);
        }
    }
    else // I haven't actually dealt with the situation of Supernovae being turned off here.  But do I even want to turn SN off?
    {
        reheated_mass = 0.0;
        ejected_mass = 0.0;
    }
    
    if (feedback_type != -1) {
        Gal[p].DiscSFR[i] += stars / dt;
            //^ stores sfr
    }
    
    DiscPre = Gal[p].DiscGas[i];
    ColdPre = Gal[p].ColdGas; //checking values before the function

    if (feedback_type != -1) {
        
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
    
        stars_sum += stars; //counting total stars formed accross annuli
    }
    
    ejected_sum += ejected_mass; //gas ejected from halo around galaxy
    
    if (feedback_type == 0) {
        metals_stars_sum += NewStarsMetals[i];
    }
    
    update_from_star_formation(p, stars, metallicity, i);
    //^ takes mass out of the gas disc to conserve mass
    
    if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
        reheated_mass = Gal[p].DiscGas[i];
    
    // These checks ensure numerical uncertainties don't blow up
    assert(fabs(Gal[p].ColdGas-ColdPre) <= 1.01*fabs(Gal[p].DiscGas[i]-DiscPre) && fabs(Gal[p].ColdGas-ColdPre) >= 0.999*fabs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
    
    DiscPre = Gal[p].DiscGas[i];
    ColdPre = Gal[p].ColdGas;
        
    //need the new matallicity in the disk instabilty case
    double metallicity_old = 1.0*metallicity;
    if (feedback_type == -1) {
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    }
        
    // Update from SN feedback
    if(reheated_mass > 0.0)
        update_from_feedback(p, centralgal, reheated_mass, metallicity, i);
    //^ taking reheated gas out of disc and into the hot gas halo
    
    // Inject new metals from SN II
    if(SupernovaRecipeOn > 0 && stars >= MIN_STARS_FOR_SN)
    {
        Gal[p].DiscGasMetals[i] += Yield * stars * (1.0 - metallicity_old); // Could just use metallicity variable here, surely
        Gal[p].MetalsColdGas += Yield * stars*(1.0 - metallicity_old);
    }
    
    if (feedback_type == -1){
    update_from_ejection(p, centralgal, ejected_sum);
    }
    
    //angular momentem adjustments
    if (feedback_type == 0)  {
        j_bin = (DiscBinEdge[i]+DiscBinEdge[i+1])*0.5;
        if(stars>=MIN_STARS_FOR_SN) {
            stars_angmom += (1 - RecycleFraction) * stars * j_bin;
            }
        else {
            stars_angmom += stars * j_bin;
        }
    }
}

    

void feedback(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double metals_stars_sum, double strdotfull, double stars_angmom, int mode, double eburst, double disc_mass_ratio[N_BINS], int feedback_type)
{
    double strdot, stars, reheated_mass, ejected_mass, fac, metallicity, area, SFE_H2, f_H2_const, Sigma_0gas, DiscGasSum, DiscStarsSum, DiscPre, ColdPre;
    double r_inner, r_outer;
    double reff, tdyn, cold_crit, H2sum; // For SFprescription==3
    
    int i;
    double StarsPre = Gal[p].StellarMass;
    check_channel_stars(p);
    
    double ejected_sum = 0.0;
    double gas_sf = 0.0;
    reheated_mass = 0.0; // initialise
    strdot = 0.0;
    
    // Loops through annuli
    for(i=0; i<N_BINS; i++)
        
    {
        //printf("reached here %d\n", feedback_type);
        if (feedback_type == 1) {
            //printf("reached here2\n");
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
        }
        else if (feedback_type == 0) {
            if(mode == 1)
                eburst = disc_mass_ratio[i];
            else
                //equation in paper, model for stars forming in burst
                eburst = BetaBurst * pow(disc_mass_ratio[i], AlphaBurst);
            
            stars = eburst * Gal[p].DiscGas[i];
            if(stars < MIN_STARFORMATION)
                stars = 0.0;
            //safety net
            if(stars > Gal[p].DiscGas[i])
                stars = Gal[p].DiscGas[i];
        }
        
        //Calls recipe for each disk,
        recipe(p, centralgal, dt, step, NewStars, NewStarsMetals, stars_sum, metals_stars_sum, strdotfull, ejected_mass, ejected_sum, reheated_mass, metallicity, stars_angmom, i, stars, feedback_type, Gal[p].DiscGas[i], 0);
    }
    if(ejected_sum>0.0)
        update_from_ejection(p, centralgal, ejected_sum);
    //^ updating masses after ejection
        
    return;
}

