#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "core_allvars.h"
#include "core_proto.h"
#define ln(x) log(x)

//new paramaters
//type: -1 : disk instability
//       0 : mergers
//       1 : starform + feedback

/* function returning the max between two numbers */
double NFW_profile(double radius, double R_vir, double M_vir, double r_2) {
    assert(radius>0);
    assert(R_vir>0);
    assert(M_vir>0);
    assert(r_2>0);
    double result = (M_vir/radius) * pow(ln((R_vir + r_2)/r_2)-(R_vir/(R_vir + r_2)),-1)
                            * ln(1 + radius/r_2);
    return result;
}

double CalcR_2(int p) {
    // Determine the distribution of dark matter in the halo =====
    double M_DM_tot = Gal[p].Mvir - Gal[p].HotGas - Gal[p].ColdGas - Gal[p].StellarMass - Gal[p].ICS - Gal[p].BlackHoleMass; // One may want to include Ejected Gas in this too
    double baryon_fraction = (Gal[p].HotGas + Gal[p].ColdGas + Gal[p].StellarMass + Gal[p].ICS + Gal[p].BlackHoleMass) / Gal[p].Mvir;
    
    //if(M_DM_tot < 0.0) M_DM_tot = 0.0;
    
    double X = log10(Gal[p].StellarMass/Gal[p].Mvir);
    
    double z, a, b, c_DM, c, r_2, rho_const;
    //R_s on wiki -> r_2 here!
    //p_o central density of profile.
    z = ZZ[Gal[p].SnapNum];
    if(z>5.0) z=5.0;
    a = 0.520 + (0.905-0.520)*exp_f(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
    b = -0.101 + 0.026*z; // Dutton & Maccio 2014
    c_DM = pow(10.0, a+b*log10(Gal[p].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
    if(Gal[p].Type==0 && Gal[p].StellarMass>0 && Gal[p].Mvir>0)
        c = c_DM * (1.0 + 3e-5*exp_f(3.4*(X+4.5))); // Di Cintio et al 2014b
    else
        c = 1.0*c_DM; // Should only happen for satellite-satellite mergers, where X cannot be trusted
    r_2 = Gal[p].Rvir / c; // Di Cintio et al 2014b
    //rho_const = M_DM_tot / (log((Gal[p].Rvir+r_2)/r_2) - Gal[p].Rvir/(Gal[p].Rvir+r_2));
    //assert(rho_const>=0.0);
    // ===========================================================
    return r_2;
}
/*
double CalcRho_0(double r_2, int p) {
    double M_DM_tot = Gal[p].Mvir - Gal[p].HotGas - Gal[p].ColdGas - Gal[p].StellarMass - Gal[p].ICS - Gal[p].BlackHoleMass; // One may want to include Ejected Gas in this too
    
    double rho_const = M_DM_tot / (log((Gal[p].Rvir+r_2)/r_2) - Gal[p].Rvir/(Gal[p].Rvir+r_2));
    return rho_const;
}
 */

//recipe essentially deals with updating all the necessay fields.
void recipe(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double metals_stars_sum, double strdotfull, double ejected_mass, double ejected_sum, double reheated_mass, double metallicity, double stars_angmom, int i, double stars, int feedback_type, double gas_sf, double V_rot)
{
    double fac, Sigma_0gas, DiscPre, ColdPre;
    double r_inner, r_outer, r_av, area, j_bin;
    check_channel_stars(p);
    r_inner = Gal[p].DiscRadii[i];
    r_outer = Gal[p].DiscRadii[i+1];
    r_av = sqrt((sqr(r_inner)+sqr(r_outer))/2.0);
    area = M_PI * (r_outer*r_outer - r_inner*r_inner);
    
    //gas_sf will simply read Gal[p].DiscGas[i]
    //for mergers and passive star formation only
    //V_rot is simply Gal[centralgal].Vvir
    
    //supernoveRecipeOn>0 -> stellar feedback
    if(SupernovaRecipeOn > 0 && ((Gal[p].DiscGas[i] > 0.0 && stars>MIN_STARS_FOR_SN) || feedback_type == -1)) //so we have enough mass to make supernovas
    {
        if(SupernovaRecipeOn == 1)
            //uses equation in the paper for supernova feedback
        {
            //finds Sigma_0gas the reference surface density
            Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / sqr(CM_PER_MPC/1e6 / UnitLength_in_cm);
            
            reheated_mass = FeedbackReheatingEpsilon * stars * pow(Sigma_0gas / (Gal[p].DiscGas[i]/area), FeedbackExponent);
            
            if(feedback_type == -1 && reheated_mass < MIN_STARFORMATION) {
                reheated_mass = 0.0;
            }
        } else if(SupernovaRecipeOn == 2) {
            //uniform reheated fraction
            reheated_mass = FeedbackReheatingEpsilon * stars;
            
        } else if(SupernovaRecipeOn == 3) { //new case with energy arguments
            /*
                c2h -> cold to hot     h2e -> hot to ejected
             
                e_in = e_supernovae = 1/2 * 630km/s ^2 * mass_stars_formed
                e_out = e_c2h + e_h2e
                related by
                    e_in = free_epsilon1 * e_out &&
                    EonM_cold = 0.5 * (v'/r)^2
                    EonM_hot = 0.5 * V_vir^2
                    E_h2e = M * gravity_well
             
             %-------------- Units, son! ---------------
             %------------------------------------------
             UnitLength_in_cm            3.08568e+24        ;WATCH OUT: Mpc/h
             UnitMass_in_g                1.989e+43        ;WATCH OUT: 10^10Msun/h
             UnitVelocity_in_cm_per_s    100000            ;WATCH OUT: km/s
            */
            
            double energy_efficiency = 0.3;   //energy_out/energy_in
            //EnergyDispersion in ParamFile defined as energy_h2e/energy_c2h
            
            double energy_in =  0.5 * 630 * 630 * stars;
                                //e_SN, Croton et al. pg8
            double energy_out = energy_efficiency * energy_in;
            
            //Energy per unit mass for each state: Cold, Hot, Ejected
            assert(V_rot>0);
            assert(Gal[p].Rvir>0);
            double energy_perunit_coldmass =
                0.5 * (V_rot*r_outer/Gal[p].Rvir) * (V_rot*r_outer/Gal[p].Rvir) + NFW_profile((r_outer+r_inner)/2, Gal[p].Rvir, Gal[p].Mvir, CalcR_2(p));
            
            double r_hot = 0.5 * Gal[p].HotGas * Gal[p].Rvir;
            double energy_perunit_hotmass =
                0.5 * V_rot * V_rot + NFW_profile((r_hot)/2, Gal[p].Rvir, Gal[p].Mvir, CalcR_2(p));
            
            double energy_perunit_ejectedmass =
                Gal[p].HotGas * 0.5 * V_rot * V_rot + NFW_profile(Gal[p].Rvir, Gal[p].Rvir, Gal[p].Mvir, CalcR_2(p));
            
            assert(energy_perunit_hotmass>energy_perunit_coldmass);
            double energy_c2h = energy_perunit_hotmass - energy_perunit_coldmass;
            double energy_h2e = energy_perunit_ejectedmass - energy_perunit_hotmass;
            
            reheated_mass = energy_out/energy_c2h;
            //assert(reheated_mass < Gal[p].DiscGas[i]);
            printf("DiscGas: %f, Stars: %f, ReheatedMass: %f, EnergyC2h: %f\n", gas_sf, stars, reheated_mass, energy_c2h);
            
        } else {
            reheated_mass = 0.0;
            assert(reheated_mass == reheated_mass && reheated_mass != INFINITY);
        }
        
        // Can't use more cold gas than is available, so balance SF and feedback
        if((stars + reheated_mass) > gas_sf && (stars + reheated_mass) > 0.0)
        {
            fac = gas_sf / (stars + reheated_mass);
            stars *= fac;
            reheated_mass *= fac;
            assert(reheated_mass == reheated_mass && reheated_mass != INFINITY);
        }
        
        if(stars<MIN_STARS_FOR_SN)
        {
            stars = MIN_STARS_FOR_SN;
            reheated_mass = gas_sf - stars;

            if(feedback_type == -1 && gas_sf < MIN_STARS_FOR_SN) {
                stars = gas_sf;
                reheated_mass = 0.0;
            }
            ejected_mass = 0.0;
            
            if(feedback_type != -1 && reheated_mass < MIN_STARFORMATION) {
                reheated_mass = 0.0; // Limit doesn't have to be the same as MIN_STARFORMATION, but needs to be something reasonable
            }
        }
        
        //equation already formed from energy arguents
        else {
             ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (V_rot * V_rot) - FeedbackReheatingEpsilon) * stars;
        }
        
        if(ejected_mass < MIN_STARFORMATION) {
            ejected_mass = 0.0;
            assert(stars+reheated_mass < 1.01*Gal[p].DiscGas[i]);
        }
    }
    else // I haven't actually dealt with the situation of Supernovae being turned off here.  But do I even want to turn SN off?
    {
        reheated_mass = 0.0;
        ejected_mass = 0.0;
    }
    
    if (feedback_type != -1) { //stores star formation rate
        Gal[p].DiscSFR[i] += stars / dt;
    }
    
    DiscPre = Gal[p].DiscGas[i];
    ColdPre = Gal[p].ColdGas; //checking values before the function

    if (feedback_type != -1) {// Update for star formation
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
    
    //take mass out of the gas disc to conserve mass
    update_from_star_formation(p, stars, metallicity, i);
    
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
    if(reheated_mass > 0.0 || feedback_type == -1)
        update_from_feedback(p, centralgal, reheated_mass, metallicity, i);
    //^ taking reheated gas out of disc and into the hot gas halo
    
    // Inject new metals from SN II
    if(SupernovaRecipeOn > 0 && stars >= MIN_STARS_FOR_SN)
    {
        Gal[p].DiscGasMetals[i] += Yield * stars * (1.0 - metallicity);
        Gal[p].MetalsColdGas += Yield * stars * (1.0 - metallicity);
    }
    
    if (feedback_type == -1) {//updating masses after ejection
    update_from_ejection(p, centralgal, ejected_sum);
    }
    
    //angular momentem adjustments
    if (feedback_type == 0)  {
        j_bin = (DiscBinEdge[i]+DiscBinEdge[i+1])*0.5;
        if(stars>=MIN_STARS_FOR_SN) {
            stars_angmom += (1 - RecycleFraction) * stars * j_bin;
        } else {
            stars_angmom += stars * j_bin;
        }
    }
}

    

void feedback(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double metals_stars_sum, double strdotfull, double stars_angmom, int mode, double eburst, double disc_mass_ratio[N_BINS], int feedback_type)
{
    double strdot, stars, reheated_mass, ejected_mass, metallicity, area, SFE_H2, f_H2_const, DiscPre, ColdPre;
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
        if (feedback_type == 1)
        {
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
                    if(Gal[p].ColdGas>0.0) {
                        strdot = strdotfull * Gal[p].DiscGas[i] / Gal[p].ColdGas; // * Gal[p].DiscH2[i] / H2sum;
                    } else {
                        strdot = 0.0;
                    }
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
        recipe(p, centralgal, dt, step, NewStars, NewStarsMetals, stars_sum, metals_stars_sum, strdotfull, ejected_mass, ejected_sum, reheated_mass, metallicity, stars_angmom, i, stars, feedback_type, Gal[p].DiscGas[i], Gal[centralgal].Vvir);
    }
    if(ejected_sum>0.0) {//updating masses after ejection
        update_from_ejection(p, centralgal, ejected_sum);
    }
    return;
}

