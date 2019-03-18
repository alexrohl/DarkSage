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

//Function to deal with energy distribution
//now implemented with inner/outer argument
//height = 11^2 / 2*PI*G*(gas + stars)
void disperse_energy_outward(double energy, double annuli_energy[N_BINS], int affected_index) {
    if (affected_index < N_BINS) {
        double new_energy = EnergyEfficiencyC2H * energy;
        double excess_energy = energy - new_energy;
        annuli_energy[affected_index] += new_energy;
        //recursive call
        disperse_energy_outward(excess_energy, annuli_energy, affected_index+1);
    }
    return;
}

void disperse_energy_inward(double energy, double annuli_energy[N_BINS], int affected_index) {
    if (affected_index >= 0) {
        double new_energy = EnergyEfficiencyC2H * energy;
        double excess_energy = energy - new_energy;
        annuli_energy[affected_index] += new_energy;
        //recursive call
        disperse_energy_inward(excess_energy, annuli_energy, affected_index-1);
    }
    return;
}

void distribute_energy(double all_stars[N_BINS], double annuli_energy_init[N_BINS], int p) {
    int i;
    double annuli_energy_final[N_BINS];
    for(i=0; i<N_BINS; i++) {
        annuli_energy_final[i] = 0.0;
        //printf("EnergyIn: %f\n", annuli_energy_init[i]);
    }
    
    //to check if we are working in inner-most or outer-most radii
    for(i=0; i<N_BINS; i++)
    {
        int inner;
        int outer;
        double h_inner;
        double h_outer;
        if (i==0) {//inner radius will receive the inner amount for now
            inner = i;
            outer = i+1;
        } else if (i < (N_BINS-1)) {
            inner = i-1;
            outer = i+1;
        } else {//outer radius will receive the outer distribution as well for now
            inner = i-1;
            outer = i;
        }
        if (all_stars[inner]+Gal[p].DiscGas[inner] <= 0.0) {
            h_inner = 0.0;
        } else {
            h_inner = 11*11 / (2*M_PI * GRAVITY * (all_stars[inner]+Gal[p].DiscGas[inner]));
        }
        if (all_stars[outer]+Gal[p].DiscGas[outer] <= 0.0) {
            h_outer = 0.0;
        } else {
            h_outer = 11*11 / (2*M_PI * GRAVITY * (all_stars[outer]+Gal[p].DiscGas[outer]));
        }
        
        double r_inner = 1.0*Gal[p].DiscRadii[i];
        double r_outer = 1.0*Gal[p].DiscRadii[i+1];
        double r_av = sqrt((sqr(r_inner)+sqr(r_outer))/2.0);
        double energy;
        if (annuli_energy_init[i] > 0.0) {
            energy = 1.0*annuli_energy_init[i];
        } else {
            energy = 0.0;
        }
        double innerVol = (M_PI/3) * (r_av - r_inner) * (r_av + 2*r_inner) * h_inner;
        double outerVol = (M_PI/3) * (r_outer - r_av) * (r_av + 2*r_outer) * h_outer;
        double totalVol = (M_PI/3) * (r_outer - r_inner) * (h_inner*(r_inner + 2*r_outer) + h_outer*(r_outer + 2*r_inner));
        double remainVol = 0.0;
        double inner_prop = 0.0;
        double remaining_prop = 0.0;
        double outer_prop = 0.0;
        if (totalVol > 0.0) {
            remainVol = totalVol - innerVol - outerVol;
            inner_prop = innerVol / totalVol;
            remaining_prop = remainVol / totalVol;
            outer_prop = outerVol / totalVol;
        } else {
            //printf("total:%f", totalVol);

        }
        assert(r_av >= r_inner);
        //printf("total:%f", totalVol);
        double area = M_PI * (r_outer*r_outer - r_inner*r_inner);
        
        assert(outer<N_BINS);
        
        //Disperse inward energy towards the center
        disperse_energy_inward(energy * inner_prop, annuli_energy_final, inner);
        //Update the annulus energy based on the remaining proportion
        annuli_energy_final[i]     += energy * remaining_prop;
        //Disperse outward energy towards edge of galaxy
        disperse_energy_outward(energy * outer_prop, annuli_energy_final, outer);
    }
    for(i=0; i<N_BINS; i++)
    {
        //printf("energy_init: %f, energy_final: %f\n", annuli_energy_init[i], annuli_energy_final[i]);
        if (!(annuli_energy_final[i] < 0.0 || isnan(annuli_energy_final[i]))) {
            annuli_energy_init[i] = 1.0*annuli_energy_final[i];
        } else {
            annuli_energy_init[i] = 0.00000;
        }
        assert(annuli_energy_init[i] >= 0);
        //printf("EnergyOut: %f\n", annuli_energy_init[i]);
    }
    //recursive call
    
}

//recipe essentially deals with updating all the necessay fields.
struct RecipeOutput recipe_dispersed(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double metals_stars_sum, double strdotfull, double ejected_mass, double ejected_sum, double reheated_mass, double metallicity, double stars_angmom, int i, double stars, int feedback_type, double gas_sf, double V_rot, double annuli_energy)
{
    double fac, Sigma_0gas, DiscPre, ColdPre;
    double r_inner, r_outer, r_av, area, j_bin, energy_h2e, energy_excess, energy_efficiency2, energy_out, energy_c2h;
    //check_channel_stars(p);
    r_inner = Gal[p].DiscRadii[i];
    r_outer = Gal[p].DiscRadii[i+1];
    r_av = sqrt((sqr(r_inner)+sqr(r_outer))/2.0);
    area = M_PI * (r_outer*r_outer - r_inner*r_inner);
    //printf("feedback type: %d, %f, %f, %f, %f\n", feedback_type, stars, stars_sum, ejected_mass, ejected_sum);
    //gas_sf will simply read Gal[p].DiscGas[i]
    //for mergers and passive star formation only
    //V_rot is simply Gal[centralgal].Vvir
    //supernoveRecipeOn>0 -> stellar feedback
    if(SupernovaRecipeOn > 0 && ((Gal[p].DiscGas[i] > 0.0 && stars>MIN_STARS_FOR_SN) || feedback_type == -1)) //so we have enough mass to make supernovas
    {
        double check = Gal[p].DiscGas[i] * V_rot * Gal[p].HotGas * r_outer * Gal[p].Rvir * Gal[p].Mvir;

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
            
        } else if(SupernovaRecipeOn == 3 && check != 0) { //new case
             /* %-------------- Units, son! ---------------
             %------------------------------------------
             UnitLength_in_cm            3.08568e+24        ;WATCH OUT: Mpc/h
             UnitMass_in_g                1.989e+43        ;WATCH OUT: 10^10Msun/h
             UnitVelocity_in_cm_per_s    100000            ;WATCH OUT: km/s
            */
            
            //EnergyEfficiencyC2H: energy_out/energy_in
            //EnergyEfficiencyH2E: proportion of excess energy going into ejecting gas

            double energy_out = EnergyEfficiencyC2H * annuli_energy;
            
            //Energy per unit mass for each state: Cold, Hot, Ejected
            double energy_perunit_coldmass =
                0.5 * (V_rot*r_av/Gal[p].Rvir) * (V_rot*r_av/Gal[p].Rvir)
                + NFW_profile(r_av, Gal[p].Rvir, Gal[p].Mvir, CalcR_2(p));
            
            double r_hot = 0.5 * Gal[p].Rvir;
            double energy_perunit_hotmass =
                0.5 * V_rot * V_rot
                + NFW_profile(r_hot, Gal[p].Rvir, Gal[p].Mvir, CalcR_2(p));
            
            double energy_perunit_ejectedmass =
                0.5 * V_rot * V_rot
                + NFW_profile(Gal[p].Rvir, Gal[p].Rvir, Gal[p].Mvir, CalcR_2(p));
            
            //assert(energy_perunit_hotmass>energy_perunit_coldmass);

            //Find the energy required to move between states
            double energy_c2h = energy_perunit_hotmass - energy_perunit_coldmass;
            double energy_h2e = energy_perunit_ejectedmass - energy_perunit_hotmass;
            assert(energy_h2e>=0.0);
            assert(energy_out>=0.0);

            if (energy_c2h <= 0.0 && energy_out > 0.0)   {
                reheated_mass = 1.0*gas_sf;
            } else {
                reheated_mass = energy_out/energy_c2h;
                //printf("high potential energy cold gas here, radius ratio: %f, reheat/stars ratio: %f\n",r_av/Gal[p].Rvir,reheated_mass/stars);
            }
            //assert(reheated_mass < Gal[p].DiscGas[i]);
            //printf("DiscGas: %f, Stars: %f, ReheatedMass: %f,EnergyOut: %f, EnergyC2h: %f\n", gas_sf, stars, reheated_mass,energy_out, energy_c2h);
            
        } else {
            reheated_mass = 0.0;
            assert(reheated_mass == reheated_mass && reheated_mass != INFINITY);
        }
        // Can't use more cold gas than is available, so balance SF and feedback
        if((stars + reheated_mass) > gas_sf && (stars + reheated_mass) > 0.0)
        {
            //printf("feedback: %d, reheatedmass: %f\n", feedback_type, reheated_mass);
            fac = gas_sf / (stars + reheated_mass);
            stars *= fac;
            reheated_mass *= fac;
            assert(reheated_mass == reheated_mass && reheated_mass != INFINITY);
            energy_excess = energy_out - reheated_mass * energy_c2h;
            if (energy_excess < 0.0) energy_excess = 0.0;
            assert(energy_excess>=0.0);
            assert(fac<1.0);
            assert(reheated_mass>0.0);
        }
        
        if(stars<MIN_STARS_FOR_SN)
        {
            stars = 1.0*MIN_STARS_FOR_SN;
            reheated_mass = gas_sf - stars;

            if(feedback_type == -1 && gas_sf < MIN_STARS_FOR_SN) {
                stars = 1.0*gas_sf;
                reheated_mass = 0.0;
            }
            ejected_mass = 0.0;
            
            if(feedback_type != -1 && reheated_mass < MIN_STARFORMATION) {
                reheated_mass = 0.0; // Limit doesn't have to be the same as MIN_STARFORMATION, but needs to be something reasonable
            }
        }
        
        
        else {
            //equation already formed from energy arguents
            if (SupernovaRecipeOn == 3) {
                //new ejected_mass argument
                //ejected_mass = energy_efficiency2 * (gas_sf/Gal[p].HotGas) * (energy_excess/energy_h2e);
                if (energy_excess > 0 && energy_h2e > 0)  {
                    ejected_mass = (energy_efficiency2 * energy_excess) / energy_h2e;
                } else {
                    ejected_mass = 0.0;
                }
            } else {
                ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (V_rot * V_rot) - FeedbackReheatingEpsilon) * stars;
                //printf("\ntrolling\n");
            }
            
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
    
    //checks
    if (energy_out > 0.001) {
        printf("TYPE: %d, EjectedMass: %f, ReheatedMass: %f, EnergyH2E: %f, EnergyC2H: %f, EnergyOut: %f, EnergyExcess: %f, FAC: %f\n", feedback_type, ejected_mass, reheated_mass, energy_h2e, energy_c2h, energy_out, energy_excess, fac);
    } else {
        //printf("STARS2: %f\n", stars);
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
    if (feedback_type != 0) {
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    }
        
    // Update from SN feedback
    if(reheated_mass > 0.0 || feedback_type == -1)
        update_from_feedback(p, centralgal, reheated_mass, metallicity, i);
    //^ taking reheated gas out of disc and into the hot gas halo
    
    // Inject new metals from SN II
    if (feedback_type == 1) {
        metallicity = get_metallicity(NewStars[i],NewStarsMetals[i]);
    }
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
    struct RecipeOutput output = {stars, ejected_mass, NewStars[i], NewStarsMetals[i]};
    //printf("ejected: %f, stars: %f\n", ejected_mass, stars);
    return output;
}

    

void feedback_dispersed(int p, int centralgal, double dt, int step, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double stars_sum, double metals_stars_sum, double strdotfull, double stars_angmom, int mode, double eburst, double disc_mass_ratio[N_BINS], int feedback_type)
{
    //struct RecipeOutput output = recipe(p, centralgal, dt, step, NewStars, NewStarsMetals, stars_sum, metals_stars_sum, strdotfull, ejected_mass, ejected_sum, reheated_mass, metallicity, stars_angmom, i, stars, feedback_type, Gal[p].DiscGas[i], Gal[centralgal].Vvir);
    double strdot, stars, reheated_mass, ejected_mass, metallicity, area, SFE_H2, f_H2_const, DiscPre, ColdPre, DiscGasSum;
    int i;
    double StarsPre = Gal[p].StellarMass;
    check_channel_stars(p);
    f_H2_const = 1.38e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*H2FractionExponent);
    SFE_H2 = 7.75e-4 * UnitTime_in_s / SEC_PER_MEGAYEAR; // This says if SfrEfficiency==1.0, then the true efficiency of SF from H2 is 7.75e-4 h Myr^-1
    double ejected_sum = 0.0;
    double gas_sf = 0.0;
    reheated_mass = 0.0; // initialise
    strdot = 0.0;
    ejected_sum = 0.0;
    stars_sum = 0.0;
    double all_stars[N_BINS];
    double annuli_energy_init[N_BINS];
    double annuli_energy_final[N_BINS];
    // Loops through annuli
    for(i=0; i<N_BINS; i++)
    {
        double r_inner = Gal[p].DiscRadii[i];
        double r_outer = Gal[p].DiscRadii[i+1];
        area = M_PI * (r_outer*r_outer - r_inner*r_inner);
        
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
        //we will store stars in the all_stars array and the energy dispersion in annuli_energy array
        all_stars[i] = 1.0*stars;
        annuli_energy_init[i] = 0.5 * 630 * 630 * stars; //e_SN, Croton et al. pg8
    }
    
    //redistribute energy here:
    int j;
    for(i=0; i<30; i++)
    {
        distribute_energy(all_stars, annuli_energy_init, p);
        /*printf("iteration: %d\n", i);
        for(j=0; j<N_BINS; j++)
        {
            printf("energy: %d, %f\n",j,annuli_energy_init[j]);
        }
         */
    }
    //energy now redistrubted.
    
    for(i=0; i<N_BINS; i++)
    {
        double energy = 1.0*annuli_energy_init[i];
        double stars = 1.0*all_stars[i];
        struct RecipeOutput output = recipe_dispersed(p, centralgal, dt, step, NewStars, NewStarsMetals, stars_sum, metals_stars_sum, 0.0, ejected_mass, ejected_sum, reheated_mass, metallicity, stars_angmom, i, stars, feedback_type, Gal[p].DiscGas[i], Gal[centralgal].Vvir, energy);
        
        //printf("ejected2: %f, stars2: %f\n", output.ejected, output.stars);
        stars = 1.0*output.stars;
        stars_sum += stars;
        ejected_sum += 1.0*output.ejected;
        NewStars[i] = 1.0*output.NewStarsDisk;
        NewStarsMetals[i] = 1.0*output.NewMetalsDisk;
    }
    //printf("Ejected Sum: %f, Stars Sum: %f\n", ejected_sum, stars_sum);
    if(ejected_sum>0.0) {//updating masses after ejection
        update_from_ejection(p, centralgal, ejected_sum);
    }
    if (feedback_type == 1) {
        double NewStarSum = 0.0;
        for(i=N_BINS-1; i>=0; i--) NewStarSum += NewStars[i];
        //printf("NewStarSum: %f\n",NewStarSum);
        // Sum stellar discs together
        if(NewStarSum>0.0)
            combine_stellar_discs(p, NewStars, NewStarsMetals);
    
        // Update the star formation rate
        Gal[p].SfrFromH2[step] += stars_sum / dt;
        Gal[p].StarsFromH2 += NewStarSum;
    
        if(Gal[p].StellarMass >= MIN_STARS_FOR_SN)
        {
            check_channel_stars(p);
            assert(Gal[p].StellarMass >= (StarsPre + NewStarSum)/1.01 && Gal[p].StellarMass <= (StarsPre + NewStarSum)*1.01);
        }
    
        DiscGasSum = get_disc_gas(p);
        assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
        assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
    }
    else if (feedback_type == 0) {
        int merger_centralgal = p;
        if(stars_sum>0)
        {
            // Update bulge spin
            int s;
            for(s=0; s<3; s++)
            {
                Gal[merger_centralgal].SpinClassicalBulge[s] = (Gal[merger_centralgal].SpinClassicalBulge[s]*Gal[merger_centralgal].ClassicalBulgeMass + Gal[merger_centralgal].SpinGas[s]*stars_angmom) / (Gal[merger_centralgal].ClassicalBulgeMass+stars_sum);
                assert(Gal[merger_centralgal].SpinClassicalBulge[s] == Gal[merger_centralgal].SpinClassicalBulge[s] && Gal[merger_centralgal].SpinClassicalBulge[s] != INFINITY && Gal[merger_centralgal].SpinClassicalBulge[s] != -INFINITY);
            }
            
            // Now adding all new stars directly to the bulge
            Gal[merger_centralgal].StellarMass += stars_sum; // Recycling fraction already taken into account when adding to stars_sum etc above
            Gal[merger_centralgal].ClassicalBulgeMass += stars_sum;
            Gal[merger_centralgal].MetalsStellarMass += metals_stars_sum;
            Gal[merger_centralgal].ClassicalMetalsBulgeMass += metals_stars_sum;
        }
        
        Gal[merger_centralgal].SfrMerge[step] += stars_sum / dt;
        Gal[merger_centralgal].StarsMergeBurst += stars_sum;
        
        check_channel_stars(merger_centralgal);
    }
    for (int i = 0; i < N_BINS; i++)
        annuli_energy_final[i] = 0.0;
    return;
}
