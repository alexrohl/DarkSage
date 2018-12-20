#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int step)
{
    double mi, ma, mass_ratio;
    double disc_mass_ratio[N_BINS], PostRetroGas[N_BINS];
    int i;
    
    for(i=0; i<N_BINS; i++)
        disc_mass_ratio[i] = 0.0;
    
    // Calculate mass ratio of merging galaxies
    if(Gal[p].StellarMass + Gal[p].ColdGas <
       Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas)
    {
        mi = Gal[p].StellarMass + Gal[p].ColdGas;
        ma = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
    }
    else
    {
        mi = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
        ma = Gal[p].StellarMass + Gal[p].ColdGas;
    }
    
    if(ma > 0)
        mass_ratio = mi / ma;
    else
        mass_ratio = 1.0;
    
    if(mass_ratio<0.0)
    {
        mass_ratio = 0.0;
        printf("Had to correct mass_ratio < 0.0");
    }
    
    add_galaxies_together(merger_centralgal, p, mass_ratio, disc_mass_ratio, PostRetroGas);
    
    collisional_starburst_recipe(disc_mass_ratio, merger_centralgal, centralgal, dt, 0, step);
    
    double BHaccrete = grow_black_hole(merger_centralgal, disc_mass_ratio);
    
    if(AGNrecipeOn>0)
        quasar_mode_wind(p, BHaccrete, centralgal);
    
    // Check whether any retrograde gas is left over
    double unstable_gas, metallicity, stars, net_stars;
    for(i=N_BINS-1; i>=0; i--)
    {
        metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
        
        if(PostRetroGas[i] < 0.99*Gal[merger_centralgal].DiscGas[i] && Gal[merger_centralgal].DiscGas[i]-PostRetroGas[i] > 1e-10)
        {
            unstable_gas = Gal[merger_centralgal].DiscGas[i] - PostRetroGas[i];
            stars = deal_with_unstable_gas(unstable_gas, merger_centralgal, i, Gal[merger_centralgal].Vvir, metallicity, centralgal, Gal[merger_centralgal].DiscRadii[i], Gal[merger_centralgal].DiscRadii[i+1]);
            
            if(stars>=MIN_STARS_FOR_SN)
                net_stars = (1 - RecycleFraction) * stars;
            else
                net_stars = stars;
            
            Gal[merger_centralgal].StellarMass += net_stars;
            Gal[merger_centralgal].MetalsStellarMass += metallicity * net_stars;
            Gal[merger_centralgal].StarsMergeBurst += net_stars;
            Gal[merger_centralgal].SfrMerge[step] += stars / dt;
            check_channel_stars(merger_centralgal);
            
            // Add the new stars from the retrograde starburst to the classical bulge
            // No longer carry any net AM into it
            if(net_stars>0)
            {
                Gal[merger_centralgal].ClassicalBulgeMass += net_stars;
                Gal[merger_centralgal].ClassicalMetalsBulgeMass += metallicity * net_stars;
            }
            
        }
    }
    
    if(mass_ratio > ThreshMajorMerger)
    {
        stars_to_bulge(merger_centralgal);
        Gal[merger_centralgal].LastMajorMerger = time;
        Gal[p].mergeType = 2;  // Mark as major merger
    }
    else
    {
        Gal[merger_centralgal].LastMinorMerger = time;
        Gal[p].mergeType = 1;  // Mark as minor merger
    }
    
    
    if(DiskInstabilityOn>0)
        check_disk_instability(merger_centralgal, centralgal, dt, step);
    else
        update_stellardisc_scaleradius(p); // will already be done within check_disk_instability otherwise
}



double grow_black_hole(int merger_centralgal, double* disc_mass_ratio)
{
    double BHaccrete, BHaccrete_tot, metallicity;//, accrete_ratio, DiscGasSum;
    int i;
    
    BHaccrete_tot = 0.0;
    
    for(i=0; i<N_BINS; i++)
    {
        if(Gal[merger_centralgal].DiscGas[i] > 0.0)
        {
            BHaccrete = BlackHoleGrowthRate * disc_mass_ratio[i] / (1.0 + sqr(280.0 / Gal[merger_centralgal].Vvir)) * Gal[merger_centralgal].DiscGas[i];
            
            if(BHaccrete > Gal[merger_centralgal].DiscGas[i]) // This could only be possible if BlackHoleGrowthRate is set to >1.0, which shouldn't happen...
            {
                BHaccrete_tot += Gal[merger_centralgal].DiscGas[i];
                Gal[merger_centralgal].ColdGas -= Gal[merger_centralgal].DiscGas[i];
                Gal[merger_centralgal].MetalsColdGas -= Gal[merger_centralgal].DiscGasMetals[i];
                Gal[merger_centralgal].DiscGas[i] = 0.0;
                Gal[merger_centralgal].DiscGasMetals[i] = 0.0;
            }
            else
            {
                BHaccrete_tot += BHaccrete;
                metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
                Gal[merger_centralgal].DiscGas[i] -= BHaccrete;
                Gal[merger_centralgal].DiscGasMetals[i] -= BHaccrete * metallicity;
                if(Gal[merger_centralgal].DiscGasMetals[i]<0.0) Gal[merger_centralgal].DiscGasMetals[i] = 0.0;
                Gal[merger_centralgal].ColdGas -= BHaccrete;
                Gal[merger_centralgal].MetalsColdGas -= BHaccrete * metallicity;
                
            }
        }
    }
    
    Gal[merger_centralgal].BlackHoleMass += BHaccrete_tot;
    return BHaccrete_tot;
}
