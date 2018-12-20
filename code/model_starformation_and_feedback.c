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
  //printf("original f here");
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
    
  update_HI_H2(p);
  //recomputes fraction of gas in form of h1 h2 hydrogren
    
  if(SFprescription==2) // Old School Prescription based on SAGE
  {
      reff = 3.0 * Gal[p].DiskScaleRadius;
      tdyn = reff / Gal[p].Vvir;
      cold_crit = 0.19 * Gal[p].Vvir * reff;
      if(Gal[p].ColdGas > cold_crit && tdyn > 0.0)
          strdotfull = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
      
      H2sum = 0.0;
      for(i=N_BINS-1; i>=0; i--) H2sum += Gal[p].DiscH2[i];
  }

  //MODULARISE
    feedback(p, centralgal, dt, step, NewStars, NewStarsMetals, stars_sum, strdotfull);
//MODULARISE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
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
      assert(Gal[p].StellarMass >= (StarsPre + NewStarSum)/1.01 && Gal[p].StellarMass <= (StarsPre + NewStarSum)*1.01);
  }


  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}



void update_from_star_formation(int p, double stars, double metallicity, int i)
{
  // In older SAGE, this updated the gas and stellar components.  It only does the gas component now due to the way in which discs are handled.

  // update gas and metals from star formation
  if(stars>=MIN_STARS_FOR_SN)
  {
      //printf("(above SN) stars = %e\n", stars);
      Gal[p].DiscGas[i] -= (1 - RecycleFraction) * stars;
      Gal[p].DiscGasMetals[i] -= metallicity * (1 - RecycleFraction) * stars;

      if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i])
        printf("update_from_star_formation report -- gas metals, gas mass = %e, %e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);
      
      Gal[p].ColdGas -= (1 - RecycleFraction) * stars;
      Gal[p].MetalsColdGas -= metallicity * (1 - RecycleFraction) * stars;
  }
  else
  {
      //printf("(below SN) stars = %e\n", stars);
      Gal[p].DiscGas[i] -= stars;
      Gal[p].DiscGasMetals[i] -= metallicity * stars;
      Gal[p].ColdGas -= stars;
      Gal[p].MetalsColdGas -= metallicity * stars;
  }
    
    
  if(Gal[p].DiscGas[i] <= 0.0){
	//printf("DiscGas in update_SF...%e\n", Gal[p].DiscGas[i]);
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;}

  if(Gal[p].ColdGas <= 0.0){
	//printf("DiscGas in update_SF...%e\n", Gal[p].DiscGas[i]);
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;}
  
}



void update_from_feedback(int p, int centralgal, double reheated_mass, double metallicity, int i)
{
  // Check first just to be sure
    if(!(reheated_mass <= Gal[p].DiscGas[i])) printf("disc gas, reheated gas = %e, %e\n", Gal[p].DiscGas[i], reheated_mass);
  assert(reheated_mass <= Gal[p].DiscGas[i]);
  assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);

  if(SupernovaRecipeOn>0)
  {
    Gal[p].ColdGas -= reheated_mass;
    Gal[p].DiscGas[i] -= reheated_mass;
    
    if(HeatedToCentral)
      Gal[centralgal].HotGas += reheated_mass;
    else
      Gal[p].HotGas += reheated_mass;

      Gal[p].EjectedSNGasMass += reheated_mass;

	if(Gal[p].DiscGas[i]>0.0)
	{
	  Gal[p].MetalsColdGas -= metallicity * reheated_mass;
      Gal[p].DiscGasMetals[i] -= metallicity * reheated_mass;
      if(HeatedToCentral)
        Gal[centralgal].MetalsHotGas += metallicity * reheated_mass;
      else
        Gal[p].MetalsHotGas += metallicity * reheated_mass;
	  assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
	}
	else
	{
	  Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[i];
	  Gal[p].DiscGasMetals[i] = 0.0;
      if(HeatedToCentral)
        Gal[centralgal].MetalsHotGas += Gal[p].DiscGasMetals[i];
      else
        Gal[p].MetalsHotGas += Gal[p].DiscGasMetals[i];
    }
  }

  assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);

  if(Gal[p].DiscGas[i] <= 0.0)
  {
    //printf("DiscGas in update_feedback...%e\n", Gal[p].DiscGas[i]);
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;
  }

  if(Gal[p].ColdGas <= 0.0)
  {
    //printf("ColdGas in update_feedback...%d, %e\n", i, Gal[p].ColdGas);
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;
  }

  if(Gal[p].HotGas <= 0.0)
  {
	//printf("HotGas in update_feedback...%e\n", Gal[p].HotGas);
	Gal[p].HotGas=0.0;
	Gal[p].MetalsHotGas=0.0;
  }
    
  Gal[p].OutflowRate += reheated_mass;
    
}


void update_from_ejection(int p, int centralgal, double ejected_mass)
{
    double metallicityHot;
    
    assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
    assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
    
    if(HeatedToCentral)
    {
    	metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
    	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    
        if(ejected_mass >= Gal[centralgal].HotGas)
        {
            Gal[centralgal].EjectedMass += Gal[centralgal].HotGas;
            Gal[centralgal].HotGas = 0.0;
            Gal[centralgal].MetalsEjectedMass += Gal[centralgal].MetalsHotGas;
            Gal[centralgal].MetalsHotGas = 0.0;
        }
        else if(ejected_mass>0 && Gal[centralgal].HotGas>0.0)
    	{
            Gal[centralgal].HotGas -= ejected_mass;
            Gal[centralgal].EjectedMass += ejected_mass;
            Gal[centralgal].MetalsHotGas -= metallicityHot * ejected_mass;
            Gal[centralgal].MetalsEjectedMass += metallicityHot * ejected_mass;
    	}

    	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    }
    else
    {
        if(ejected_mass >= Gal[p].HotGas)
        {
            Gal[p].EjectedMass += Gal[p].HotGas;
            Gal[p].MetalsEjectedMass += Gal[p].MetalsHotGas;
            Gal[p].HotGas = 0.0;
            Gal[p].MetalsHotGas = 0.0;
        }
        else if(ejected_mass>0 && Gal[p].HotGas>0)
        {
            metallicityHot = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
            Gal[p].EjectedMass += ejected_mass;
            Gal[p].MetalsEjectedMass += ejected_mass * metallicityHot;
            Gal[p].HotGas -= ejected_mass;
            Gal[p].MetalsHotGas -= ejected_mass * metallicityHot;
        }
        
    }
    
    //printf("Eject, Metals = %e, %e\n", Gal[p].EjectedMass, Gal[p].MetalsEjectedMass);
    assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
    assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
}


void combine_stellar_discs(int p, double NewStars[N_BINS], double NewStarsMetals[N_BINS])
{
	double sdisc_spin_mag, J_sdisc, J_new, J_retro, J_sum, cos_angle_sdisc_comb, cos_angle_new_comb, DiscStarSum;
	double SDiscNewSpin[3];
	double Disc1[N_BINS], Disc1Metals[N_BINS], Disc2[N_BINS], Disc2Metals[N_BINS];
	int i;
	
    //printf("disc stars from combine discs 1\n");
    DiscStarSum = get_disc_stars(p);
    
	for(i=0; i<N_BINS; i++)
    {
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
		assert(NewStarsMetals[i] <= NewStars[i]);
    }
	
	// Try not to get confused, where "new" here implies the newly formed stars.  In the cooling recipe, "new" meant the combined disc, here instead denoted "comb".
	
	J_new = 0.0;
	for(i=N_BINS-1; i>=0; i--)
		J_new += NewStars[i] * (DiscBinEdge[i]+DiscBinEdge[i+1])/2.0; // The assumption that DiscBinEdge is now proportional to radius has broken down
	
	// Determine projection angles for combining discs
	if(Gal[p].StellarMass > Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass)
	{
		// Ensure the stellar disc spin magnitude is normalised
		sdisc_spin_mag = sqrt(sqr(Gal[p].SpinStars[0]) + sqr(Gal[p].SpinStars[1]) + sqr(Gal[p].SpinStars[2]));
		assert(sdisc_spin_mag==sdisc_spin_mag);
		if(sdisc_spin_mag>0.0)
		{
			for(i=0; i<3; i++)
			{
				Gal[p].SpinStars[i] /= sdisc_spin_mag; 
                assert(Gal[p].SpinStars[i]==Gal[p].SpinStars[i]);
			}
		}
	
		J_sdisc = get_disc_ang_mom(p, 1);
	
		// Obtain new spin vector of stellar disc
        for(i=0; i<3; i++){
			SDiscNewSpin[i] = Gal[p].SpinGas[i]*J_new + Gal[p].SpinStars[i]*J_sdisc;
            assert(SDiscNewSpin[i]==SDiscNewSpin[i]);}
        
		// Normalise the new spin
		sdisc_spin_mag = sqrt(sqr(SDiscNewSpin[0]) + sqr(SDiscNewSpin[1]) + sqr(SDiscNewSpin[2]));
        if(sdisc_spin_mag>0.0)
        {
            for(i=0; i<3; i++)
                SDiscNewSpin[i] /= sdisc_spin_mag;
        }
		
        sdisc_spin_mag = sqrt(sqr(SDiscNewSpin[0]) + sqr(SDiscNewSpin[1]) + sqr(SDiscNewSpin[2]));
        if(sdisc_spin_mag<0.99 || sdisc_spin_mag>1.01)
        {
            printf("SpinStars somehow became %e\n", sdisc_spin_mag);
            printf("with J_sdisc, J_new = %e, %e\n", J_sdisc, J_new);
            printf("DiscStars, DiscGas = %e, %e\n", get_disc_stars(p), get_disc_gas(p));
        }
        assert(sdisc_spin_mag >= 0.99 && sdisc_spin_mag <= 1.01);
        
		cos_angle_sdisc_comb = Gal[p].SpinStars[0]*SDiscNewSpin[0] + Gal[p].SpinStars[1]*SDiscNewSpin[1] + Gal[p].SpinStars[2]*SDiscNewSpin[2];
		cos_angle_new_comb = Gal[p].SpinGas[0]*SDiscNewSpin[0] + Gal[p].SpinGas[1]*SDiscNewSpin[1] + Gal[p].SpinGas[2]*SDiscNewSpin[2];
	}
	else
	{
		cos_angle_sdisc_comb = 1.0;
		cos_angle_new_comb = 1.0;
		J_sdisc = 0.0;
        for(i=0; i<3; i++)
            SDiscNewSpin[i] = Gal[p].SpinGas[i];
	}
	
	// Combine the discs
	if(cos_angle_sdisc_comb<1.0)
    {
        //printf("projecting both\n");
		project_disc(Gal[p].DiscStars, cos_angle_sdisc_comb, p, Disc1);
		project_disc(Gal[p].DiscStarsMetals, cos_angle_sdisc_comb, p, Disc1Metals);
		project_disc(NewStars, cos_angle_new_comb, p, Disc2);
		project_disc(NewStarsMetals, cos_angle_new_comb, p, Disc2Metals);
		
        Gal[p].StellarMass = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
        Gal[p].MetalsStellarMass = Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass;
		for(i=N_BINS-1; i>=0; i--)
		{
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			Gal[p].DiscStars[i] = Disc1[i] + Disc2[i];
			Gal[p].DiscStarsMetals[i] = Disc1Metals[i] + Disc2Metals[i];
			if(Gal[p].DiscStars[i]==0.0 && Gal[p].DiscStarsMetals[i] < 1e-20) Gal[p].DiscStarsMetals[i] = 0.0;
            Gal[p].StellarMass += Gal[p].DiscStars[i];
            Gal[p].MetalsStellarMass += Gal[p].DiscStarsMetals[i];
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
    else
	{
        //printf("adding directly with cos_angle_sdisc_comb = %e\n", cos_angle_sdisc_comb);
        DiscStarSum = get_disc_stars(p);
        double NewStarSum = 0.0;
        for(i=N_BINS-1; i>=0; i--) NewStarSum += NewStars[i];
        //printf("DiscStarSum, StellarMass, NewStarSum = %e, %e, %e\n", DiscStarSum, Gal[p].StellarMass, NewStarSum);
		for(i=N_BINS-1; i>=0; i--)
		{
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
			Gal[p].DiscStars[i] += NewStars[i];
			Gal[p].DiscStarsMetals[i] += NewStarsMetals[i];
            Gal[p].StellarMass += NewStars[i];
            Gal[p].MetalsStellarMass += NewStarsMetals[i];
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
	
    //printf("disc stars from combine discs 2\n");
    DiscStarSum = get_disc_stars(p);
    
	// Readjust disc to deal with any retrograde stars
	if(cos_angle_sdisc_comb<0.0)
		J_retro = J_sdisc*fabs(cos_angle_sdisc_comb);
	else if(cos_angle_new_comb<0.0)
		J_retro = J_new*fabs(cos_angle_new_comb);
	else
		J_retro = 0.0;
	J_sum = J_sdisc*fabs(cos_angle_sdisc_comb) + J_new*fabs(cos_angle_new_comb);
		
	if(J_retro>0.0)
	{
        //printf("Dealing with retrograde stars\n");
		project_disc(Gal[p].DiscStars, (J_sum - 2.0*J_retro)/J_sum, p, Disc1);
		project_disc(Gal[p].DiscStarsMetals, (J_sum - 2.0*J_retro)/J_sum, p, Disc1Metals);
		for(i=0; i<N_BINS; i++)
		{
			Gal[p].DiscStars[i] = Disc1[i];
			Gal[p].DiscStarsMetals[i] = Disc1Metals[i];
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		}
	}
	
	// Set the new spin direction of the stellar disc
	for(i=0; i<3; i++){
		Gal[p].SpinStars[i] = SDiscNewSpin[i];
		assert(Gal[p].SpinStars[i]==Gal[p].SpinStars[i]);}
    
    //printf("disc stars from combine discs 3\n");
    DiscStarSum = get_disc_stars(p);
    if(DiscStarSum > 1.01*(Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass) || DiscStarSum < (Gal[p].StellarMass-Gal[p].SecularBulgeMass-Gal[p].ClassicalBulgeMass)/1.01)
        printf("Stellar Disc, bulge, total = %e, %e, %e\n", DiscStarSum, Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass, Gal[p].StellarMass);
    
    if(DiscStarSum>0.0) assert(DiscStarSum+Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass <= 1.01*Gal[p].StellarMass && DiscStarSum+Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass >= Gal[p].StellarMass/1.01);

	for(i=0; i<N_BINS; i++){
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);}

    update_stellardisc_scaleradius(p);
}


void project_disc(double DiscMass[N_BINS], double cos_angle, int p, double *NewDisc)
{
	double high_bound, ratio_last_bin;
	int i, j, j_old, k;
    
	cos_angle = fabs(cos_angle); // This function will not deal with retrograde motion so needs an angle less than pi/2
	
	j_old = 0;

	for(i=0; i<N_BINS; i++)
	{
		high_bound = DiscBinEdge[i+1] / cos_angle;
		j = j_old;
		
		while(DiscBinEdge[j]<=high_bound)
		{
			j++;
			if(j==N_BINS) break;
		} 
		j -= 1;
		
		NewDisc[i] = 0.0;
		for(k=j_old; k<j; k++) 
		{
			NewDisc[i] += DiscMass[k];
			DiscMass[k] = 0.0;
		}
		if(i!=N_BINS-1)
		{
			if(j!=N_BINS-1){
				ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]));
				assert(ratio_last_bin<=1.0);}
			else if(high_bound < Gal[p].Rvir/Gal[p].Vvir){
				ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir/Gal[p].Vvir-DiscBinEdge[j]));
				assert(ratio_last_bin<=1.0);}
			else
				ratio_last_bin = 1.0;
			NewDisc[i] += ratio_last_bin * DiscMass[j];
			DiscMass[j] -= ratio_last_bin * DiscMass[j];
		}
		else
		{
			NewDisc[i] = DiscMass[i];
		}
		assert(NewDisc[i]>=0.0);

		j_old = j;
	}
}


void update_HI_H2(int p)
{
    double area, f_H2, f_H2_HI, Pressure, f_sigma;
    int i;
    double angle = acos(Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2])*180.0/M_PI;
    double P_0 = 5.93e-12 / UnitMass_in_g * UnitLength_in_cm * UnitTime_in_s * UnitTime_in_s;
    
    if(Gal[p].Vvir>0.0)
    {
        for(i=0; i<N_BINS; i++)
        {
            area = M_PI * (sqr(Gal[p].DiscRadii[i+1]) - sqr(Gal[p].DiscRadii[i]));
            
            if(H2prescription==1)
            {
                double s, Zp, chi, c_f, Sigma_comp0, Tau_c;
                Zp = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]) / 0.02; // Might also want solar metal fraction to be variable too
                
                if(Zp>0.01 && Zp<1) // Fu et al. 2013
                    c_f = ClumpFactor*pow(Zp, -ClumpExponent);
                else if(Zp>=1)
                    c_f = ClumpFactor;
                else // prescription from Fu not defined here, so assume clumping factor can't exceed this value ~25
                    c_f = ClumpFactor*pow(0.01, -ClumpExponent);
                
                Sigma_comp0 = c_f * Gal[p].DiscGas[i]/area; // see Krumholz & Dekel 2012, originally comes from McKee & Krumholz 2010
                Tau_c = 320 * Zp * Sigma_comp0 * UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm * Hubble_h;
                chi = 3.1 * (1+ 3.1*pow(Zp,0.365)) / 4.1;
                s = log(1 + 0.6*chi + 0.01*chi*chi) / (0.6*Tau_c);
                if(s<2)
                {
                    f_H2 = 1.0 - 0.75*s/(1+0.25*s); // Not actual H2/cold, but rather H2/(H2+HI)
                    f_H2_HI = 1.0 / (1.0/f_H2 - 1.0);
                }
                else
                    f_H2_HI = 0.0;
                
            }
            else
            {
                if(angle <= ThetaThresh)
                {
                    f_sigma =  1.1e6/UnitVelocity_in_cm_per_s / (0.5*Gal[p].Vvir*exp(-(Gal[p].DiscRadii[i]+Gal[p].DiscRadii[i+1])/4.0/Gal[p].StellarDiscScaleRadius)); // Ratio of gas vel dispersion to stars', assuming gas is always 11 km/s
                    Pressure = 0.5*M_PI*G * Gal[p].DiscGas[i] * (Gal[p].DiscGas[i] + f_sigma*Gal[p].DiscStars[i]) / sqr(area);
                }
                else
                    Pressure = 0.5*M_PI*G * sqr(Gal[p].DiscGas[i]/area);
                f_H2_HI = H2FractionFactor * pow(Pressure/P_0, H2FractionExponent);

            }

            
            if(f_H2_HI > 0.0)
            {
                assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
                f_H2 = 0.75 * 1.0/(1.0/f_H2_HI + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i]) / 1.3; //Changes f_H2 from being H2/HI to H2/Cold Gas
                //if(Gal[p].Type==1) f_H2 *= 1.3;
                Gal[p].DiscH2[i] = f_H2 * Gal[p].DiscGas[i];
                Gal[p].DiscHI[i] = Gal[p].DiscH2[i] / f_H2_HI;
            }
            else
            {
                Gal[p].DiscH2[i] = 0.0;
                Gal[p].DiscHI[i] = 0.75*(Gal[p].DiscGas[i]-Gal[p].DiscGasMetals[i])/1.3; // All properly cold hydrogen must be in the form of HI if there's no H2.
                //if(Gal[p].Type==1) Gal[p].DiscHI[i] *= 1.3;
            }
        }
    }
}



