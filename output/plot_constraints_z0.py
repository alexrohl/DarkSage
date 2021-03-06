# Make plots for z=0 galaxies that are typically used to constrain Dark Sage.

from pylab import *
import os
import routines as r
import random

# Warnings are annoying
import warnings
warnings.filterwarnings("ignore")


###### USER NEEDS TO SET THESE THINGS ######
dir_name = raw_input() #if null then default
indir = 'results/millennium/' # directory where the Dark Sage data are

sim = 0 # which simulation Dark Sage has been run on -- if it's new, you will need to set its defaults below.
#   0 = Mini Millennium, 1 = Full Millennium, 2 = SMDPL

fpre = 'model_z0.000' # what is the prefix name of the z=0 files
files = range(8) # list of file numbers you want to read

Nannuli = 30 # number of annuli used for discs in Dark Sage
FirstBin = 1.0 # copy from parameter file -- sets the annuli's sizes
ExponentBin = 1.4
###### ============================== ######



##### SIMULATION DEFAULTS #####
if sim==0:
    h = 0.73
    Lbox = 62.5/h * (len(files)/8.)**(1./3)
elif sim==1:
    h = 0.73
    Lbox = 500.0/h * (len(files)/512.)**(1./3)
elif sim==2:
    h = 0.6777
    Lbox = 400.0/h * (len(files)/1000.)**(1./3)
# add here 'elif sim==3:' etc for a new simulation
else:
    print 'Please specify a valid simulation.  You may need to add its defaults to this code.'
    quit()
######  ================= #####



##### READ DARK SAGE DATA #####
DiscBinEdge = np.append(0, np.array([FirstBin*ExponentBin**i for i in range(Nannuli)])) / h
G = r.darksage_snap(indir+fpre, files, Nannuli=Nannuli)
######  ================= #####



##### SET PLOTTING DEFAULTS #####
fsize = 26
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman', 'legend.numpoints': 1, 'legend.columnspacing': 1, 'legend.fontsize': fsize-4, 'xtick.top': True, 'ytick.right': True})

NpartMed = 100 # minimum number of particles for finding relevant medians for minima on plots

if dir_name != '':
    outdir = indir+dir_name+'/' # where the plots will be saved
else:
    outdir = indir + 'plots/'


if not os.path.exists(outdir): os.makedirs(outdir)
######  =================== #####





##### PLOT 1: MASS FUNCTIONS #####

# Stellar mass function
fig, ax = plt.subplots(1, 3, sharey=True)
SM = np.array(G['StellarMass']*1e10/h, dtype=np.float32)
BTT = (G['InstabilityBulgeMass'] + G['MergerBulgeMass']) / G['StellarMass']
f_SM = (G['LenMax']==NpartMed)*(SM>0)*np.isfinite(SM)
SM_med = round(np.median(np.log10(SM[f_SM])), 3) if len(f_SM[f_SM])>0 else 8.9

r.massfunction(SM, Lbox, range=[SM_med-0.1, 12.1], ls='--', ax=ax[0], label=r'{\sc Dark Sage}, $N_{\rm p}\!\geq\!20$')
r.massfunction(SM[G['LenMax']>=100], Lbox, range=[SM_med-0.1, 12.1], ax=ax[0], label=r'{\sc Dark Sage}, $N_{\rm p,max}\!\geq\!100$')
r.massfunction(SM[(BTT<=0.5)*(G['LenMax']>=100)], Lbox, range=[SM_med-0.1, 12.1], c='b', lw=1, ax=ax[0])
r.massfunction(SM[BTT<=0.5], Lbox, range=[SM_med-0.1, 12.1], c='b', lw=1, ls='--', ax=ax[0])
r.massfunction(SM[(BTT>0.5)*(G['LenMax']>=100)], Lbox, range=[SM_med-0.1, 12.1], c='r', lw=1, ax=ax[0])
r.massfunction(SM[BTT>0.5], Lbox, range=[SM_med-0.1, 12], c='r', lw=1, ls='--', ax=ax[0])
r.stellar_massfunction_obsdata(h, ax[0])
SMF_bd, logM_bd = r.schechter(3.67e-3*(h/0.7)**3, 10**(10.74)/(h/0.7)**2, -0.525, logM=np.arange(SM_med-0.1, 12.1,0.1))
SMF_dd, logM_dd = r.schechter(0.855e-3*(h/0.7)**3, 10**(10.70)/(h/0.7)**2, -1.39, logM=np.arange(SM_med-0.1, 12.1,0.1))
ax[0].plot(logM_dd, SMF_dd, 'b-', lw=8, alpha=0.2, label=r'Moffett et al.~(2016) disc-dominated')
ax[0].plot(logM_bd, SMF_bd, 'r-', lw=8, alpha=0.2, label=r'Moffett et al.~(2016) bulge-dominated')
ax[0].legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.025, -0.03))
ax[0].set_xlim(SM_med, 12)
ax[0].set_xlabel(r'$\log_{10}(m_*~[{\rm M}_{\odot}])$')
ax[0].set_ylabel(r'$\Phi~[{\rm Mpc}^{-3}~{\rm dex}^{-1}]$')
ax[0].set_xticks([9,9.5,10,10.5,11,11.5])
ax[0].set_yscale('log')


# HI mass function
HIM = np.array(np.sum(G['DiscHI'],axis=1)*1e10/h, dtype=np.float32)
f_HIM = (G['LenMax']==NpartMed) * (HIM>0) * np.isfinite(HIM)
HIM_med = round(np.median(np.log10(HIM[f_HIM])), 3) if len(f_HIM[f_HIM])>0 else 9.05

r.HIH2_massfunction_obsdata(h=h, HI=True, H2=False, Z=True, M=True, ax=ax[1])
HIMF_J18, logM_J18 = r.schechter(4.5e-3*(h/0.7)**3, (10**9.94)*(0.7/h)**2, -1.25, logM=np.arange(HIM_med-0.1, 11.5,0.1))
ax[1].plot(logM_J18, HIMF_J18, '-', color='chocolate', lw=6, alpha=0.4, label=r'Jones et al.~(2018)')
ax[1].legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.025, -0.03))
r.massfunction(HIM, Lbox, range=[HIM_med-0.1, 11.5], ls='--', ax=ax[1])
r.massfunction(HIM[G['LenMax']>=100], Lbox, range=[HIM_med-0.1, 11.5], ax=ax[1])
ax[1].set_xlabel(r'$\log_{10}(m_{\rm H\,{\LARGE {\textsc i}}}~[{\rm M}_{\odot}])$')
ax[1].set_xlim(HIM_med, 10.85)
ax[1].set_ylabel('')
ax[1].set_xticks(np.arange(9.2,10.85,0.3))


# H2 mass function
H2M = np.array(np.sum(G['DiscH2'],axis=1)*1e10/h, dtype=np.float32)
f_H2M = (G['LenMax']==NpartMed) * (H2M>0) * np.isfinite(H2M)
H2M_med = round(np.median(np.log10(H2M[f_H2M])), 3) if len(f_H2M[f_H2M])>0 else 8.5

r.HIH2_massfunction_obsdata(h=h, HI=False, H2=True, K=True, OR=False, B=True, ax=ax[2])
ax[2].legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.025, -0.03))
r.massfunction(H2M, Lbox, range=[H2M_med-0.1, 11.5], ls='--', ax=ax[2])
r.massfunction(H2M[G['LenMax']>=100], Lbox, range=[H2M_med-0.1, 11.5], ax=ax[2])
ax[2].set_xlabel(r'$\log_{10}(m_{\rm H_2}~[{\rm M}_{\odot}])$')
ax[2].set_xlim(H2M_med, 10.8)
ax[2].set_ylim(1e-6,3e-2)
ax[2].set_ylabel('')
ax[2].set_xticks(np.arange(8.8,10.8,0.4))

fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
r.savepng(outdir+'1-MassFunctions', xsize=1700, ysize=512, fig=fig)
##### ====================== #####



##### PLOT 2: HI FRACTION and METALLICITY #####

Nmin = 20
# HI fraction
fig, ax = plt.subplots(2, 1, sharex=True)
logM_B15, logHIfrac_B15, err = r.Brown_HI_fractions(h)
f = (SM>=10**8.93) * (G['LenMax']>=100)
bins, mean_SM, mean_HIfrac = r.meanbins(SM[f], (HIM/SM)[f], 10**logM_B15)
bins=10**np.arange(SM_med-0.1,12.1,0.1)
x_av, y_high, y_med, y_low, y_mean = r.percentiles(SM, HIM/SM, bins=bins, addMean=True, Nmin=Nmin)
ax[0].plot(np.log10(x_av), np.log10(y_mean), 'k--', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p} \! \geq \! 20$')
x_av, y_high, y_med, y_low, y_mean = r.percentiles(SM[G['LenMax']>=100], (HIM/SM)[G['LenMax']>=100], bins=bins, addMean=True, Nmin=Nmin)
ax[0].plot(np.log10(x_av), np.log10(y_mean), 'k-', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p,max} \! \geq \! 100$')
ax[0].plot(logM_B15, logHIfrac_B15, 's', color='purple', ms=10, label=r'Brown et al.~(2015)', alpha=0.5, zorder=3)
ax[0].plot(np.log10(mean_SM), np.log10(mean_HIfrac), 'k*', ms=10, label=r'{\sc Dark Sage}, matched bins', zorder=2)
ax[0].set_ylim(-1.75,0.5)
ax[0].set_yticks(np.arange(-1.6,0.5,0.4))
ax[0].set_ylabel(r'$\log_{10}\langle m_{\rm H\,{\LARGE {\textsc i}}} / m_* \rangle$')
ax[0].legend(loc='upper right', frameon=False)


# Mass--metallicity
lOH = 9 + np.log10(G['MetalsColdGas'] / G['ColdGas'] / 0.02)
SM_mean, lOH_high, lOH_mid, lOH_low = r.percentiles(SM, lOH, bins=bins, Nmin=Nmin)
ax[1].plot(np.log10(SM_mean), lOH_mid, 'k--', lw=3)
ax[1].plot(np.log10(SM_mean), lOH_low, 'k--', lw=1.5)
ax[1].plot(np.log10(SM_mean), lOH_high, 'k--', lw=1.5)
SM_mean, lOH_high, lOH_mid, lOH_low = r.percentiles(SM[G['LenMax']>=100], lOH[G['LenMax']>=100], bins=bins, Nmin=Nmin)
ax[1].plot(np.log10(SM_mean), lOH_mid, 'k-', lw=3, label=r'{\sc Dark Sage} median')
ax[1].plot(np.log10(SM_mean), lOH_low, 'k-', lw=1.5)
ax[1].plot(np.log10(SM_mean), lOH_high, 'k-', lw=1.5, label=r'16$^{\rm th}$ \& 84$^{\rm th}$ \%iles')
x_obs, y_low, y_high = r.Tremonti04(h)
ax[1].fill_between(x_obs, y_high, y_low, color='orchid', alpha=0.4)
ax[1].plot([0,1], [0,1], '-', color='orchid', alpha=0.4, lw=8, label=r'Tremonti et al.~(2004)')
ax[1].axis([SM_med, 11.7, 8, 9.4])
ax[1].set_xlabel(r'$\log_{10}(m_*~[{\rm M}_{\odot}])$')
ax[1].set_ylabel(r'$12 + \log_{10}({\rm O/H})$')
ax[1].set_yticks(np.arange(8.25,9.5,0.25))
ax[1].legend(loc='lower right', frameon=False, ncol=2, bbox_to_anchor=(1.025,0))

fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
r.savepng(outdir+'2-HIfrac_MassMet', xsize=768, ysize=700)
##### =================================== #####




##### PLOT 3: BLACK HOLE -- BULGE MASS #####
plt.clf()
BM = (G['InstabilityBulgeMass'] + G['MergerBulgeMass']) * 1e10/h
BM_med = np.log10(np.median(BM[(G['LenMax']==NpartMed)*(BM>0)]))
BHM = G['BlackHoleMass'] * 1e10/h
bins = 10**np.arange(BM_med, 12.5, 0.1)
BM_mean, BHM_high, BHM_mid, BHM_low, BHM_mean = r.percentiles(BM, BHM, bins=bins, addMean=True, Nmin=Nmin)
floor = 1
BHM_low[BHM_low<=floor] = floor
r.BH_bulge_obs(h)
plt.plot(np.log10(BM_mean), np.log10(BHM_mid), 'k--', lw=3, label=r'{\sc Dark Sage} median, $N_{\rm p}\!\geq\!20$')
plt.plot(np.log10(BM_mean), np.log10(BHM_high), 'k--', lw=1.5)
plt.plot(np.log10(BM_mean), np.log10(BHM_low), 'k--', lw=1.5)
BM_mean, BHM_high, BHM_mid, BHM_low, BHM_mean = r.percentiles(BM[G['LenMax']>=100], BHM[G['LenMax']>=100], bins=bins, addMean=True, Nmin=Nmin)
floor = 1
BHM_low[BHM_low<=floor] = floor
plt.plot(np.log10(BM_mean), np.log10(BHM_mid), 'k-', lw=3, label=r'Median, $N_{\rm p,max}\!\geq\!100$')
plt.plot(np.log10(BM_mean), np.log10(BHM_high), 'k-', lw=1.5)
plt.plot(np.log10(BM_mean), np.log10(BHM_low), 'k-', lw=1.5, label=r'16$^{\rm th}$ \& 84$^{\rm th}$ \%iles')
plt.xlabel(r'$\log_{10}(m_{\rm bulge}~[{\rm M}_{\odot}])$')
plt.ylabel(r'$\log_{10}(m_{\rm BH}~[{\rm M}_{\odot}])$')
plt.axis([BM_med, 12.2, 5.5, 9.8])
plt.legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.05,0.95), ncol=2)
r.savepng(outdir+'3-BHBM', xsize=768, ysize=400)
##### ================================ #####


##### PLOT 4: BARYONIC TULLY-FISHER #####
plt.clf()
BaryM = (G['ColdGas']/1.3 + G['StellarMass']) * 1e10/h
BaryM_med = np.log10(np.median(BaryM[G['LenMax']==NpartMed]))
Vel_Profiles = DiscBinEdge[1:] / (G['DiscRadii'][:,1:]*1e3/h)
Vmax = np.max(Vel_Profiles, axis=1)
filt = (G['ColdGas']/1.3 > G['StellarMass']) * (BaryM>=BaryM_med-0.11) # Stark data only had gas-dominated galaxies
## Code can get slow on these lines if doing a large simulation
V3re = np.array([np.interp(3*G['StellarDiscScaleRadius'][i], G['DiscRadii'][i,:], DiscBinEdge/(G['DiscRadii'][i,:]*1e3/h)) for i in np.where(filt)[0]])
V4re = np.array([np.interp(4*G['StellarDiscScaleRadius'][i], G['DiscRadii'][i,:], DiscBinEdge/(G['DiscRadii'][i,:]*1e3/h)) for i in np.where(filt)[0]])
##
Vend = Vel_Profiles[:,-1][filt]
cond = ((abs(np.log10(V3re/Vend))<=np.log10(1.15)) + (abs(np.log10(V4re/Vend))<=np.log10(1.1)))
filt[np.where(filt)[0][~cond]] = False
BaryM_mean, Vmax_high, Vmax_mid, Vmax_low, Vmax_mean  = r.percentiles(BaryM[filt], Vmax[filt],  bins=10**np.arange(BaryM_med-0.1,12,0.1), addMean=True, Nmin=Nmin)
plt.plot(np.log10(Vmax_mid), np.log10(BaryM_mean), 'k--', lw=3, label=r'{\sc Dark Sage} median, $N_{\rm p}\!\geq\!20$')
plt.plot(np.log10(Vmax_high), np.log10(BaryM_mean), 'k--', lw=1.5)
plt.plot(np.log10(Vmax_low), np.log10(BaryM_mean), 'k--', lw=1.5)
filt = filt * (G['LenMax']>=100)
BaryM_mean, Vmax_high, Vmax_mid, Vmax_low, Vmax_mean  = r.percentiles(BaryM[filt], Vmax[filt],  bins=10**np.arange(BaryM_med-0.1,12,0.1), addMean=True, Nmin=Nmin)
plt.plot(np.log10(Vmax_mid), np.log10(BaryM_mean), 'k-', lw=3, label=r'{\sc Dark Sage} median, $N_{\rm p,max}\!\geq\!100$')
plt.plot(np.log10(Vmax_high), np.log10(BaryM_mean), 'k-', lw=1.5)
plt.plot(np.log10(Vmax_low), np.log10(BaryM_mean), 'k-', lw=1.5, label=r'{\sc Dark Sage} 16$^{\rm th}$ \& 84$^{\rm th}$ \%iles')
#
x_obs = np.linspace(1.25,2.2,100)
y_obs_arr = np.array([[4.01*x_obs + 2.05], [4.01*x_obs + 1.53], [3.87*x_obs + 2.05], [3.87*x_obs + 1.53]]) # Random uncertainty only here
y_obs = 3.94*x_obs + 1.79
y_obs_min = np.min(y_obs_arr, axis=0)[0] + 2*np.log10(0.75/h)
y_obs_max = np.max(y_obs_arr, axis=0)[0] + 2*np.log10(0.75/h)
plt.fill_between(x_obs, y_obs_max, y_obs_min, color='darkmagenta', alpha=0.3)
plt.plot(x_obs, y_obs, '-', color='darkmagenta', lw=2)
plt.plot([0,1], [0,1], '-', color='darkmagenta', lw=8, alpha=0.3, label=r'Stark et al.~(2009)')
plt.legend(loc='upper left', frameon=False)
plt.xlabel(r'$\log_{10}\left(V_{\rm max}~[{\rm km\,s}^{-1}]\right)$')
plt.ylabel(r'$\log_{10}\left(m_* + X^{-1}\,m_{\rm H\,{\LARGE{\textsc i}}+H_2}~[{\rm M}_{\odot}]\right)$')
plt.xticks(np.arange(1.9,2.6,0.1))
plt.axis([1.85, 2.65, BaryM_med, 11.7])
r.savepng(outdir+'4-BTF', xsize=768, ysize=400)
##### ============================= #####



##### PLOT 4: HI and H2 PROFILES #####
filt = (Vmax>=175) * (Vmax<=235) * (G['StellarMass']/h > 1.0) * (G['ColdGas']/h > 10**-0.8) * (G['Type']==0) * (BTT<0.5)
area = np.pi*(G['DiscRadii'][:,1:]**2 - G['DiscRadii'][:,:-1]**2) * 1e12/h**2
R_av = np.sqrt((G['DiscRadii'][:,1:]**2 + G['DiscRadii'][:,:-1]**2) / 2.) * 1e3/h # kpc
Sigma_HI = (G['DiscHI']*1e10/h) / area
Sigma_H2 = (G['DiscH2']*1e10/h) / area

Nprof = 100 # Number of profiles to plot
args = random.sample(np.where(filt)[0], Nprof) if len(filt[filt])>Nprof else np.where(filt)[0]

fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)
for a in args:
    ax[0].plot(R_av[a,:], Sigma_HI[a,:], 'k.-', lw=1, alpha=0.2)
    ax[1].plot(R_av[a,:], Sigma_H2[a,:], 'k.-', lw=1, alpha=0.2)
ax[0].plot([0,1], [0,1], 'k.-', lw=1, alpha=0.5, label=r'{\sc Dark Sage}') # legend
ax[0].set_yscale('log', nonposy='clip')
ax[0].errorbar([0,0], [0,1], [1,1], ecolor='darkcyan', color='darkcyan', label=r'Leroy et al. (2008)')
ax[0].set_ylabel(r'$\Sigma_{\rm H\,{\LARGE{\textsc i}}}$ [M$_{\bigodot}$ pc$^{-2}$]')
ax[0].legend(loc='best', frameon=False)
r.Leroygals(HI=True, HighVvir=True, LowVvir=True, ax=ax[0], h=h, c='darkcyan', alpha=0.8)

ax[1].set_yscale('log', nonposy='clip')
ax[1].set_ylabel(r'$\Sigma_{\rm H_2}$ [M$_{\bigodot}$ pc$^{-2}$]')
ax[1].set_xlabel(r'$r$ [kpc]')
r.Leroygals(H2=True, HighVvir=True, LowVvir=True, ax=ax[1], h=h, c='darkcyan', alpha=0.8)
ax[1].axis([0,25,5e-1,1e3])


fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
r.savepng(outdir+'5-HIH2Surface', xsize=768, ysize=700)
##### ========================== #####

