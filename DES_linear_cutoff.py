import numpy as np
import matplotlib.pyplot as plt
import pylab as pil

# Utilities Libraries
from argparse import ArgumentParser

## Setting up argument parser
parser = ArgumentParser(description='cutoff data based on ', usage='python DE_linear_cut.py -t threshold')

parser.add_argument('-t', '--threshold',  type = float, help='chi2 threshold')
parser.add_argument('-c', '--cut_type',  type = str, help='name of the cut')
parsed = parser.parse_args()

threshold = parsed.threshold
cut_type = parsed.cut_type

# reading the files:
invcov      = np.loadtxt('./data/DES_inv_cov.dat')               # inverse DES covariance matrix
dat         = np.loadtxt('./data/DES_vec_data.dat')              # measurements
dl          = np.loadtxt('./data/DES_theory_vec_linear.dat')     # linear theory predictions
dnl         = np.loadtxt('./data/DES_theory_vec_nonlinear.dat')  # nonlinear theory predictions
used_items  = np.loadtxt('./data/DES_used_items.dat')            # used items for the original dataset
thetas      = np.loadtxt('./data/DES_1YR_final_theta_bins.dat')  # angular bins

NLlike  = (dnl - dat).T @ invcov @ (dnl - dat)
Llike   = (dl - dat).T @ invcov @ (dl - dat)
print('NL, L likes:', NLlike, Llike)

cov = np.linalg.inv(invcov)

#------------------------------------------------------------
#> function to remove the data
def remove_data(threshold, data, tnl, tl, invcov):
    # 
    # this functions removes a data-point at a time from DES dataset based on its contributuion on the chi2
    # INPUT ARGUMENTS:
    #   threshold: the difference in chi2 between linear and nonlinear theory that we want to achieve
    #   data    : DES measurements
    #   tnl     : theoretical nonlinear predictions
    #   tl      : theoretical linear predictions
    #   invcov  : inverse covariance matrix 
    #
    # RETURNS: 
    #   ind_to_remove   : list of indexes to remove from the data vector
    #   chi_evolv       : list of chisquares at each step (for plotting)
    
    D = tnl - tl
    delta_chi2 = D.T @invcov @ D
       
    print('inital       delta_chi2:', delta_chi2)
    print('objective    delta_chi2:', threshold)
    
    delta_chi2_now  = delta_chi2
    ind_to_remove   = []
    chi2_evol       = []
    
    chi2_evol.append(delta_chi2_now)
    k=0
    while(delta_chi2_now > threshold and k <  data.shape[0]+1):

        contr = []  # list of contribution of each data point
        for i in range(data.shape[0]):
            # remove the data point i from the data and calculate its contribution
            D_trial = D.copy()
            D_trial[i] = 0.0
            delta_chi2_i = D_trial.T @ (invcov @ D_trial)
            
            # always decreasing chi2: this is not correct
            #if delta_chi2_i < delta_chi2_now:
            #    contr.append(np.abs(delta_chi2_i-delta_chi2_now))
            #else:
            #    contr.append(0)
                
            contr.append(np.abs(delta_chi2_i-delta_chi2_now))

          
        # remove the index with maximum contribution to delta chi^2
        j = contr.index(np.max(contr))
               
        # collect all the indices to remove
        ind_to_remove.append(j)
        
        # remove the data point
        D[j] = 0.0
        
        # recalculate the new delta chi^2
        delta_chi2_now = np.abs(D.T @ invcov @D)
           
        # collect the improvement on delta chi^2
        chi2_evol.append(delta_chi2_now)

        k+=1
        
    return ind_to_remove, chi2_evol


# now remove the data
remove_index, chi2_steps = remove_data(threshold=threshold, data=dat,invcov=invcov, tnl=dnl,tl=dl)
print('number of datapoints removed:', len(remove_index))

# plot the delta chi^2 improvement
plt.plot(range(len(chi2_steps)), chi2_steps)
plt.ylabel(r'$\Delta \chi^2$')
plt.xlabel('step')
plt.axhline(threshold, color = 'black', linestyle = '--', linewidth = 0.75)
pil.savefig('./pdf/chi2_improvement_'+cut_type+'.pdf', bbox_inches='tight')
pil.savefig('./img/chi2_improvement_'+cut_type+'.png', bbox_inches='tight')
plt.clf()


# now delete the datapoints
to_use_items = np.delete(used_items, remove_index, axis = 0)

# open file and write headerc (for CosmoMC)
file = open('./output/DES_1YR_final_'+cut_type+'_cut.dat', 'w')
file.write('#  type bin1 bin2 theta_min theta_max \n')

cutoff_array_max = np.zeros((4,5,5))
cutoff_array_min = np.zeros((4,5,5))
cutoff_array_max.fill(250)
cutoff_array_min.fill(250)

# XIP
xip = to_use_items[ np.where((to_use_items[:,0] == 1.0))]
#replace the bin number with angles
for i in range(xip.shape[0]):
    xip[i,3] = thetas[int(xip[i,3])-1]
    
# find theta min and max for each bin
for i in range(1,5):
    for j in range(1,5):
        #print(i,j)
        xip_this = xip[np.where((xip[:,1] == i) &  (xip[:,2] == j)) ]
        if xip_this.shape[0] > 0:
            theta_max = xip_this[:,3].max()
            theta_min = xip_this[:,3].min()
            cutoff_array_max[0,i-1,j-1] = theta_max
            cutoff_array_min[0,i-1,j-1] = theta_min
            print('xip  ',i,j,np.floor(theta_min), np.ceil(theta_max))
            file.write('xip  '+str(i)+'  '+str(j)+'  '+str(np.floor(theta_min)) + '  '+str(np.ceil(theta_max))+'\n' )

# XIM
xim = to_use_items[ np.where((to_use_items[:,0] == 2.0))]
for i in range(xim.shape[0]):
    xim[i,3] = thetas[int(xim[i,3])-1]
    
for i in range(1,5):
    for j in range(1,5):
        #print(i,j)
        xim_this = xim[np.where((xim[:,1] == i) &  (xim[:,2] == j)) ]
        if xim_this.shape[0]>0:
            theta_max = xim_this[:,3].max()
            theta_min = xim_this[:,3].min()
            cutoff_array_max[1,i-1,j-1] = theta_max
            cutoff_array_min[1,i-1,j-1] = theta_min
            print('xim  ',i,j,np.floor(theta_min), np.ceil(theta_max))
            file.write('xim  '+str(i)+'  '+str(j)+'  '+str(np.floor(theta_min)) + '  '+str(np.ceil(theta_max))+'\n' )

# Gamma_t
gammat = to_use_items[ np.where((to_use_items[:,0] == 3.0))]

for i in range(gammat.shape[0]):
    gammat[i,3] = thetas[int(gammat[i,3])-1]
    
for i in range(1,6):
    for j in range(1,5):
        #print(i,j)
        gammat_this = gammat[np.where((gammat[:,1] == i) &  (gammat[:,2] == j)) ]
        if gammat_this.shape[0]>0:
            theta_max = gammat_this[:,3].max()
            theta_min = gammat_this[:,3].min()
            cutoff_array_max[2,i-1,j-1] = theta_max
            cutoff_array_min[2,i-1,j-1] = theta_min
            print('gammat  ',i,j,np.floor(theta_min), np.ceil(theta_max))
            file.write('gammat  '+str(i)+'  '+str(j)+'  '+str(np.floor(theta_min)) + '  '+str(np.ceil(theta_max))+'\n' )

# w(theta)
wtheta = to_use_items[ np.where((to_use_items[:,0] == 4.0))]
for i in range(wtheta.shape[0]):
    wtheta[i,3] = thetas[int(wtheta[i,3])-1]

for i in range(1,6):
    for j in range(1,6):
        #print(i,j)
        wtheta_this = wtheta[np.where((wtheta[:,1] == i) &  (wtheta[:,2] == j)) ]
        if wtheta_this.shape[0]>0:
            theta_max = wtheta_this[:,3].max()
            theta_min = wtheta_this[:,3].min()
            cutoff_array_max[3,i-1,j-1] = theta_max
            cutoff_array_min[3,i-1,j-1] = theta_min
            print('wtheta  ',i,j,np.floor(theta_min), np.ceil(theta_max))
            file.write('wtheta  '+str(i)+'  '+str(j)+'  '+str(np.floor(theta_min)) + '  '+str(np.ceil(theta_max))+'\n' )


file.close()

props = dict(boxstyle='round', facecolor='white', alpha=0.5)


### producing the plots
for meas_type in range(1,5):
    
    # get the data for this kinf of measurement, xip, xim, gammat or w
    this_meas = used_items[np.where((used_items[:,0] == int(meas_type)))]
    
    # get the angles for this data
    this_meas_thetas_bins = this_meas[:,3].copy()
    for n in range(this_meas_thetas_bins.shape[0]):
        this_meas_thetas_bins[n] = thetas[int(this_meas_thetas_bins[n])-1]
        
    # get the values of the 2point correlation functions
    this_meas_dat = dat[np.where((used_items[:,0] == int(meas_type)))]
    
    # compute the max-min of correlation * theta (useful for y limits)
    min_dat = np.multiply(this_meas_dat,this_meas_thetas_bins ).min()
    max_dat = np.multiply(this_meas_dat,this_meas_thetas_bins ).max()
    
    # extract the bins present for this measurement
    bin1_list = np.unique(this_meas[:,1])
    bin2_list = np.unique(this_meas[:,2])
    
    # create the canvas for the plot (divided in subplots for every bins combination)
    if (meas_type != 4):
        f, axarr = plt.subplots(bin2_list.shape[0], bin1_list.shape[0], sharex = True,  figsize=(20, 15))
        f.subplots_adjust( hspace=0,wspace = 0 )
    else:
        f, axarr = plt.subplots(1, bin2_list.shape[0], sharex = True,   figsize= (20,3))
        f.subplots_adjust( hspace=0,wspace = 0 )
        
    # loop over the bins to fill the subplots
    for b1 in bin1_list:
        for b2 in bin2_list:
            # data for this bins
            this_meas_bins = this_meas[np.where((this_meas[:,1] == b1) & (this_meas[:,2] == b2) )]
            
            # angles for this bins
            this_thetas_bins = this_meas_bins[:,3].copy()
            this_thetas = []
            # compute the theta angles
            for tbin in this_thetas_bins:
                this_thetas.append(thetas[int(tbin) - 1])
                
            # get the indices to extract the theory and data of this bins (and multiply by angle)
            this_ind = np.where((used_items[:,0]==int(meas_type)) & (used_items[:,1] == b1) & (used_items[:,2] == b2))
            this_dat = dat[np.where((used_items[:,0]==int(meas_type)) & (used_items[:,1] == b1) & (used_items[:,2] == b2))]
            this_dat = np.multiply(this_dat, np.array(this_thetas))
            this_tnl = dnl[np.where((used_items[:,0]==int(meas_type)) & (used_items[:,1] == b1) & (used_items[:,2] == b2))]
            this_tnl = np.multiply(this_tnl, np.array(this_thetas))
            this_tl  = dl[np.where((used_items[:,0]==int(meas_type)) & (used_items[:,1] == b1) & (used_items[:,2] == b2))]
            this_tl = np.multiply(this_tl, np.array(this_thetas))
            
            # compute the errorbars
            this_err = []
            for i in range(len(this_ind[0])):
                this_err.append(2.0*np.sqrt(cov[int(this_ind[0][i]), int(this_ind[0][i])]))
               
            # error * angles 
            this_err = np.multiply(this_err, np.array(this_thetas))
            
            # no data in this bins choice
            if len(this_thetas) == 0 and meas_type != 4:
                axarr[int(b2)-1,int(b1)-1].set_xlim(2.7,250)
                axarr[int(b2)-1,int(b1)-1].set_visible(False)
                
            # there is data in this bins
            if len(this_thetas) > 0 and meas_type != 4 :
                # plot theory and data (each multiplied by angle)
                axarr[int(b2)-1,int(b1)-1].semilogx(this_thetas, this_tnl, label = 'theory nonlinear')
                axarr[int(b2)-1,int(b1)-1].semilogx(this_thetas, this_tl, label = 'theory linear')
                axarr[int(b2)-1,int(b1)-1].errorbar(this_thetas, this_dat, yerr=this_err, marker = 'o', linestyle = '', color = 'black')
                
                # set the x-y lims 
                axarr[int(b2)-1,int(b1)-1].set_xlim(np.floor(np.min(this_thetas)),250)
                axarr[int(b2)-1,int(b1)-1].set_ylim(-max_dat*0.5, max_dat*1.5)
                axarr[int(b2)-1,int(b1)-1].set_xscale('log')
                
                # set the bins labels
                top = axarr[int(b2)-1,int(b1)-1].get_ylim()[1]*0.95
                right = axarr[int(b2)-1,int(b1)-1].get_xlim()[1]*0.95
                axarr[int(b2)-1,int(b1)-1].text(0.95, 0.95 , (str(int(b2))+str(int(b1))),
                                            horizontalalignment='right',
                                            verticalalignment='top', bbox=props,
                                            transform=axarr[int(b2)-1,int(b1)-1].transAxes)
                
                # shaded regions representing the cuts
                axarr[int(b2)-1,int(b1)-1].axvspan(1, cutoff_array_min[meas_type-1, int(b1)-1, int(b2)-1]*0.9, alpha=0.5, color='grey', label='cutoff')
                axarr[int(b2)-1,int(b1)-1].axvspan(cutoff_array_max[meas_type-1, int(b1)-1, int(b2)-1]*1.1,300, alpha=0.5, color='grey')
                
                # remove yticks on the subplots that are not on the left.
                # for the subplots on the left set also the labels
                if (b1 > 1):
                    axarr[int(b2)-1,int(b1)-1].set_yticks([])
                else:
                    if (meas_type == 1 ):
                        axarr[int(b2)-1,int(b1)-1].set_ylabel(r'$\theta \, \xi_+$  (arcmin)', fontsize = 15)
                    elif (meas_type == 2 ):
                        axarr[int(b2)-1,int(b1)-1].set_ylabel(r'$\theta \, \xi_-$  (arcmin)', fontsize = 15)
                    elif (meas_type == 3 ):
                        axarr[int(b2)-1,int(b1)-1].set_ylabel(r'$\theta \, \gamma_t$  (arcmin)', fontsize = 15)
                  
                # remove xticks on the subplots that are not at the bottom.
                # for the subplots at the bottom set also the labels
                if ((b2 < bin2_list.max()) and (meas_type != 4)):
                    axarr[int(b2)-1,int(b1)-1].set_xticks([])
                else:
                    axarr[int(b2)-1,int(b1)-1].set_xlabel(r'$\theta$ (arcmin)', fontsize = 15)
                    
                if (b1 == 1 and b2 ==1):
                    axarr[0,0].legend(loc='upper left', fontsize = 12)
                    
            # if w(theta) the plot style is different
            elif len(this_thetas) > 0 and meas_type == 4:
                # plot theory and data
                axarr[int(b1)-1].semilogx(this_thetas, this_tnl, label = 'theory nonlinear')
                axarr[int(b1)-1].semilogx(this_thetas, this_tl, label = 'theory linear')
                axarr[int(b1)-1].errorbar(this_thetas, this_dat, yerr=this_err, marker = 'o', linestyle = '', color = 'black')

                # set x-y limits
                axarr[int(b1)-1].set_xlim(np.floor(np.min(this_thetas)),250)
                axarr[int(b1)-1].set_ylim(-max_dat*0.5, max_dat*1.5)
                axarr[int(b1)-1].set_xscale('log')
                
                # shaded regions representing the cuts
                axarr[int(b1)-1].axvspan(1, cutoff_array_min[meas_type-1, int(b1)-1, int(b2)-1]*0.9, alpha=0.5, color='grey')
                axarr[int(b1)-1].axvspan(cutoff_array_max[meas_type-1, int(b1)-1, int(b2)-1]*1.1,300, alpha=0.5, color='grey')
              
                # bins labels
                top = axarr[int(b1)-1].get_ylim()[1]*0.95
                right = axarr[int(b1)-1].get_xlim()[1]*0.95
                axarr[int(b1)-1].text(0.95, 0.95, (str(int(b1))+str(int(b2))),
                                    horizontalalignment='right',
                                    verticalalignment='top', bbox=props, 
                                    transform=axarr[int(b1)-1].transAxes )
                
                # the xlabels are set for each plot here
                axarr[int(b1)-1].set_xlabel(r'$\theta$ (arcmin)', fontsize = 15)
                
                # remove yticks on the subplots that are not on the left.
                # for the subplots on the left set also the labels
                if (b1 > 1):
                    #print(b1)
                    axarr[int(b1)-1].set_yticks([])
                else:
                    axarr[int(b1)-1].set_ylabel(r'$\theta \, w (\theta) $  (arcmin)' , fontsize = 15)
                if (b1 == 1):
                    axarr[0].legend(loc='lower left', fontsize = 12)
                
    # save the figure and clear the figure now for the next plot    
    pil.savefig('./pdf/m'+str(int(meas_type))+cut_type+'.pdf', bbox_inches='tight')
    pil.savefig('./img/m'+str(int(meas_type))+cut_type+'.png', bbox_inches='tight')
    plt.clf()


