
########## FOE ##########


########## Endstation ##########

## Ion chamber: FMB Oxford I404
from math import exp

def curr_to_flux(amp):
    '''Converts Ion Chamber current [A] to flux [ph/s] for FMB Oxford IC filled with gas N2 at 1 atm'''
    E = getE(q=1)	## E in [keV]
    V_ion = 0.036	## ionization energy of N2 gas in [keV]
    IC_len = 6.0	## active length of Ion Chamber in [cm]
    qe = 1.602e-19	## electron charge in [C]

    ## Absorption length [cm] of gas N2 (1 atm, 1.131 g/L) vs E [keV]
    # based on polynomial fit to the calculated abs length data from: henke.lbl.gov/optical_constants/atten2.html 
    # see /home/xf11bm/masa/atten_len_N2* 
    abs_len = 355.21 - 112.26*E + 11.200*E*E - 0.10611*E*E*E	

    N_abs = amp*V_ion/(qe*E)
    flux = N_abs / (1.0 - exp(-IC_len/abs_len))
    return(flux)    

def get_bim3(q=0):
    '''Returns flux at ion chamber in [ph/s] (q=1 for quiet)'''
    bim3_v1 = caget('XF:11BMB-BI{IM:3}:IC1_MON')
    bim3_v2 = caget('XF:11BMB-BI{IM:3}:IC2_MON')
    bim3_h1 = caget('XF:11BMB-BI{IM:3}:IC3_MON')
    bim3_h2 = caget('XF:11BMB-BI{IM:3}:IC4_MON')
    flux_v = curr_to_flux(bim3_v1+bim3_v2)
    flux_h = curr_to_flux(bim3_h1+bim3_h2)
    if q==0:
        print('BIM3 Ion Chamber (Endstation)')
        print('Vertical:')
        print(' Signal 1: %.3e A' % bim3_v1)
        print(' Signal 2: %.3e A' % bim3_v2)
        print(' Signal - total: %.3e A' % (bim3_v1+bim3_v2))
        print('   Flux - total: %.3e ph/s' % flux_v)
        print('Horizontal:')
        print(' Signal 1: %.3e A' % bim3_h1)
        print(' Signal 2: %.3e A' % bim3_h2)
        print(' Signal - total: %.3e A' % (bim3_h1+bim3_h2))
        print('   Flux - total: %.3e ph/s' % flux_h)
    return(flux_h)

