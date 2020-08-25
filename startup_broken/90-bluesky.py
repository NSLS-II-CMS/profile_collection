def detselect(detector_object, suffix="_stats1_total"):
    """Switch the active detector and set some internal state"""

    if isinstance(detector_object, (list, tuple)):
        gs.DETS = detector_object
        gs.PLOT_Y = detector_object[0].name + suffix
        gs.TABLE_COLS = [gs.PLOT_Y] 

    else:
        gs.DETS =[detector_object]
        gs.PLOT_Y = detector_object.name + suffix
        gs.TABLE_COLS = [gs.PLOT_Y] 



##### I/O devices 
from epics import (caget, caput)

def pneumatic(inout,pv_r,pv_in,pv_out,ss1,quiet):
    if inout == 1:
        caput(pv_in,1)
        ss2 = 'has been inserted'
    elif inout == 0:
        caput(pv_out,1)
        ss2 = 'has been retracted'
    else:
        if caget(pv_r)==1:
            ss2 = 'is IN'
        else:
            ss2 = 'is OUT'
    if quiet==0:
        print(ss1+' '+ss2) 


## Fluorescence screen 1 (FOE)
def io_fs1(inout,q=0):
    pv_r = 'XF:11BMA-BI{FS:1}Pos-Sts'
    pv_in = 'XF:11BMA-BI{FS:1}Cmd:In-Cmd'
    pv_out = 'XF:11BMA-BI{FS:1}Cmd:Out-Cmd'
    ss1 = 'Fluorescence screen 1 (FOE)'
    pneumatic(inout,pv_r,pv_in,pv_out,ss1,q)

## Fluorescence screen 3 (Endstation)
def io_fs3(inout,q=0):
    pv_r = 'XF:11BMB-BI{FS:3}Pos-Sts'
    pv_in = 'XF:11BMB-BI{FS:3}Cmd:In-Cmd'
    pv_out = 'XF:11BMB-BI{FS:3}Cmd:Out-Cmd'
    ss1 = 'Fluorescence screen 3 (Endstation)'
    pneumatic(inout,pv_r,pv_in,pv_out,ss1,q)

## Fluorescence screen 4 (Endstation)
def io_fs4(inout,q=0):
    pv_r = 'XF:11BMB-BI{FS:4}Pos-Sts'
    pv_in = 'XF:11BMB-BI{FS:4}Cmd:In-Cmd'
    pv_out = 'XF:11BMB-BI{FS:4}Cmd:Out-Cmd'
    ss1 = 'Fluorescence screen 3 (Endstation)'
    pneumatic(inout,pv_r,pv_in,pv_out,ss1,q)

## BIM 5 - RIGI (Endstation)
def io_bim5(inout,q=0):
    pv_r = 'XF:11BMB-BI{IM:5}Pos-Sts'
    pv_in = 'XF:11BMB-BI{IM:5}Cmd:In-Cmd'
    pv_out = 'XF:11BMB-BI{IM:5}Cmd:Out-Cmd'
    ss1 = 'BIM 5 - RIGI (Endstation)'
    pneumatic(inout,pv_r,pv_in,pv_out,ss1,q)


## Attenuation filter box
def io_atten(pos,inout,q=0):
    if pos >= 1 and pos <= 8:
        pv_r = 'XF:11BMB-OP{Fltr:' + str(int(pos)) + '}Pos-Sts'        
        pv_in = 'XF:11BMB-OP{Fltr:' + str(int(pos)) + '}Cmd:In-Cmd'        
        pv_out = 'XF:11BMB-OP{Fltr:' + str(int(pos)) + '}Cmd:Out-Cmd'        
        ss1 = 'Atten filter ' + str(int(pos))
        pneumatic(inout,pv_r,pv_in,pv_out,ss1,q)
    else:
        print('Attenuator position must be an integer between 1 and 8')

from math import (exp, log)

def get_atten_trans():
    E = getE(q=1)		# Current E [keV]
    
    if E < 6.0 or E > 18.0:
        print('Transmission data not available at the current X-ray enegy.')

    else:
        N = []
        for i in np.arange(8):
            N.append(caget('XF:11BMB-OP{Fltr:' + str(int(i)+1) + '}Pos-Sts'))

        N_Al = N[0] + 2*N[1] + 4*N[2] + 8*N[3]
        N_Nb = N[4] + 2*N[5] + 4*N[6] + 8*N[7]

        d_Nb = 0.1	# Thickness [mm] of one Nb foil 
        d_Al = 0.25	# Thickness [mm] of one Al foil 

        # Absorption length [mm] based on fits to LBL CXRO data for 6 < E < 19 keV
        l_Nb = 1.4476e-3 - 5.6011e-4 * E + 1.0401e-4 * E*E + 8.7961e-6 * E*E*E
        l_Al = 5.2293e-3 - 1.3491e-3 * E + 1.7833e-4 * E*E + 1.4001e-4 * E*E*E

        # transmission factors
        tr_Nb = exp(-N_Nb*d_Nb/l_Nb) 
        tr_Al = exp(-N_Al*d_Al/l_Al) 
        tr_tot = tr_Nb*tr_Al

        print('%dx 0.25mm Al (%.1e) and %dx 0.10mm Nb (%.1e)' % (N_Al, tr_Al, N_Nb, tr_Nb))
        print('Combined transmission is %.1e' % tr_tot)

        return(tr_tot)

def set_atten_trans(tr):
    E = getE(q=1)		# Current E [keV]
    
    if E < 6.0 or E > 18.0:
        print('Transmission data not available at the current X-ray enegy.')

    elif tr > 1.0 or tr < 1.0e-10:
        print('Requested attenuator transmission is not valid.')

    else:
        d_Nb = 0.1	# Thickness [mm] of one Nb foil 
        d_Al = 0.25	# Thickness [mm] of one Al foil 

        # Absorption length [mm] based on fits to LBL CXRO data for 6 < E < 19 keV
        l_Nb = 1.4476e-3 - 5.6011e-4 * E + 1.0401e-4 * E*E + 8.7961e-6 * E*E*E
        l_Al = 5.2293e-3 - 1.3491e-3 * E + 1.7833e-4 * E*E + 1.4001e-4 * E*E*E
 
        d_l_Nb = d_Nb/l_Nb 
        d_l_Al = d_Al/l_Al 

        # Number of foils to be inserted (picks a set that gives smallest deviation from requested transmission)
        dev=[]
        for i in np.arange(16):
            for j in np.arange(16):
                dev_ij = abs(tr - exp(-i*d_l_Nb)*exp(-j*d_l_Al))
                dev.append(dev_ij)
                if (dev_ij == min(dev)):
                    N_Nb = i			# number of Nb foils selected
                    N_Al = j			# number of Al foils selected

        N=[]
        state = N_Al
        for i in np.arange(4):
            N.append(state % 2)
            state = int(state/2)

        state = N_Nb
        for i in np.arange(4):
            N.append(state % 2)
            state = int(state/2)

        for i in np.arange(8):
            io_atten(i+1,N[i],q=1)

        sleep(1.); return(get_atten_trans())


### 1-stage valves
def single_valve(cmd,pv_r,pv_op,pv_cl,ss1,quiet):
    if cmd=='st':
        st = caget(pv_r)
        if st == 1:
            ss2 = 'valve is open'
        if st == 0:
            ss2 = 'valve is closed'
    if cmd=='o' or cmd=='open':
        caput(pv_op,1)
        ss2 = 'valve has been opened'
    if cmd=='c' or cmd=='close':
        caput(pv_cl,1)
        ss2 = 'valve has been closed'
    if quiet==0:
        print(ss1+' '+ss2) 

### 2-stage valves
def dual_valve(cmd,pv_r_soft,pv_op_soft,pv_cl_soft,pv_r,pv_op,pv_cl,ss1,quiet):
    if cmd=='st':
        st_h = caget(pv_r)
        st_s = caget(pv_r_soft)
        if st_h == 1 and st_s == 1:
            ss2 = 'main open and soft open'
        if st_h == 0 and st_s == 1:
            ss2 = 'main closed and soft open'
        if st_h == 1 and st_s == 0:
            ss2 = 'main open and soft closed'
        if st_h == 0 and st_s == 0:
            ss2 = 'main closed and soft closed'
    if cmd=='o' or cmd=='open':
        caput(pv_cl_soft,1)
        sleep(0.2)
        caput(pv_op,1)
        ss2 = 'main has been opened (soft closed)'
    if cmd=='so' or cmd=='soft':
        caput(pv_cl,1)
        sleep(0.5)
        caput(pv_op_soft,1)
        ss2 = 'soft has been opened (main closed)'
    if cmd=='c' or cmd=='close':
        caput(pv_cl,1)
        sleep(0.2)
        caput(pv_cl_soft,1)
        ss2 = 'valve has been closed'
    if quiet==0:
        print(ss1+' '+ss2) 


## Isolation valve - incident path
def iv_inc(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for hard open, 'so' or 'soft' for soft open, 'c' or 'close' for closing'''
    pv_r_soft = 'XF:11BMB-VA{Mir:KB-IV:1_Soft}Pos-Sts'
    pv_op_soft = 'XF:11BMB-VA{Mir:KB-IV:1_Soft}Cmd:Opn-Cmd'
    pv_cl_soft = 'XF:11BMB-VA{Mir:KB-IV:1_Soft}Cmd:Cls-Cmd'
    pv_r = 'XF:11BMB-VA{Mir:KB-IV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{Mir:KB-IV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{Mir:KB-IV:1}Cmd:Cls-Cmd'
    ss1 = 'Isolation valve for incident path: '
    dual_valve(cmd,pv_r_soft,pv_op_soft,pv_cl_soft,pv_r,pv_op,pv_cl,ss1,q)

## Isolation valve - sample/detector chamber
def iv_chm(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for hard open, 'so' or 'soft' for soft open, 'c' or 'close' for closing'''
    pv_r_soft = 'XF:11BMB-VA{Chm:Det-IV:1_Soft}Pos-Sts'
    pv_op_soft = 'XF:11BMB-VA{Chm:Det-IV:1_Soft}Cmd:Opn-Cmd'
    pv_cl_soft = 'XF:11BMB-VA{Chm:Det-IV:1_Soft}Cmd:Cls-Cmd'
    pv_r = 'XF:11BMB-VA{Chm:Det-IV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{Chm:Det-IV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{Chm:Det-IV:1}Cmd:Cls-Cmd'
    ss1 = 'Isolation valve for sample/WAXS chamber: '
    dual_valve(cmd,pv_r_soft,pv_op_soft,pv_cl_soft,pv_r,pv_op,pv_cl,ss1,q)

## Isolation valve - SAXS flightpath
def iv_pipe(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for hard open, 'so' or 'soft' for soft open, 'c' or 'close' for closing'''
    pv_r_soft = 'XF:11BMB-VA{BT:SAXS-IV:1_Soft}Pos-Sts'
    pv_op_soft = 'XF:11BMB-VA{BT:SAXS-IV:1_Soft}Cmd:Opn-Cmd'
    pv_cl_soft = 'XF:11BMB-VA{BT:SAXS-IV:1_Soft}Cmd:Cls-Cmd'
    pv_r = 'XF:11BMB-VA{BT:SAXS-IV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{BT:SAXS-IV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{BT:SAXS-IV:1}Cmd:Cls-Cmd'
    ss1 = 'Isolation valve for SAXS pipe: '
    dual_valve(cmd,pv_r_soft,pv_op_soft,pv_cl_soft,pv_r,pv_op,pv_cl,ss1,q)


## Vent valve - chamber upstream
def vv_us(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for hard open, 'so' or 'soft' for soft open, 'c' or 'close' for closing'''
    pv_r_soft = 'XF:11BMB-VA{Chm:Smpl-VV:1_Soft}Pos-Sts'
    pv_op_soft = 'XF:11BMB-VA{Chm:Smpl-VV:1_Soft}Cmd:Opn-Cmd'
    pv_cl_soft = 'XF:11BMB-VA{Chm:Smpl-VV:1_Soft}Cmd:Cls-Cmd'
    pv_r = 'XF:11BMB-VA{Chm:Smpl-VV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{Chm:Smpl-VV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{Chm:Smpl-VV:1}Cmd:Cls-Cmd'
    ss1 = 'Upstream vent valve for sample/WAXS chamber: '
    dual_valve(cmd,pv_r_soft,pv_op_soft,pv_cl_soft,pv_r,pv_op,pv_cl,ss1,q)

## Vent valve - chamber downstream
def vv_ds(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for hard open, 'so' or 'soft' for soft open, 'c' or 'close' for closing'''
    pv_r_soft = 'XF:11BMB-VA{Chm:Det-VV:1_Soft}Pos-Sts'
    pv_op_soft = 'XF:11BMB-VA{Chm:Det-VV:1_Soft}Cmd:Opn-Cmd'
    pv_cl_soft = 'XF:11BMB-VA{Chm:Det-VV:1_Soft}Cmd:Cls-Cmd'
    pv_r = 'XF:11BMB-VA{Chm:Det-VV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{Chm:Det-VV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{Chm:Det-VV:1}Cmd:Cls-Cmd'
    ss1 = 'Downstream vent valve for sample/WAXS chamber: '
    dual_valve(cmd,pv_r_soft,pv_op_soft,pv_cl_soft,pv_r,pv_op,pv_cl,ss1,q)


## Gate valve (Endstation) - upstream/small
def gv_us(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for open, 'c' or 'close' for closing'''
    pv_r = 'XF:11BMB-VA{Slt:4-GV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{Slt:4-GV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{Slt:4-GV:1}Cmd:Cls-Cmd'
    ss1 = 'Upstream/small gate valve (Endstation): '
    single_valve(cmd,pv_r,pv_op,pv_cl,ss1,q)

## Gate valve (Endstation) - downstream/large
def gv_ds(cmd='st',q=0):
    '''cmd: 'st' for status, 'o' or 'open' for open, 'c' or 'close' for closing'''
    pv_r = 'XF:11BMB-VA{Chm:Det-GV:1}Pos-Sts'
    pv_op = 'XF:11BMB-VA{Chm:Det-GV:1}Cmd:Opn-Cmd'
    pv_cl = 'XF:11BMB-VA{Chm:Det-GV:1}Cmd:Cls-Cmd'
    ss1 = 'Downstream/large gate valve (Endstation): '
    single_valve(cmd,pv_r,pv_op,pv_cl,ss1,q)


##### Endstation pumps
## Pump for flightpaths
def pump_fp(onoff, q=0):
    pv_r = 'XF:11BMB-VA{BT:SAXS-Pmp:1}Sts:Enbl-Sts'
    pv_w = 'XF:11BMB-VA{BT:SAXS-Pmp:1}Cmd:Enbl-Cmd'
    if onoff == 1:
        caput(pv_w,0)
        sleep(0.2)
        caput(pv_w,1)	
        ss='Flightpath pump has been turned ON'
    elif onoff == 0:
        caput(pv_w,0)	
        ss='Flightpath pump has been turned OFF'
    else:
        if caget(pv_r)==1:
            ss='Flightpath pump is ON'
        else:
            ss='Flightpath pump is OFF'
    if q==0:
        print(ss)

## Pump for sample/WAXS chamber
def pump_chm(onoff, q=0):
    pv_r = 'XF:11BMB-VA{Chm:Det-Pmp:1}Sts:Enbl-Sts'
    pv_w = 'XF:11BMB-VA{Chm:Det-Pmp:1}Cmd:Enbl-Cmd'
    if onoff == 1:
        caput(pv_w,0)
        sleep(0.2)
        caput(pv_w,1)	
        ss='Chamber pump has been turned ON'
    elif onoff == 0:
        caput(pv_w,0)	
        ss='Chamber pump has been turned OFF'
    else:
        if caget(pv_r)==1:
            ss='Chamber pump is ON'
        else:
            ss='Chamber pump is OFF'
    if q==0:
        print(ss)





