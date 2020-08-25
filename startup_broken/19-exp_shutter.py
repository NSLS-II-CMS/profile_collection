

##### Experimental shutters #####
# These shutters are controlled by sending a 5V pulse via QEM output on the Delta Tau controller MC06 
# (the same unit that controls slits S5). Both the opening and closing of the shutter are triggered 
# by the rise of the pulse. 
#
# Note: 
# - PV for the QEM output on MC06 is: 
#	XF:11BMB-CT{MC:06}Asyn.AOUT
# - This PV is located under Slit 5/Asyn --> asynRecord/More... --> asynOctet interface I/O --> ASCII 
# - 'M112=1' sets the state to high
# - 'M112=0' sets the state to low
# - 'M111=1' launches the change in state
# - A sleep time of ~2 ms between successive caput commands is needed to get proper response; 1 ms is too short
#####

#global xshutter_state
xshutter_state=0		## TODO: read the shutter state and set this accordingly

## Open shutter
def xshutter_trigger():
    sleep_time = 0.005 
    caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M112=1')
    sleep(sleep_time)
    caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M111=1')
    sleep(sleep_time)
    caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M112=0')
    sleep(sleep_time)
    caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M111=1')


def xshutter(inout,q=0):
    global xshutter_state

    if inout=='o' or inout=='open' or inout==1:
        if xshutter_state==0:
            xshutter_trigger()   
            xshutter_state = 1
            if q==0:
                print('Experimental shutter opened')
                return(xshutter_state)
        elif xshutter_state==1:
            print('Experimental shutter is already open; no changes made')
        else:
            print('xshutter_state is neither 0 nor 1; no changes made')

    if inout=='c' or inout=='close' or inout==0:
        if xshutter_state==1:
            xshutter_trigger()   
            xshutter_state = 0
            if q==0:
                print('Experimental shutter closed')
                return(xshutter_state)
        elif xshutter_state==0:
            print('Experimental shutter is already closed; no changes made')
        else:
            print('xshutter_state is neither 0 nor 1; no changes made')


