#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4




################################################################################
#  Classes for controlling the robotics and automation on the beamline.
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
#
# Config file or globalstore to save information.
################################################################################





class SampleExchangeRobot(Stage):
    
    def __init__(self, name='SampleExchangeRobot', base=None, use_gs=True, **kwargs):
        
        super().__init__(name=name, base=base, **kwargs)
        
        self._sample = None
        
        # The region can be:
        #  'safe' : arm won't collid with anything, it is near the (+,+,+) limit of its travel.
        #  'parking' : arm is close to the parking lot (movement may hit a sample)
        #  'stage' : arm is close to the sample stage/stack (movement may collide with stack, on-axis camera, or downstream window)
        #  'undefined' : position is unknown (do not assume it is safe to move!)
        self._region = 'undefined'
        
   
       
        # self.yabs(-82.0) # Good height for 'slotted approach'
        # self.yabs(-77.0) # Good height for 'grip' (grip-screws sitting at bottom of wells)
        # self.yabs(-67.0) # Good height for 'hover' (sample held above stage)
        self._delta_y_hover = 5.0
        self._delta_y_slot = 4.0

        #'SAFE' position of gripper
        self._position_safe = [0, -104.9, 0.0, 0.0, +90] # x, y, z, r, phi
        
        #self._position_sample_gripped = [-100, -104.9, -94.8, 18.0, +90] # x, y, z, r, phi
        #self._position_hold = [0, -104.9, -94.8, 0.0, +90] # x, y, z, r, phi
        
        #default position of gripper to pick up from Stage
        #self._position_sample_gripped = [ -100.49986249999999, -102.89986875, -93.7, 17.5, 90.0 ] # x, y, z, r, phi
        #self._position_hold = [ 0, -106.89986875, -94.04984999999999, 0, 90.0 ] # x, y, z, r, phi
        
        #tested without SmarAct stage. smx=50; smy=-2.37
        #self._position_sample_gripped = [ -99, -107, -94, 0.0, 91 ] # x, y, z, r, phi
        #self._position_hold = [ 0, -107, -94, 0, 91 ] # x, y, z, r, phi

        #tested with the gripper with spring. smx=50; smy=-2.37
        self._position_sample_gripped = [ -98, -103, -94.5, 0.0, 91 ] # x, y, z, r, phi
        self._position_hold = [ 0, -103, -94.5, 0, 91 ] # x, y, z, r, phi


        #defacult position of gripper to pick up from Garage(1,1)
        #self._position_garage = [-96, -200, -129.5, 0.0, 0.0] # x, y, z, r, phi
        #self._position_garage = [ -97.5, -201.000121875, -129.9003625, -0.412427, 0.0 ] # x, y, z, r, phi
        #self._position_garage = [ -98.999675, -198, -130.000375, -0.621273, 0.0 ] # x, y, z, r, phi
        #self._position_garage = [ -98, -200, -128, 0.0, 1 ] # x, y, z, r, phi

        #tested with the gripper with spring. smx=50; smy=-2.37
        self._position_garage = [ -96, -197.5, -127, 0.0, 1 ] # x, y, z, r, phi

        #default position for stage
        #self._position_stg_exchange = [+30.0, -2.37] # smx, smy
        #self._position_stg_safe = [-30.0, -2.37] # smx, smy
        #self._position_stg_measure = [] # smx, smy

        #default position for stage without SmarAct motor
        self._position_stg_exchange = [+50.0, -2.37] # smx, smy
        self._position_stg_safe = [-30.0, -2.37] # smx, smy
        self._position_stg_measure = [] # smx, smy
        
        #garage structure parameter 
        self._delta_garage_x = 44.45 # 1.75 inch = 44.45 mm
        self._delta_garage_y = 63.5 # 2.5 inch = 63.5 mm
        
        if use_gs and hasattr(gs, 'robot'):
            if '_position_sample_gripped' in gs.robot:
                self._position_sample_gripped = gs.robot['_position_sample_gripped']
            if '_position_hold' in gs.robot:
                self._position_hold = gs.robot['_position_hold']
            if '_position_garage' in gs.robot:
                self._position_garage = gs.robot['_position_garage']
        else:
            gs.robot = {}
        
        for axis_name, axis in self._axes.items():
            axis._move_settle_max_time = 30.0
        

        
    def _set_axes_definitions(self):
        '''Internal function which defines the axes for this stage. This is kept
        as a separate function so that it can be over-ridden easily.'''
        
        # The _axes_definitions array holds a list of dicts, each defining an axis
        self._axes_definitions = [ {'name': 'x',
                            'motor': armx,
                            'enabled': True,
                            'scaling': +1.0,
                            'units': 'mm',
                            'hint': 'positive moves left/outboard',
                            },
                            {'name': 'r',
                            'motor': armr,
                            #'motor': strans,
                            'enabled': True,
                            'scaling': +1.0,
                            'units': 'mm',
                            'hint': 'positive moves radial arm outwards',
                            },        
                            {'name': 'y',
                            'motor': army,
                            'enabled': True,
                            'scaling': +1.0,
                            'units': 'mm',
                            'hint': 'positive moves arm up',
                            },
                            {'name': 'z',
                            'motor': armz,
                            'enabled': True,
                            'scaling': +1.0,
                            'units': 'mm',
                            'hint': 'positive moves arm downstream',
                            },
                            {'name': 'phi',
                            'motor': armphi,
                            'enabled': True,
                            'scaling': +1.0,
                            'units': 'mm',
                            'hint': 'positive moves arm downstream',
                            },
                            
                        ]      
                            


                                
                            
    def home_y(self, verbosity=3, delays=0.5, retries=5, max_wait=2.0):
        
        if retries<1:
            print('ERROR: home_y failed (too many retries).')
            return False
        
        # Activate homing
        caput('XF:11BMB-ES{SM:1-Ax:Y}Start:Home-Cmd', 1)
        
        # Make sure homing actually started
        start_time = time.time()
        while caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Homing-Sts')!=1 and time.time()-start_time<max_wait:
            sleep(.01)
            if verbosity>=5:
                print( 'phase, status, homing = {}, {}, {}'.format( caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Home-Sts'), caget('XF:11BMB-ES{SM:1-Ax:Y}Pgm:Home-Sts'), caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Homing-Sts') ) )
            
        
        if time.time()-start_time>max_wait:
            # Retry
            self.home_y(verbosity=verbosity, delays=delays, retries=retries-1, max_wait=max_wait)
            
            
        # Wait for motion to be complete
        sleep(delays)
        while army.moving or caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Homing-Sts')!=0:
            sleep(delays)
            
        
        if abs(self.yabs(verbosity=0)-0)>0.1:
            print("ERROR: y didn't home (position = {})".format(self.yabs(verbosity=0)))
            return False
        
        if caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Home-Sts')!=7 or caget('XF:11BMB-ES{SM:1-Ax:Y}Pgm:Home-Sts')!=0 or caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Homing-Sts')!=0:
            print( 'phase, status, homing = {}, {}, {}'.format( caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Home-Sts'), caget('XF:11BMB-ES{SM:1-Ax:Y}Pgm:Home-Sts'), caget('XF:11BMB-ES{SM:1-Ax:Y}Sts:Homing-Sts') ) )
            print('ERROR: y homing failed.')
            return False
        
        return True
            
        
            
            
        
    def home(self, verbosity=3, delays=0.5):
        '''Home the axes, so that one can now trust the position values.'''
        
        if not self.checkSafe():
            return
        
        if self._sample is not None:
            print("ERROR: You shouldn't home the robot with a sample in the gripper.")
            return
        
        if verbosity>=3:
            print('Homing robot')
        

        # army to positive limit (moves arm to top of vertical range); set this to be zero
        success = self.home_y(verbosity=verbosity, delays=delays)
        if not success:
            return False
        
        # armx to positive limit; set this to be zero
        self.xabs(0)
        caput('XF:11BMB-ES{SM:1-Ax:X}Start:Home-Cmd', 1)
        sleep(delays)
        while armx.moving:
            sleep(delays)
        
        # Rotate arm so that it doesn't collide when doing a +z scan
        self.phiabs(+90) # gripper pointing -x (towards sample stack)
        sleep(delays)
        while armphi.moving:
            sleep(delays)
            
        # Home armphi 
        caput('XF:11BMB-ES{SM:1-Ax:Yaw}Start:Home-Cmd', 1)
        sleep(delays)
        while armphi.moving and caget('XF:11BMB-ES{SM:1-Ax:Yaw}Sts:Homing-Sts')!=0:
            sleep(delays)
        
        # Rotate arm so that it doesn't collide when doing a +z scan
        self.phiabs(+90) # gripper pointing -x (towards sample stack)
        sleep(delays)
        while armphi.moving:
            sleep(delays)
        
        
        # armz to positive limit (moves arm to downstream of range); set this to be zero
        caput('XF:11BMB-ES{SM:1-Ax:Z}Start:Home-Cmd', 1)
        sleep(delays)
        while armz.moving:
            sleep(delays)
        
        
        #caput('XF:11BMB-ES{SM:1-Ax:ArmR}Mtr.HOMF',1) # armr home forward
        caput('XF:11BMB-ES{SM:1-Ax:ArmR}Mtr.HOMR',1) # armr home reverse
        sleep(delays)
        while self._axes['r'].motor.moving:
            sleep(delays)
        
        self._region = 'safe'
        return True
        
        
    def checkSafe(self, check_stage=True):

        if self._region is not 'safe':
            print("ERROR: Robot arm must start in the 'safe' region of the chamber (current region is '{}'). Move the robot to the safe region (and/or set _region to 'safe').".format(self._region))
            return False
        
        #smx_safe, smy_safe = self._position_stg_safe
        #if check_stage and smx.user_readback.value > (smx_safe+0.1):
            #print("ERROR: smx ({}) is in an unsafe position.".format(smx.user_readback.value))
            #return False
        
        return True
        
        
        
    def setReferenceSampleGripped(self):
        """The position when the sample (on the stage) is gripped."""
        
        x = self.xpos(verbosity=0)
        y = self.ypos(verbosity=0)
        z = self.zpos(verbosity=0)
        r = self.rpos(verbosity=0)
        phi = self.phipos(verbosity=0)
        self._position_sample_gripped = x, y, z, r, phi # x, y, z, r, phi

        print("self._position_sample_gripped = [ {}, {}, {}, {}, {} ] # x, y, z, r, phi".format(x, y, z, r, phi))

        if hasattr(gs, 'robot'):
            gs.robot['_position_sample_gripped'] = self._position_sample_gripped

        
        self._position_hold = 0, y, z, 0, phi # x, y, z, r, phi
        
        print("self._position_hold = [ {}, {}, {}, {}, {} ] # x, y, z, r, phi".format(0, y, z, 0, phi))

        if hasattr(gs, 'robot'):
            gs.robot['_position_hold'] = self._position_hold

        
    def setReferenceGarage(self):
        """The position when the lower-left (1,1) sample of the garage is gripped."""
        
        x = self.xpos(verbosity=0)
        y = self.ypos(verbosity=0)
        z = self.zpos(verbosity=0)
        r = self.rpos(verbosity=0)
        phi = self.phipos(verbosity=0)
        
        self._position_garage = x, y, z, r, phi # x, y, z, r, phi
        
        print("self._position_garage = [ {}, {}, {}, {}, {} ] # x, y, z, r, phi".format(x, y, z, r, phi))

        if hasattr(gs, 'robot'):
            gs.robot['_position_garage'] = self._position_garage
        
        
    def motionSlot(self, direction):
        self.yr(direction*self._delta_y_slot)

    def motionHover(self, direction):
        self.yr(direction*self._delta_y_hover)
        
    def sequenceGotoSafe(self, verbosity=3):
        
        x, y, z, r, phi = self._position_safe
        
        #if abs( self.phipos(verbosity=verbosity) - 90 ) < 0.1:
            ## phi = 90deg, prongs pointed at stage
            #pass
        #elif abs( self.phipos(verbosity=verbosity) - 0 ) < 0.1:
            ## phi = 0deg, prongs pointed at stage
            #pass
        #else:
            #pass
            
        self.phiabs(phi)
        self.yabs(y)
        self.xabs(x)
        self.zabs(z)
            
        self._region = 'safe'


    def sequenceGotoSampleStageSlotted(self, x_motion=True, verbosity=3):
        
        if self._sample is not None:
            print("ERROR: There is already a sample being gripped by robot arm (sample {}.".format(self._sample.name))
            return
        
        if not self.checkSafe():
            return
        
        if verbosity>=3:
            print('TBD')
        
        
        # Move sample stage
        x, y = self._position_stg_exchange # smx, smy
        mov([smx, smy], [x,y])
        
        
        x, y, z, r, phi = self._position_sample_gripped
        
        # Pre-align the arm in (y,z)
        self.phiabs(phi, verbosity=verbosity)
        self.zabs(z, verbosity=verbosity)
        self.yabs(y-self._delta_y_slot, verbosity=verbosity)
        
        self._region = 'stage'
        if x_motion:
            # Move arm so that it is slotted over the sample
            self.xabs(x, verbosity=verbosity)
            #self.rabs(r, verbosity=verbosity)
        
        
    def sequencePutSampleOntoStage(self, verbosity=3):
        
        if self._sample is None:
            print("ERROR: No sample currently being gripped by robot arm.")
            return
        
        if not self.checkSafe(check_stage=False):
            return
        
        if verbosity>=3:
            print('Putting sample onto stage')
        
        
        # Move sample stage
        x, y = self._position_stg_exchange # smx, smy
        mov([smx, smy], [x,y])

        
        x, y, z, r, phi = self._position_sample_gripped
        
        # Pre-align the arm in (y,z)
        self.phiabs(phi, verbosity=verbosity)
        self.zabs(z, verbosity=verbosity)
        self.yabs(y+self._delta_y_hover, verbosity=verbosity)
        
        self._region = 'stage'
        # Push the sample out so that it is hovering above the stage
        # r is removed without SmarAct motor
        #mov([armx, self._axes['r'].motor], [x, r])
        mov(armx,x)
        
        # Move sample down (-y)
        self.yr(-self._delta_y_hover, verbosity=verbosity) # Now in contact with stage
        
        # De-grip
        self.yr(-self._delta_y_slot, verbosity=verbosity)
        self._sample = None
        
        # Move away from stage
        x, y, z, r, phi = self._position_hold
        #mov([armx, self._axes['r'].motor], [x, r])
        # r is removed without SmarAct motor
        mov(armx, x)
        
        self.sequenceGotoSafe(verbosity=verbosity)
        
        

    def sequenceGetSampleFromStage(self, verbosity=3):
        
        if self._sample is not None:
            print("ERROR: There is already a sample being gripped by robot arm (sample {}.".format(self._sample.name))
            return
        
        if not self.checkSafe():
            return
        
        if verbosity>=3:
            print('Getting sample from stage')
        
        
        # Move sample stage
        x, y = self._position_stg_exchange # smx, smy
        mov([smx, smy], [x,y])
        
        
        x, y, z, r, phi = self._position_sample_gripped
        
        # Pre-align the arm in (y,z)
        self.phiabs(phi, verbosity=verbosity)
        self.zabs(z, verbosity=verbosity)
        self.yabs(y-self._delta_y_slot, verbosity=verbosity)
        
        self._region = 'stage'
        # Move arm so that it is slotted over the sample
        self.xabs(x, verbosity=verbosity)
        # self.rabs(r, verbosity=verbosity)
        
        # Grip sample
        self.yr(+self._delta_y_slot, verbosity=verbosity)
        
        self._sample = 'exists' # TODO: Fix
        
        # Pick sample up (+y)
        self.yr(+self._delta_y_hover, verbosity=verbosity)
        
        # Move away from stage
        x, y, z, r, phi = self._position_hold
        #mov([armx, self._axes['r'].motor], [x, r])
        # r is removed without SmarAct motor
        mov(armx, x)
        
        self.sequenceGotoSafe(verbosity=verbosity)
        
        
        
    def sequenceGetSampleFromGarage(self, shelf_num, spot_num, verbosity=3):
        
        if shelf_num<1 or shelf_num>4:
            print("ERROR: Invalid shelf {}".format(shelf_num))
            return
        if spot_num<1 or spot_num>3:
            print("ERROR: Invalid spot {}".format(spot_num))
            return
        
        if self._sample is not None:
            print("ERROR: There is already a sample being gripped by robot arm (sample {}.".format(self._sample.name))
            return
        
        if verbosity>=3:
            print('Getting sample from garage ({}, {})'.format(shelf_num, spot_num))
        
        
        x, y, z, r, phi = self._position_garage
        
        self.phiabs(phi)
        
        # Align ourselves with this parking spot
        success = self.sequencePrepGarageXY(shelf_num, spot_num, verbosity=verbosity)
        if not success:
            return
        x = self.xpos(verbosity=0)
        y = self.ypos(verbosity=0)
        
        # Lower so that the slot is aligned
        self.yr(-self._delta_y_slot, verbosity=verbosity)
        
        # Move towards parking lot
        self._region = 'parking'
        self.zabs(z, verbosity=verbosity)
        
        # Grip sample
        self.yr(+self._delta_y_slot, verbosity=verbosity)
        
        self._sample = 'exists' # TODO: Fix

        # Pick sample up (+y)
        self.yr(+self._delta_y_hover, verbosity=verbosity)

        # Move away from parking
        self.zabs(0, verbosity=verbosity)
        self.xabs(0, verbosity=verbosity)
        self.yabs(self._position_safe[1], verbosity=verbosity)
        self.sequenceGotoSafe(verbosity=verbosity)
        

    def sequencePutSampleInGarage(self, shelf_num, spot_num, verbosity=3):
        
        if shelf_num<1 or shelf_num>4:
            print("ERROR: Invalid shelf {}".format(shelf_num))
            return
        if spot_num<1 or spot_num>3:
            print("ERROR: Invalid spot {}".format(spot_num))
            return
        
        if self._sample is None:
            print("WARNING: There is no sample being gripped by robot arm.")
        
        if verbosity>=3:
            print('Putting sample into garage ({}, {})'.format(shelf_num, spot_num))
        
        
        x, y, z, r, phi = self._position_garage
        
        self.phiabs(phi)
        
        # Align ourselves with this parking spot
        success = self.sequencePrepGarageXY(shelf_num, spot_num, verbosity=verbosity)
        if not success:
            return
        x = self.xpos(verbosity=0)
        y = self.ypos(verbosity=0)
        
        # Hover
        self.yr(+self._delta_y_hover, verbosity=verbosity)
        
        # Move towards parking lot
        self._region = 'parking'
        self.zabs(z, verbosity=verbosity)
        
        # Deposit sample
        self.yr(-self._delta_y_hover, verbosity=verbosity)
        self.yr(-self._delta_y_slot, verbosity=verbosity)
        
        self._sample = None

        # Move away from parking
        self.zabs(0, verbosity=verbosity)
        self.xabs(0, verbosity=verbosity)
        self.sequenceGotoSafe(verbosity=verbosity)


    def sequencePrepGarageXY(self, shelf_num, spot_num, verbosity=3):
        
        if shelf_num<1 or shelf_num>4:
            print("ERROR: Invalid shelf {}".format(shelf_num))
            return
        if spot_num<1 or spot_num>3:
            print("ERROR: Invalid spot {}".format(spot_num))
            return

        if not self.checkSafe():
            return
        
        if verbosity>=3:
            print('  Going to garage ({}, {})'.format(shelf_num, spot_num))
        
        if self.zpos(verbosity=0)<-10:
            print("ERROR: z ({}) position unsafe.".format(self.zpos(verbosity=0)))

        if abs(self.phipos(verbosity=0))>1:
            print("ERROR: phi ({}) position unsafe.".format(self.phipos(verbosity=0)))
        
        x, y, z, r, phi = self._position_garage
        
        x += (spot_num-1)*self._delta_garage_x
        y += (shelf_num-1)*self._delta_garage_y
        
        if verbosity>=4:
            print('    Going to (x,y) = ({}, {})'.format(x, y))
        
        
        # Do y first to avoid catching cables
        self.yabs(y)
        self.xabs(x)
        
        xactual = self.xpos(verbosity=0)
        yactual = self.ypos(verbosity=0)
        
        if abs(x-xactual)>0.2:
            print('ERROR: x did not reach requested position (request = {}, actual = {})'.format(x, xactual))
            return False
        
        if abs(y-yactual)>0.2:
            print('ERROR: y did not reach requested position (request = {}, actual = {})'.format(y, yactual))
            return False
        
        return True
        
        
    def loadSample(self, shelf_num, spot_num, verbosity=3):
        
        # Check if a sample is on stage
        # Unload sample if necessary
        
        self.sequenceGetSampleFromGarage(shelf_num, spot_num, verbosity=verbosity)
        self.sequencePutSampleOntoStage()
    
    
    def _testing_calibrationStage(self, verbosity=3):
        
        if self._sample is not None:
            print("ERROR: Calibration cannot be done with sample on robot arm.")
            return
        
        if not self.checkSafe(check_stage=False):
            return
        
        if verbosity>=3:
            print('Approaching to stage')
        
        # Move sample stage
        x, y = self._position_stg_exchange # smx, smy
        mov([smx, smy], [x,y])
        
        
        x, y, z, r, phi = self._position_sample_gripped
        
        # Pre-align the arm in (y,z)
        self.phiabs(phi, verbosity=verbosity)
        self.zabs(z, verbosity=verbosity)
        self.yabs(y-self._delta_y_slot, verbosity=verbosity)
        
        # Move to a position x is away from sample bar
        self.xabs(-40, verbosity=verbosity)
        self.rabs(0, verbosity=verbosity)
         
    def _testing_calibrationGarage(self, verbosity=3):
        
        #use Garage(1,1) to calibration the gripper position
        shelf_num=1
        spot_num=1
  
        if self._sample is not None:
            print("ERROR: Calibration cannot be done with sample on robot arm. ")
            return
        
        if verbosity>=3:
            print('Approaching to garage ({}, {})'.format(shelf_num, spot_num))
        
        
        x, y, z, r, phi = self._position_garage
        
        self.phiabs(phi)
        
        # Align ourselves with this parking spot
        success = self.sequencePrepGarageXY(shelf_num, spot_num, verbosity=verbosity)
        if not success:
            return
        x = self.xpos(verbosity=0)
        y = self.ypos(verbosity=0)
        
        # Lower so that the slot is aligned
        self.yr(-self._delta_y_slot, verbosity=verbosity)
        self.zabs(-60)
        
       

    def _stress_test(self, cycles=2, verbosity=5):

        if not self.checkSafe():
            return
        
        self.home()
        
        if not self.checkSafe():
            return
        
        
        for i in range(cycles):
            if verbosity>=2:
                print('Stress test cycle {}'.format(i))
                
            for shelf_num in range(1, 4+1):
                for spot_num in range(1, 3+1):

                    if verbosity>=2:
                        print('Stress test garage ({}, {})'.format(shelf_num, spot_num))

                    
                    self.sequenceGetSampleFromGarage(shelf_num, spot_num, verbosity=verbosity)
                    self.sequencePutSampleOntoStage(verbosity=verbosity)
                    
                    sleep(3)
                                        
                    self.sequenceGetSampleFromStage(verbosity=verbosity)
                    self.sequencePutSampleInGarage(shelf_num, spot_num, verbosity=verbosity)
                    

    def run(self, cycles=1, verbosity=5):

        if not self.checkSafe():
            return
        
        #self.home()
        
        
        for i in range(cycles):
            if verbosity>=2:
                print('Run test cycle {}'.format(i))

            for hol in Garage_holders:
                
                    [shelf_num, spot_num] = hol.GaragePosition

            
                    if verbosity>=2:
                        print('Run test garage ({}, {})'.format(shelf_num, spot_num))
    
                    
                    self.sequenceGetSampleFromGarage(shelf_num, spot_num, verbosity=verbosity)
                    self.sequencePutSampleOntoStage(verbosity=verbosity)

                    hol.listSamples()
                    sleep(3)
                    hol.doSamples()
                    
                    
                    self.sequenceGetSampleFromStage(verbosity=verbosity)
                    self.sequencePutSampleInGarage(shelf_num, spot_num, verbosity=verbosity)        
        
    def listGarage(self, verbosity=3):

        for hol in Garage_holders:
               
                [shelf_num, spot_num] = hol.GaragePosition
                print('In Garage ({}, {})'.format(shelf_num, spot_num))
                hol.listSamplesPositions()        
        
        

        

class Queue(object):
    """
    Holds the current state of the sample queue, allowing samples settings
    to be 'extracted'; or even allowing a particular sample to be physically
    loaded.
    """
    
    pass




# Note: This will break until class is updated to not use gs at all.
# robot = SampleExchangeRobot(use_gs=False)

#robot._region='safe'
