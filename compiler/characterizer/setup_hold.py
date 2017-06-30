import sys
import globals
import tech
import stimuli
import debug
import charutils as ch
import ms_flop

OPTS = globals.get_opts()


class setup_hold():
    """
    Functions to calculate the setup and hold times of the SRAM
    (Bisection Methodology)
    """

    def __init__(self):
        # This must match the spice model order
        self.pins = ["data", "dout", "dout_bar", "clk", "vdd", "gnd"]
        self.model_name = "ms_flop"
        self.model_location = OPTS.openram_tech + "sp_lib/ms_flop.sp"
        self.period = tech.spice["feasible_period"]
        self.vdd = tech.spice["supply_voltage"]
        self.gnd = tech.spice["gnd_voltage"]

        debug.info(2,"Feasible period from technology file: {0} ".format(self.period))

                


    def write_stimulus(self, mode, target_time, correct_value):
        """Creates a stimulus file for SRAM setup/hold time calculation"""

        # creates and opens the stimulus file for writing
        temp_stim = OPTS.openram_temp + "stim.sp"
        self.sf = open(temp_stim, "w")

        self.write_header(correct_value)

        # instantiate the master-slave d-flip-flop
        self.sf.write("* Instantiation of the Master-Slave D-flip-flop\n")
        stimuli.inst_model(stim_file=self.sf,
                           pins=self.pins,
                           model_name=self.model_name)

        # create a buffer for the inputs
        self.sf.write("* Buffer subckt\n")
        stimuli.create_buffer(stim_file=self.sf,
                              buffer_name="buffer",
                              size=[1, 1])

        self.write_data(mode=mode,
                        target_time=target_time,
                        correct_value=correct_value)

        self.write_clock()

        self.write_measures(mode=mode, 
                            correct_value=correct_value)
                         

        stimuli.write_control(self.sf,2.5*self.period)

        self.sf.close()

    def write_header(self, correct_value):
        """ Write the header file with all the models and the power supplies. """
        self.sf.write("* Stimulus for setup/hold: data {0} period {1}n\n".format(correct_value, self.period))

        # include files in stimulus file
        self.model_list = tech.spice["fet_models"] + [self.model_location]
        stimuli.write_include(stim_file=self.sf,
                              models=self.model_list)

        # add vdd/gnd statements
        self.sf.write("* Global Power Supplies\n")
        stimuli.write_supply(self.sf)


    def write_data(self, mode, target_time, correct_value):
        """ Create the buffered data signals for setup/hold analysis """
        self.sf.write("* Generation of the data and clk signals\n")
        incorrect_value = stimuli.get_inverse_value(correct_value)
        if mode=="HOLD":
            start_value = correct_value
            end_value = incorrect_value
        else:
            start_value = incorrect_value
            end_value = correct_value

        stimuli.gen_pulse(stim_file=self.sf,
                          sig_name="data",
                          v1=start_value,
                          v2=end_value,
                          offset=target_time,
                          period=2*self.period,
                          t_rise=self.constrained_input_slew,
                          t_fall=self.constrained_input_slew)

    def write_clock(self):
        """ Create the buffered clock signal for setup/hold analysis """
        stimuli.gen_pulse(stim_file=self.sf,
                          sig_name="clk",
                          offset=self.period,
                          period=self.period,
                          t_rise=self.related_input_slew,
                          t_fall=self.related_input_slew)



    def write_measures(self, mode, correct_value):
        """ Measure statements for setup/hold with right phases. """

        if correct_value == self.vdd:
            dout_rise_or_fall = "RISE"
        else:
            dout_rise_or_fall = "FALL"

        if mode == "SETUP":
            din_rise_or_fall = dout_rise_or_fall
        else:
            if correct_value == self.vdd:
                din_rise_or_fall = "FALL"
            else:
                din_rise_or_fall = "RISE"

            
        incorrect_value = stimuli.get_inverse_value(correct_value)

        self.sf.write("* Measure statements for pass/fail verification\n")
        self.sf.write(".IC v({0})={1}\n".format("dout", incorrect_value))

        trig_name = "clk"
        targ_name = "dout"
        trig_val = targ_val = 0.5 * self.vdd
        stimuli.gen_meas_delay(stim_file=self.sf,
                               meas_name="clk2q_delay",
                               trig_name=trig_name,
                               targ_name=targ_name,
                               trig_val=trig_val,
                               targ_val=targ_val,
                               trig_dir="RISE",
                               targ_dir=dout_rise_or_fall,
                               td=0)

        targ_name = "data"
        stimuli.gen_meas_delay(stim_file=self.sf,
                               meas_name="setup_hold_time",
                               trig_name=trig_name,
                               targ_name=targ_name,
                               trig_val=trig_val,
                               targ_val=targ_val,
                               trig_dir="RISE",
                               targ_dir=din_rise_or_fall,
                               td=0)
        



    def bidir_search(self, correct_value, mode):
        """ This will perform a bidirectional search for either setup or hold times.
        It starts with the feasible priod and looks a half period beyond or before it
        depending on whether we are doing setup or hold. 
        """

        # NOTE: The feasible bound is always feasible. This is why they are different for setup and hold.
        # The clock will always be offset by a period fromt he start, so we want to look before and after
        # this time by half a period. However, the setup/hold is measured according to the BUFFERED clock, not the ideal clock.
        
        if mode == "SETUP":
            feasible_bound = 0.5*self.period
            infeasible_bound = 2*self.period
        else:
            infeasible_bound = 0.5*self.period
            feasible_bound = 2*self.period

        # Initial check if reference feasible bound time passes for correct_value, if not, we can't start the search!
        self.write_stimulus(mode=mode, 
                            target_time=feasible_bound, 
                            correct_value=correct_value)
        stimuli.run_sim()
        ideal_clk_to_q = ch.convert_to_float(ch.parse_output("timing", "clk2q_delay"))
        setuphold_time = ch.convert_to_float(ch.parse_output("timing", "setup_hold_time"))
        debug.info(2,"*** {0} CHECK: {1} Ideal Clk-to-Q: {2} Setup/Hold: {3}".format(mode, correct_value,ideal_clk_to_q,setuphold_time))

        if type(ideal_clk_to_q)!=float or type(setuphold_time)!=float:
            debug.error("Initial hold time fails for data value feasible bound {0} Clk-to-Q {1} Setup/Hold {2}".format(feasible_bound,ideal_clk_to_q,setuphold_time),2)

        if mode == "SETUP": # SETUP is clk-din, not din-clk
            setuphold_time *= -1e9
        else:
            setuphold_time *= 1e9
            
        debug.info(2,"Checked initial {0} time {1}, data at {2}, clock at {3} ".format(mode,
                                                                                       setuphold_time,
                                                                                       feasible_bound,
                                                                                       self.period))

        while True:
            target_time = (feasible_bound + infeasible_bound)/2
            self.write_stimulus(mode=mode, 
                                target_time=target_time, 
                                correct_value=correct_value)

            debug.info(2,"Correct value: {0} Target time: {1} Infeasible: {2} Feasible: {3}".format(correct_value,
                                                                                       target_time,
                                                                                       infeasible_bound,
                                                                                       feasible_bound))


            stimuli.run_sim()
            clk_to_q = ch.convert_to_float(ch.parse_output("timing", "clk2q_delay"))
            setuphold_time = ch.convert_to_float(ch.parse_output("timing", "setup_hold_time"))
            if type(clk_to_q)==float and (clk_to_q<1.1*ideal_clk_to_q) and type(setuphold_time)==float:
                if mode == "SETUP": # SETUP is clk-din, not din-clk
                    setuphold_time *= -1e9
                else:
                    setuphold_time *= 1e9

                debug.info(2,"PASS Clk-to-Q: {0} Setup/Hold: {1}".format(clk_to_q,setuphold_time))
                passing_setuphold_time = setuphold_time
                feasible_bound = target_time
            else:
                debug.info(2,"FAIL Clk-to-Q: {0} Setup/Hold: {1}".format(clk_to_q,setuphold_time))
                infeasible_bound = target_time

            #raw_input("Press Enter to continue...")
            if ch.relative_compare(feasible_bound, infeasible_bound, error_tolerance=0.001):
                debug.info(3,"CONVERGE {0} vs {1}".format(feasible_bound,infeasible_bound))
                break
            

        debug.info(2,"Converged on {0} time {1}.".format(mode,passing_setuphold_time))
        return passing_setuphold_time


    def setup_LH_time(self):
        """Calculates the setup time for low-to-high transition for a DFF
        """
        return self.bidir_search(self.vdd, "SETUP")


    def setup_HL_time(self):
        """Calculates the setup time for high-to-low transition for a DFF
        """
        return self.bidir_search(self.gnd, "SETUP")
    
    def hold_LH_time(self):
        """Calculates the hold time for low-to-high transition for a DFF
        """
        return self.bidir_search(self.vdd, "HOLD")

    def hold_HL_time(self):
        """Calculates the hold time for high-to-low transition for a DFF
        """
        return self.bidir_search(self.gnd, "HOLD")


    def analyze(self,related_slews, constrained_slews):
        """main function to calculate both setup and hold time for the
        d-flip-flop returns a dictionary that contains 4 times for both
        setup/hold times for high_to_low and low_to_high transition
        dictionary keys: setup_time_one (low_to_high), setup_time_zero
        (high_to_low), hold_time_one (low_to_high), hold_time_zero
        (high_to_low)
        """
        LH_setup = []
        HL_setup = []
        LH_hold = []
        HL_hold = []
        for self.related_input_slew in related_slews:
            for self.constrained_input_slew in constrained_slews:
                debug.info(1, "Clock slew: {0} Data slew: {1}".format(self.related_input_slew,self.constrained_input_slew))
                LH_setup_time = self.setup_LH_time()
                debug.info(1, "  Setup Time for low_to_high transistion: {0}".format(LH_setup_time))
                HL_setup_time = self.setup_HL_time()
                debug.info(1, "  Setup Time for high_to_low transistion: {0}".format(HL_setup_time))
                LH_hold_time = self.hold_LH_time()
                debug.info(1, "  Hold Time for low_to_high transistion: {0}".format(LH_hold_time))
                HL_hold_time = self.hold_HL_time()
                debug.info(1, "  Hold Time for high_to_low transistion: {0}".format(HL_hold_time))
                LH_setup.append(LH_setup_time)
                HL_setup.append(HL_setup_time)
                LH_hold.append(LH_hold_time)
                HL_hold.append(HL_hold_time)
                
        times = {"setup_times_LH": LH_setup,
                 "setup_times_HL": HL_setup,
                 "hold_times_LH": LH_hold,
                 "hold_times_HL": HL_hold
                 }
        return times

    def analytical_model(self,related_slews, constrained_slews):
        """ Just return the fixed setup/hold times from the technology.
        """
        LH_setup = []
        HL_setup = []
        LH_hold = []
        HL_hold = []
        
        for self.related_input_slew in related_slews:
            for self.constrained_input_slew in constrained_slews:
                # convert from ps to ns
                LH_setup.append(tech.spice["msflop_setup"]/1e3)
                HL_setup.append(tech.spice["msflop_setup"]/1e3)
                LH_hold.append(tech.spice["msflop_hold"]/1e3)
                HL_hold.append(tech.spice["msflop_hold"]/1e3)
                
        times = {"setup_times_LH": LH_setup,
                 "setup_times_HL": HL_setup,
                 "hold_times_LH": LH_hold,
                 "hold_times_HL": HL_hold
                 }
        return times
