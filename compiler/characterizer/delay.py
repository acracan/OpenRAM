import sys
import re
import globals
import debug
import tech
import math
import stimuli
import charutils as ch
import utils

OPTS = globals.get_opts()

class delay():
    """
    Functions to measure the delay of the SRAM at a given address and
    data bit.
    """

    def __init__(self,sram,spfile):
        self.name = sram.name
        self.num_words = sram.num_words
        self.word_size = sram.word_size
        self.addr_size = sram.addr_size
        self.sram_sp_file = spfile
        
        self.vdd = tech.spice["supply_voltage"]
        self.gnd = tech.spice["gnd_voltage"]

    def check_arguments(self):
        """Checks if arguments given for write_stimulus() meets requirements"""
        try:
            int(self.probe_address, 2)
        except ValueError:
            debug.error("Probe Address is not of binary form: {0}".format(self.probe_address),1)

        if len(self.probe_address) != self.addr_size:
            debug.error("Probe Address's number of bits does not correspond to given SRAM",1)

        if not isinstance(self.probe_data, int) or self.probe_data>self.word_size or self.probe_data<0:
            debug.error("Given probe_data is not an integer to specify a data bit",1)


    def write_stimulus(self, period):
        """Creates a stimulus file for simulations to probe a certain bitcell, given an address and data-position of the data-word 
        (probe-address form: '111010000' LSB=0, MSB=1)
        (probe_data form: number corresponding to the bit position of data-bus, begins with position 0) 
        """
        self.check_arguments()

        # obtains list of time-points for each rising clk edge
        self.obtain_cycle_times(period)

        # creates and opens stimulus file for writing
        temp_stim = "{0}/stim.sp".format(OPTS.openram_temp)
        self.sf = open(temp_stim, "w")
        self.sf.write("* Stimulus for period of {0}n\n\n".format(period))

        # include files in stimulus file
        model_list = tech.spice["fet_models"] + [self.sram_sp_file]
        stimuli.write_include(stim_file=self.sf, models=model_list)

        # add vdd/gnd statements

        self.sf.write("* Global Power Supplies\n")
        stimuli.write_supply(stim_file=self.sf,
                             vdd_name=tech.spice["vdd_name"],
                             gnd_name=tech.spice["gnd_name"],
                             vdd_voltage=self.vdd,
                             gnd_voltage=self.gnd)

        # instantiate the sram
        self.sf.write("* Instantiation of the SRAM\n")
        stimuli.inst_sram(stim_file=self.sf,
                          abits=self.addr_size, 
                          dbits=self.word_size, 
                          sram_name=self.name)

        # create a buffer and an inverter
        self.sf.write("* Buffers and inverter Initialization\n")
        # FIXME: We should replace the clock buffer with the same
        # 2x buffer for control signals. This needs the buffer to be
        # added to the control logic though.
        stimuli.create_buffer(stim_file=self.sf,
                              buffer_name="clk1_buffer",
                              size=[1, 4])

        stimuli.create_buffer(stim_file=self.sf,
                              buffer_name="clk2_buffer",
                              size=[8, 16])

        stimuli.create_buffer(stim_file=self.sf,
                              buffer_name="buffer",
                              size=[1, 2])

        stimuli.create_inverter(stim_file=self.sf)

        # add a buffer for each signal and an inverter for WEb
        signal_list = []
        for i in range(self.word_size):
            signal_list.append("D[{0}]".format(i))
        for j in range(self.addr_size):
            signal_list.append("A[{0}]".format(j))
        for k in tech.spice["control_signals"]:
            signal_list.append(k)
        self.sf.write("* Buffers for each generated signal\n")
        stimuli.inst_buffer(stim_file=self.sf,
                           buffer_name="buffer",
                           signal_list=signal_list)
        stimuli.inst_buffer(stim_file=self.sf,
                           buffer_name="clk1_buffer",
                           signal_list=["clk"])
        stimuli.inst_buffer(stim_file=self.sf,
                           buffer_name="clk2_buffer",
                           signal_list=["clk_buf"])

        # add access transistors for data-bus
        self.sf.write("* Transmission Gates for data-bus and control signals\n")
        stimuli.inst_accesstx(stim_file=self.sf, dbits=self.word_size)

        # generate data and addr signals
        self.sf.write("* Generation of data and address signals\n")
        for i in range(self.word_size):
            if i == self.probe_data:
                stimuli.gen_data(stim_file=self.sf,
                                 clk_times=self.cycle_times,
                                 sig_name="D[{0}]".format(i),
                                 period=period,
                                 slew=self.slew)
            else:
                stimuli.gen_constant(stim_file=self.sf,
                                     sig_name="D[{0}]".format(i),
                                     v_val=self.gnd)

        stimuli.gen_addr(self.sf,
                         clk_times=self.cycle_times,
                         addr=self.probe_address,
                         period=period,
                         slew=self.slew)

        # generate control signals
        self.sf.write("* Generation of control signals\n")
        stimuli.gen_csb(self.sf, self.cycle_times, period, self.slew)
        stimuli.gen_web(self.sf, self.cycle_times, period, self.slew)
        stimuli.gen_oeb(self.sf, self.cycle_times, period, self.slew)
        stimuli.gen_web_trans(self.sf, self.cycle_times, period, self.slew)

        self.sf.write("* Generation of global clock signal\n")
        stimuli.gen_pulse(stim_file=self.sf,
                          sig_name="CLK",
                          v1=self.gnd,
                          v2=self.vdd,
                          offset=period,
                          period=period,
                          t_rise = self.slew,
                          t_fall = self.slew)
                          
        self.write_measures()

        # run until the last cycle time
        stimuli.write_control(self.sf,self.cycle_times[-1])

        self.sf.close()

    def write_measures(self):
        # meas statement for delay and power measurements
        self.sf.write("* Measure statements for delay and power\n")

        trig_name = tech.spice["clk"] + "_buf_buf"
        targ_name = "{0}".format("DATA[{0}]".format(self.probe_data))
        trig_val = targ_val = 0.5 * self.vdd
        # add measure statments for delay0
        stimuli.gen_meas_delay(stim_file=self.sf,
                               meas_name="DELAY0",
                               trig_name=trig_name,
                               targ_name=targ_name,
                               trig_val=trig_val,
                               targ_val=targ_val,
                               trig_dir="FALL",
                               targ_dir="FALL",
                               td=self.cycle_times[self.read0_cycle])

        stimuli.gen_meas_delay(stim_file=self.sf,
                               meas_name="DELAY1",
                               trig_name=trig_name,
                               targ_name=targ_name,
                               trig_val=trig_val,
                               targ_val=targ_val,
                               trig_dir="FALL",
                               targ_dir="RISE",
                               td=self.cycle_times[self.read1_cycle])
        
        # add measure statements for power
        t_initial = self.cycle_times[self.write0_cycle]
        t_final = self.cycle_times[self.write0_cycle+1]
        stimuli.gen_meas_power(stim_file=self.sf,
                               meas_name="WRITE0_POWER",
                               t_initial=t_initial,
                               t_final=t_final)

        t_initial = self.cycle_times[self.write1_cycle]
        t_final = self.cycle_times[self.write1_cycle+1]
        stimuli.gen_meas_power(stim_file=self.sf,
                               meas_name="WRITE1_POWER",
                               t_initial=t_initial,
                               t_final=t_final)
        
        t_initial = self.cycle_times[self.read0_cycle]
        t_final = self.cycle_times[self.read0_cycle+1]
        stimuli.gen_meas_power(stim_file=self.sf,
                               meas_name="READ0_POWER",
                               t_initial=t_initial,
                               t_final=t_final)

        t_initial = self.cycle_times[self.read1_cycle]
        t_final = self.cycle_times[self.read1_cycle+1]
        stimuli.gen_meas_power(stim_file=self.sf,
                               meas_name="READ1_POWER",
                               t_initial=t_initial,
                               t_final=t_final)
        



    def find_feasible_period(self,initial_period):
        """Uses an initial period and finds a feasible period before we
        run the binary search algorithm to find min period. We check if
        the given clock period is valid and if it's not, we continue to
        double the period until we find a valid period to use as a
        starting point. """

        feasible_period = initial_period
        time_out = 8
        while True:
            debug.info(1, "Finding feasible period: {0}ns".format(feasible_period))
            time_out -= 1

            if (time_out <= 0):
                debug.error("Timed out, could not find a feasible period.",2)

            (success, feasible_delay0, feasible_delay1)=self.try_feasible_period(feasible_period)                
            if not success:
                feasible_period = 2 * feasible_period
                continue

            debug.info(1, "Starting Binary Search Algorithm with feasible_period: {0}ns feasible_delay0/1 {1}ns/{2}ns".format(feasible_period,feasible_delay0,feasible_delay1))
            return (feasible_period, feasible_delay0, feasible_delay1)


    def try_feasible_period(self, period):
        """ This tries to simulate a period and checks if the result
        works. If so, it returns True and the feasible output delay."""

        # Checking from not data_value to data_value
        self.write_stimulus(period)
        stimuli.run_sim()
        delay0 = ch.convert_to_float(ch.parse_output("timing", "delay0"))
        delay1 = ch.convert_to_float(ch.parse_output("timing", "delay1"))        
        
        # if it failed or the read was longer than a period
        if type(delay0)!=float or type(delay1)!=float:
            return (False,0,0)
        else:
            delay0 *= 1e9
            delay1 *= 1e9
            debug.info(2,"Feasible period {0}, delay0={1}n delay1={2}ns".format(period,delay0,delay1))
        #key=raw_input("press return to continue")

        # The delay is from the negative edge for our SRAM
        return (True,delay0,delay1)



    def find_min_period(self,feasible_period,feasible_delay0, feasible_delay1):
        """Searches for the smallest period with output delays being within 5% of 
        long period. """

        previous_period = ub_period = feasible_period
        lb_period = 0.0

        # Binary search algorithm to find the min period (max frequency) of design
        time_out = 25
        while True:
            time_out -= 1
            if (time_out <= 0):
                debug.error("Timed out, could not converge on minimum period.",2)

            target_period = 0.5 * (ub_period + lb_period)
            debug.info(1, "MinPeriod Search: {0}ns (ub: {1} lb: {2})".format(target_period,
                                                                             ub_period,
                                                                             lb_period))

            if self.try_period(target_period, feasible_delay0, feasible_delay1):
                ub_period = target_period
            else:
                lb_period = target_period

            if ch.relative_compare(ub_period, lb_period, error_tolerance=0.05):
                # ub_period is always feasible
                return ub_period


    def try_period(self, period, feasible_delay0, feasible_delay1):
        """ This tries to simulate a period and checks if the result
        works. If it does and the delay is within 5% still, it returns True."""

        # Checking from not data_value to data_value
        self.write_stimulus(period)
        stimuli.run_sim()
        delay0 = ch.convert_to_float(ch.parse_output("timing", "delay0"))
        delay1 = ch.convert_to_float(ch.parse_output("timing", "delay1"))
        if type(delay0)==float:
            delay0 *= 1e9
        if type(delay1)==float:
            delay1 *= 1e9
        debug.info(2,"Period {0}, delay0={1}ns, delay1={2}ns".format(period,delay0, delay1))
        # if it failed or the read was longer than a period
        if type(delay0)!=float or type(delay1)!=float:
            return False
        else:
            if ch.relative_compare(delay1*1e9,feasible_delay1,error_tolerance=0.05):
                return False
            elif ch.relative_compare(delay0*1e9,feasible_delay0,error_tolerance=0.05):
                return False


        #key=raw_input("press return to continue")

        return True
    
    def set_probe(self,probe_address, probe_data):
        """ Probe address and data can be set separately to utilize other
        functions in this characterizer besides analyze."""
        self.probe_address = probe_address
        self.probe_data = probe_data

    def analyze(self,probe_address, probe_data, slew=tech.spice["rise_time"], load=1):
        """main function to calculate the min period for a low_to_high
        transistion and a high_to_low transistion returns a dictionary
        that contains all both the min period and associated delays
        Dictionary Keys: min_period1, delay1, min_period0, delay0
        """
        self.slew = slew
        self.load = load
        self.set_probe(probe_address, probe_data)

        (feasible_period, feasible_delay1, feasible_delay0) = self.find_feasible_period(tech.spice["feasible_period"])
        
        min_period = self.find_min_period(feasible_period, feasible_delay1, tech.spice["supply_voltage"])
        if min_period == None:
            return None
        debug.info(1, "Min Period: {0}n with a delay of {1}".format(min_period, feasible_delay1))
        
        read0_power=ch.convert_to_float(ch.parse_output("timing", "read0_power"))
        write0_power=ch.convert_to_float(ch.parse_output("timing", "write0_power"))
        read1_power=ch.convert_to_float(ch.parse_output("timing", "read1_power"))
        write1_power=ch.convert_to_float(ch.parse_output("timing", "write1_power"))

        data = {"min_period": min_period, 
                "delay1": feasible_delay1, 
                "delay0": feasible_delay0,
                "read_power": (read0_power+read1_power)/2,
                "write_power": (write0_power+write1_power)/2
                }
        return data


    def obtain_cycle_times(self, period):
        """Returns a list of key time-points [ns] of the waveform (each rising edge)
        of the cycles to do a timing evaluation. The last time is the end of the simulation
        and does not need a rising edge."""

        # idle cycle, no operation
        t_current = period 
        self.cycle_times = []
        
        # cycle0: W data 1 address 1111 to initialize cell to a value
        self.cycle_times.append(t_current)
        t_current += period

        # cycle1: W data 0 address 1111 (to ensure a write of value works)
        self.cycle_times.append(t_current)
        self.write0_cycle=1
        t_current += period
        
        # cycle2: W data 1 address 0000 (to clear the data bus cap)
        self.cycle_times.append(t_current)
        t_current += period

        # cycle3: R data 0 address 1111 to check W0 works
        self.cycle_times.append(t_current)
        self.read0_cycle=3
        t_current += period

        # cycle4: W data 1 address 1111 (to ensure a write of value works)
        self.cycle_times.append(t_current)
        self.write1_cycle=4
        t_current += period

        # cycle5: W data 0 address 0000 (to clear the data bus cap)
        self.cycle_times.append(t_current)
        t_current += period
        
        # cycle6: R data 1 address 1111 to check W1 works
        self.cycle_times.append(t_current)
        self.read1_cycle=6
        t_current += period

        # cycle7: wait a clock period to end the simulation
        self.cycle_times.append(t_current)
        t_current += period


