###############################################################################
#
# Rohde & Schwarz GmbH & Co. KG
# 1MA282 Bluetooth for IoT
# In-Band Emissions, at 1Ms/s and 2Ms/s
# Test Cases: TP/TRPM-LE/CA/BV-03-C, TP/TRM-LE/CA/BV-08-C
# Version 1.11
#
###############################################################################

import wx
import csv as csv
import math as math
import datetime as datetime
import matplotlib as mlib
import matplotlib.pyplot as plt

# Variables #

test = 0                                                                        # Test Type
devc = 0                                                                        # Device Type

cfrq = 2.401E9                                                                  # Center Frequency
span = 1.000E6                                                                  # Frequency Span
rbw  = 100E3                                                                    # Resolution Bandwidth
vbw  = 300E3                                                                    # Video Bandwidth
swp  = 10                                                                       # Number of Sweeps
swpt = 0.1                                                                      # Sweep Time in s
att  = 10.0                                                                     # Input Attenuation in dB
fchn = 0                                                                        # Number of first Channel
lchn = 80                                                                       # Number of last Channel

out  = ""                                                                       # Output String
rslt = []                                                                       # Result List

# Define Test Specifications #

if(test == 0):                                                                  # 0 - In-Band Emission, uncoded Data at 1 Ms/s
    test_name = "In-Band Emissions, uncoded Data at 1 Ms/s"
    test_id = "TP/TRM-LE/CA/BV-03-C"
    test_lim  = [ [ 1.0E6, 'N.D.' ], [ 2.0E6, -20 ], [ 3.0E6, -30] ]
if(test == 1):                                                                  # 1 - In-Band Emissions at 2 Ms/s
    test_name = "In-Band Emissions at 2 Ms/s"
    test_id = "TP/TRM-LE/CA/BV-08-C"
    test_lim  = [ [ 3.0E6, 'N.D.' ], [ 5.0E6, -20 ], [ 6.0E6, -30] ]

# Define Operating Frequencies #

if(devc == 0):                                                                  # 0 - Peripheral and Central Devices
    devc_name = "Peripheral and Central Device"
    devc_frq = [ 2.406E9, 2.440E9, 2.476E9 ]
if(devc == 1):                                                                  # 1 - Broadcaster and Observer Devices
    devc_name = "Broadcaster and Observer Device"
    devc_frq = [ 2.402E9, 2.426E9, 2.480E9 ]

# User Information #

print("----------------------------------------------------------------------")
print("1MA282 - Bluetooth for IoT")
print(test_name)
print(test_id)  
print("----------------------------------------------------------------------")

# R&S FSV Configuration #

print("----------------------------------------------------------------------")
print("R&S FSV Configuration.")
print("----------------------------------------------------------------------")
    
FSV.write("*RST")                                                               # Reset Device
FSV.write("*CLS")                                                               # Clear Status
FSV.ask("*OPT?")                                                                # Option Query
FSV.ask("*IDN?")                                                                # Query ID String

print("")

FSV.write("SYST:DISP:UPD ON")                                                   # Display Update ON during Remote Control
FSV.write("INST:SEL SANALYZER")                                                 # Activate Spectrum Analyzer Mode
FSV.write("INP:ATT " + str(att))                                                # Define Input Attenuation
FSV.write("FREQ:CENT " + str(cfrq))                                             # Define Center Frequency
FSV.write("FREQ:SPAN " + str(span))                                             # Define Frequency Span
FSV.write("SENS:BAND " + str(rbw))                                              # Define Resolution Bandwidth
FSV.write("SENS:BAND:VID " + str(vbw))                                          # Define Video Bandwidth 
FSV.write("CALC:UNIT:POW W")                                                    # Set Power Unit to W
FSV.write("SENS:WIND:DET AVER")                                                 # Define Trace Detector
FSV.write("DISP:TRAC:MODE MAXH")                                                # Define Trace Mode
FSV.write("SWE:TIME " + str(swpt))                                              # Define Sweep Time
FSV.write("SWE:COUN " + str(swp))                                               # Set Number of Sweeps
    
print("----------------------------------------------------------------------")

# Measurement Section #

for k in range(0, 3):

    caption = '1MA282 - Bluetooth for IoT'
    message = str(test_name) + '\n\nSet Transmit Frequency for ' + str(devc_name) + ' to ' + str(devc_frq[k]/1E6)+ ' MHz.\nPress OK to continue Measurement...'
    usrinpt = flib.wx.MessageDialog(_int.app.frame, message, caption, style=flib.wx.OK|flib.wx.CANCEL).ShowModal()
    
    print("----------------------------------------------------------------------")
    print("Measurement for " + str(devc_name) + ".")
    print("Transmit Frequency @ " + str(devc_frq[k]/1E6) + " MHz.")
    print("----------------------------------------------------------------------")

    # Verify, if User has set Operating Freq. and continue Measurement. Abort, if not!
    if(usrinpt != flib.wx.ID_OK):
        print("Aborted!")
        print("----------------------------------------------------------------------")
        continue
	
    FSV.write("FREQ:CENT " +  str(devc_frq[k]))                                 # Define Center Frequency
    FSV.write("SENS:ADJ:LEV")                                                   # Set RefLevel at Transmit Frequency
    FSV.ask("*OPC?")                                                            # Wait for Analyzer...

    for n in range(fchn, lchn + 1):
        
        temp  = [ ]
        ptx_n = 0.0
        
        # Calculate Frequency n in Hz
        freq_n = cfrq + n * 1E6

        print("----------------------------------------------------------------------")
        print("Frequency [" + str(n) + "]: " + str(freq_n / 1E6) + " MHz")
        print("----------------------------------------------------------------------")

        # Define Center Frequency
        FSV.write("FREQ:CENT " + str(freq_n))

        # Initiate Measurement and wait for Instrument
        FSV.write("INIT:CONT OFF")
        FSV.write("INIT")
        FSV.ask("*OPC?")
  
        # Turn ON Marker
        FSV.write("CALC:MARK ON")

        print("")

        for i in range(0, 10):
            
            # Calculate Frequency i in Hz
            freq_i = freq_n - 450e3 + i * 100e3
            
            # Move Marker to Frequency i
            FSV.write("CALC:MARK:X " + str(freq_i))

            # Query Power Level in W at Frequency i
            FSV.write("CALC:MARK:Y?") 
            
            # Read Power Level in W as Raw-String
            pow_i = float(FSV.read())
            
            # Calculate PTX in W
            ptx_n = ptx_n + pow_i
            
            # Store all Values in temporary List
            # Temp[List]: Transmit Freq., N, Freq_N, I, Freq_I, PowerLevel_I in W, Power_Level_I in dBm, PTX_N,I in W, PTX_N,I in dBm
            temp.append([devc_frq[k], n, freq_n, i, freq_i, pow_i, 10 * math.log(pow_i * math.pow(10, 3), 10), ptx_n, 10 * math.log(ptx_n * math.pow(10, 3), 10)])
       
        for j in range (0, len(temp)):
       
            # Add accumulated PTX Value to temporary List
            # Temp[List] extended by PTX_N in W and PTX_N in dBm
            temp[j].extend([ptx_n, 10 * math.log(ptx_n * math.pow(10, 3), 10)])
            
            # Add temporary List to Result List
            rslt.append(temp[j])

    print("----------------------------------------------------------------------")

# Writing Measurement Results to CSV #
    
datetime  = datetime.datetime.now()
directory = plib.application.get_script_dir() + str("\\") + str(test_id.replace("/", "_")) + str("\\")
filename  = datetime.strftime('%Y%m%d%H%M%S') + str("_")  + str(devc_name.replace(" ", "_"))

if not os.path.exists(directory):
    os.makedirs(directory)

output = str(directory) +  str(filename) + str(".csv")
        
with open(output, 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONE)
    writer.writerow(["FRQ_TX[Hz]", "N", "FRQ_N[Hz]", "I", "FRQ_I[Hz]", "POW_I[W]", "POW_I[dBm]", "PTX_N[W]", "PTX_N[dBm]", "PTX[W]", "PTX[dBm]"])
    writer.writerows(rslt)
	
# Evaluation Section #

# Get all relevant Values from Result List
rftx = [x[0]  for x in rslt][::10]
rfrq = [x[2]  for x in rslt][::10]
rptx = [x[10] for x in rslt][::10]

for o in range(len(rptx)/(lchn-fchn+1)):  
    
    # Initialize Result
    res = True

    # Devide Values based on Number of measured Channels and current Operating Frequency!
    ftx = rftx[o*(lchn-fchn+1):(o+1)*(lchn-fchn+1)]
    frq = rfrq[o*(lchn-fchn+1):(o+1)*(lchn-fchn+1)]
    ptx = rptx[o*(lchn-fchn+1):(o+1)*(lchn-fchn+1)]

    # User Information
    out += "Transmit Frequency @ " + str(float(ftx[0])/1E6) + " MHz.\n"
    out += '{:>47}'.format("Limit:")
    out += '{:>13}'.format("Margin:")
    out += '{:>10}'.format("Check:")
    out += "\n\n"

    for m in range(len(ptx)):
        
        # Set Limits for PTX
        if(frq[m] <= ftx[0] - test_lim[2][0] or  frq[m] >= ftx[0] + test_lim[2][0]): 
            r_lim = test_lim[2][1]
        if(frq[m] >= ftx[0] - test_lim[1][0] and frq[m] <= ftx[0] + test_lim[1][0]):
            r_lim = test_lim[1][1]
        if(frq[m] >= ftx[0] - test_lim[0][0] and frq[m] <= ftx[0] + test_lim[0][0]):
            r_lim = test_lim[0][1]
			
        # Check PTX for breaking Limits
        if(r_lim == "N.D."):
            res = res and True
            r_out = "------"
        if(r_lim != "N.D." and ptx[m] <= r_lim):
            res = res and True
            r_out = "PASSED"
        if(r_lim != "N.D." and ptx[m] >= r_lim):
            res = res and False
            r_out = "FAILED"
         
        # Generate Output String
        out += "PTX " + '{:>4}'.format("[" + str(fchn+m) + "]")
        out += " @ " + str(frq[m] / 1E6) + " MHz:"
        out += '{:>12}'.format(format(ptx[m], '.2f') + " dBm")
        out += '{:>15}'.format((format(r_lim, '.2f') + " dBm") if (r_lim != "N.D.") else "-- dBm")
        out += '{:>12}'.format((format(ptx[m]-r_lim, '.2f') + " dB") if (r_lim != "N.D.") else "-- dB")
        out += '{:>9}'.format(str(r_out)) + "\n"

    out += "\n----------------------------------------------------------------------\n"

    # Function for plotting Measurement Results as Bar Chart

    def do_plot():
        
        # Calculate Values for Plot: Frequencies in MHz, Invertet-PTX in dBm
        plt_frq = [x / 1E6 for x in frq]
        plt_ptx = [1E3 + x for x in ptx]
        
        # Calculate Values for Axis: xMin, xMax, yMin, yMax
        xmin = int(min(frq)/1E6-1)
        xmax = int(max(frq)/1E6+1)
        ymin = int(math.floor(min(ptx)/10))*10-5
        ymax = int(math.ceil(max(ptx)/10))*10+5
   
        # Toolbar: None
        mlib.rcParams['toolbar'] = 'None'

        # Configure Window
        fig = plt.figure()
        fig.canvas.set_window_title('1MA282 - Bluetooth for IoT')
        fig.suptitle(str(test_name) + "\n" + str(devc_name) + " @ " + str(float(ftx[0])/1E6) + " MHz", fontsize=14)
        fig.subplots_adjust(top=0.85)
        
        # General Plot Configuration
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('Frequency in [MHz]')
        ax.set_ylabel('PTX in [dBM]')  
        ax.axis([xmin, xmax, ymin, ymax])

        # Bar Plot
        ax.bar(plt_frq, plt_ptx, width=0.7, bottom=-1E3, align='center', alpha=0.5)
        
        # Horizontal Limit Lines - 20dBm
        ax.hlines(test_lim[1][1], min(plt_frq)-1, (ftx[0]-test_lim[0][0]-0.5E6)/1E6, color='red', linewidth=2)
        ax.hlines(test_lim[1][1], (ftx[0]+test_lim[0][0]+0.5e6)/1E6, max(plt_frq)+1, color='red', linewidth=2)              
        
        # Horizontal Limit Lines - 30dBm
        ax.hlines(test_lim[2][1], min(plt_frq)-1, (ftx[0]-test_lim[2][0]+0.5E6)/1E6, color='red', linewidth=2)
        ax.hlines(test_lim[2][1], (ftx[0]+test_lim[2][0]-0.5E6)/1E6, max(plt_frq)+1, color='red', linewidth=2)

        # Show Test Result
        if(res): ax.set_title('Test Result: PASSED!', fontsize=13)
        else: ax.set_title('Test Result: FAILED!', fontsize=13)

        # Show Plot
        plt.show()

    # Execute Plotting Function
    wx.CallAfter(do_plot)
    
    # Wait for generating Bar Chart
    time.sleep(0.3)

# User Information #

if(out != ""):
    print("----------------------------------------------------------------------")
    print("Evaluation for " + str(devc_name) + ".")
    print("Measured Results.")
    print("----------------------------------------------------------------------")
    print(out[:-1])
    
print("----------------------------------------------------------------------")
print("Executing Script finished!")
print("----------------------------------------------------------------------")