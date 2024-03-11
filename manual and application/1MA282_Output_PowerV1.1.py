###############################################################################
#
# Rohde & Schwarz GmbH & Co. KG
# 1MA282 Bluetooth for IoT
# Measuring Output Power
# Test Case: TP/TRM-LE/CA/BV-01-C
# Version 1.1
#
###############################################################################

import wx
import math as math

# Variables #

devc = 0                                                                        # Device Type
  
cfrq = 2.400E9                                                                  # Center Frequency
span = 0e0                                                                      # Frequency Span
rbw  = 3E6                                                                      # Resolution Bandwidth
vbw  = 3E6                                                                      # Video Bandwidth
att  = 20                                                                       # Input Attenuation in dB
rlev = 10                                                                       # Reference Level in dBm
tlvl = -30                                                                      # Trigger Level in dBm
coun = 100                                                                      # Number of Measurement Points
rnge = 0.8                                                                      # Measurement Range in 1/100 %

out  = ""                                                                       # Output String
rslt = []                                                                       # Result List

# Define Operating Frequencies #

if(devc == 0):                                                                  # 0 - Peripheral and Central Devices
    devc_name = "Peripheral and Central Device"
    devc_frq = [ 2.406E9, 2.440E9, 2.476E9 ]
if(devc == 1):                                                                  # 1 - Broadcaster and Observer Devices
    devc_name = "Broadcaster and Observer Device"
    devc_frq = [ 2.402E9, 2.426E9, 2.480E9 ]

# User Information #

print("----------------------------------------------------------------------------------")
print("1MA282 - Bluetooth for IoT")
print("Measuring Output Power")
print("TP/TRM-LE/CA/BV-01-C")  
print("----------------------------------------------------------------------------------")

# R&S FSV Configuration #

print("----------------------------------------------------------------------------------")
print("R&S FSV Configuration.")
print("----------------------------------------------------------------------------------")
    
FSV.write("*RST")                                                               # Reset Device
FSV.write("*CLS")                                                               # Clear Status
FSV.ask("*OPT?")                                                                # Option Query
FSV.ask("*IDN?")                                                                # Query ID String

time.sleep(0.5)

print("")

FSV.write("SYST:DISP:UPD ON")                                                   # Display Update ON during Remote Control.
FSV.write("INST:SEL SANALYZER")                                                 # Activate Spectrum Analyzer Mode.
FSV.write("INP:ATT " + str(att))                                                # Define Input Attenuation.
FSV.write("FREQ:CENT " + str(cfrq))                                             # Define Center Frequency.
FSV.write("FREQ:SPAN " + str(span))                                             # Define Frequency Span.
FSV.write("SENS:BAND " + str(rbw))                                              # Define Resolution Bandwidth.
FSV.write("SENS:BAND:VID " + str(vbw))                                          # Define Video Bandwidth.
FSV.write("SENS:WIND:DET POS")                                                  # Define Trace Detector.
FSV.write("DISP:TRAC:MODE WRIT")                                                # Define Trace Mode.
FSV.write("DISP:TRAC:Y:RLEV " + str(rlev))                                      # Define Reference Level.

print("----------------------------------------------------------------------------------")

# Measurement Section #

for k in range(0, 3): 

    # User Information:
    print("----------------------------------------------------------------------------------")
    print("Measurement for " + str(devc_name) + ".")
    print("Transmit Frequency @ " + str(devc_frq[k]/1E6) + " MHz.")
    print("----------------------------------------------------------------------------------")
        
    # Reset Analyzer Display:
    FSV.write("INIT:CONT ON")
    FSV.write("CALC:UNIT:POW DBM")
    FSV.write("CALC:MARK:AOFF")
    FSV.write("CALC:TLIN1:STAT OFF")
    FSV.write("CALC:TLIN2:STAT OFF")
     
    # Define Center Frequency:
    FSV.write("FREQ:CENT " +  str(devc_frq[k]))
    
    # Reset Trigger:
    FSV.write("TRIG:SOUR IMM")
    
    # Wait for Analyzer:
    time.sleep(0.25)
    
    # Set Triger to IF Power and define Trigger Level:
    FSV.write("TRIG:SOUR IFP")
    FSV.write("TRIG:LEV:IFP " + str(tlvl))

    # User Interaction:
    caption = '1MA282 - Bluetooth for IoT'
    message = 'Measuring Output Power\n\n' + 'Set Transmit Frequency for ' + str(devc_name) + ' to ' + str(devc_frq[k]/1E6)+ ' MHz. Adjust Sweep Time at Analyzer to display a complete Packet and use Trigger Offset to move Burst to the Center.\n\nPress OK to continue Measurement...'
    usrinpt = flib.wx.MessageDialog(_int.app.frame, message, caption, style=flib.wx.OK|flib.wx.CANCEL).ShowModal()
    
    # Verify, if User has set Operating Freq. and continue Measurement. Abort, if not!
    if(usrinpt != flib.wx.ID_OK):
        print("----------------------------------------------------------------------------------")
        print("Aborted!")
        print("----------------------------------------------------------------------------------")
        continue

    # User Information:
    print("")

    # Initiate Measurement and wait for Instrument:
    FSV.write("INIT:CONT OFF")
    FSV.write("INIT")
    FSV.ask("*OPC?")

    # Turn ON Marker 1:  
    FSV.write("CALC:MARK ON")
    
    # Query Power Level in dBm:   
    FSV.write("CALC:MARK:Y?") 
    
    # Read Power Level as Raw-String:
    peak = float(FSV.read())
    
    # User Information:
    print("")

    # Read Sweep Time and Trigger Offset:
    FSV.write("SWE:TIME?")
    swpt = float(FSV.read())
    FSV.write("TRIG:HOLD?")
    hold = float(FSV.read())

    # Activate Marker 1: 
    FSV.write("CALC:MARK1 ON")
    FSV.write("CALC:MARK1:X " + str(hold+0.25*(swpt)))
    if 'ledg' in globals():
        FSV.write("CALC:MARK1:X " + str(ledg))
    
    # Activate Marker 2:
    FSV.write("CALC:MARK2 ON")
    FSV.write("CALC:MARK2:X " + str(hold+0.75*(swpt)))
    if 'redg' in globals():
        FSV.write("CALC:MARK2:X " + str(redg))
    
    # User Interaction:
    caption = '1MA282 - Bluetooth for IoT'
    message = 'Measuring Output Power\n\n' + 'Output Power is calculated over ' + str(int(rnge*100)) + '% of the Duration of the Burst. Use Marker to define\nthe Measurement Range by positioning them at the lower and upper Edge of the Signal.\n\nPress OK to continue the Measurement...'
    usrinpt = flib.wx.MessageDialog(_int.app.frame, message, caption, style=flib.wx.OK|flib.wx.CANCEL).ShowModal()
    
    # Verify, if User has set Marker and continue Measurement. Abort, if not!
    if(usrinpt != flib.wx.ID_OK):
        print("----------------------------------------------------------------------------------")
        print("Aborted!")
        print("----------------------------------------------------------------------------------")
        continue
    
    # Query Burst Edges:  
    FSV.write("CALC:MARK1:X?")
    ledg = float(FSV.read())
    FSV.write("CALC:MARK1:STAT OFF")  
    FSV.write("CALC:MARK2:X?")
    redg = float(FSV.read())                                               
    FSV.write("CALC:MARK2:STAT OFF")
    
    # Calculate Burst Length and Limits:
    brst = redg - ledg
    llim = ledg+((1-rnge)/2*(brst))
    rlim = redg-((1-rnge)/2*(brst))
    
    # Turn Marker ON:
    FSV.write("CALC:MARK ON")  
    
    # Define Time Lines:
    FSV.write("CALC:TLIN1:STAT ON")
    FSV.write("CALC:TLIN2:STAT ON")
    FSV.write("CALC:TLIN1 " + str(llim))
    FSV.write("CALC:TLIN2 " + str(rlim))
    
    # Set Power Unit to Watt:
    FSV.write("CALC:UNIT:POW W")
    
    # User Information:
    print("")

    # Initialize Variables:
    lavg = pavg = pow = val = 0

    # Measurements:
    for n in range(0, coun+1):
        
        # Query Power Value at M[Freq]:
        frq = llim+n*(rlim-llim)/coun
        FSV.write("CALC:MARK:X " + str(frq)) 
        FSV.write("CALC:MARK:Y?")
        val = float(FSV.read())
        
        # Store Results in List:
        list = []
        list.append(frq)
        list.append(val)
        rslt.append(list)

        # Total Power:
        pow = pow + val

    # Average Power:
    pavg = pow / (coun + 1)
    
    # Calculate Output Power: 
    lavg = 10 * math.log((pavg/(1*math.pow(10,-3))),10)
    
    print("----------------------------------------------------------------------------------")

    # Initialize Result
    res = True

    # Generate Output String:
    out += "Transmit Frequency @ " + str(devc_frq[k]/1E6) + " MHz.\n"
    out += '{:>50}'.format("Min.Limit:")
    out += '{:>12}'.format("Max.Limit:")
    out += '{:>10}'.format("Margin:")
    out += '{:>10}'.format("Check:")
    out += "\n\n"
    
    out += '{:>25}'.format("Average Power [P_AVG]: ")
    out += '{:>11}'.format(format(lavg, '.2f') + " dBm")
    out += '{:>14}'.format("-20.00 dBm")
    out += '{:>12}'.format("20.00 dBm")

    if(lavg < -20 or lavg > +20):   
        if(lavg < -20):
            out += '{:>11}'.format((format(lavg+20, '.2f') + " dB"))
        if(lavg > +20):
            out += '{:>11}'.format((format(lavg-20, '.2f') + " dB"))           
        
        res = res and False
        out += '{:>9}'.format("FAILED")       
    
    else:
        
        if(lavg == 0):
            out += '{:>11}'.format((format(0, '.2f') + " dB"))
        if(lavg < 0):
            out += '{:>11}'.format((format(lavg+20, '.2f') + " dB"))
        if(lavg > 0):
            out += '{:>11}'.format((format(lavg-20, '.2f') + " dB"))                   
        
        res = res and True
        out += '{:>9}'.format("PASSED")

    out += "\n"
    out += '{:>25}'.format("Peak Power Value [P_PK]:")
    out += '{:>11}'.format(format(peak, '.2f') + " dBm")
    out += '{:>14}'.format("------ dBm")
    out += '{:>12}'.format(format(lavg+3, '.2f') + " dBm")
    
    out += '{:>11}'.format((format(peak-(lavg+3), '.2f') + " dB"))
    
    if(peak > lavg + 3):
        out += '{:>10}'.format("FAILED\n")
    else:
        out += '{:>10}'.format("PASSED\n")
        
    out += "\n"
    out += '{:>63}'.format("Test Result:")
    out += '{:>18}'.format("PASS" if(res) else "FAIL")
    out += "\n"
    out += "\n----------------------------------------------------------------------------------\n"

# Delete Global Variables #

if 'brst' in globals():
    del brst
if 'ledg' in globals() and 'redg' in globals():
    del ledg, redg
if 'llim' in globals() and 'rlim' in globals():
    del llim, rlim
if 'lavg' in globals() and 'pavg' in globals():
    del lavg, pavg

# Evaluation Print #

if(out != ""):
    print("----------------------------------------------------------------------------------")
    print("Evaluation for " + str(devc_name) + ".")
    print("Measured Results.")
    print("----------------------------------------------------------------------------------")
    print(out[:-1])
    print("Additional Information:")
    print("")
    print(" Maximum Limit for Average Power is set to 10 dBm if IUT is compliant to Core")
    print(" Specification V4.2 or earlier and not compliant to Core Specification Addendum 5.")
    print("")
    print("----------------------------------------------------------------------------------")
        
print("----------------------------------------------------------------------------------")
print("Executing Script finished!")
print("----------------------------------------------------------------------------------")