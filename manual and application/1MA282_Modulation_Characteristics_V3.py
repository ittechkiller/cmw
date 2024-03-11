###############################################################################
#
# Rohde & Schwarz GmbH & Co. KG
# 1MA282 Bluetooth for IoT
# Modulation CHracteristics & Carrier Frequency Offset and Drift, at 1Ms/s and 2Ms/s
# Test Cases: 
#   RF-PHY/TRM-LE/CA/BV-05-C - Modulation Characteristics, uncoded data at 1 Ms/s
#   RF-PHY/TRM-LE/CA/BV-06-C - Carrier frequency offset and drift, uncoded data at 1 Ms/s
#   RF-PHY/TRM-LE/CA/BV-09-C - Stable Modulation Characteristics, uncoded data at 1 Ms/s
#   RF-PHY/TRM-LE/CA/BV-10-C - Modulation Characteristics at 2 Ms/s
#   RF-PHY/TRM-LE/CA/BV-11-C - Stable Modulation Characteristics at 2 Ms/s
#   RF-PHY/TRM-LE/CA/BV-12-C - Carrier frequency offset and drift at 2 Ms/s
#   RF-PHY/TRM-LE/CA/BV-13-C - Modulation Characteristics, LE Coded (S=8)
#   RF-PHY/TRM-LE/CA/BV-14-C - Carrier frequency offset and drift, LE Coded (S=8)
#
# Version 3.0
# - replaced numpy.fromstring with numpy.frombuffer
# - corrected delta f2 99.9% calculation
# - added narrow channel filter on raw IQ data (matching CMW filtering)
# - for filtering the raw IQ data is downloaded from the instrument, filtered 
#   and uploaded via an iq.tar file
###############################################################################


# Variables #

Freq    = 2440E6 #BT channel
SymRate = 2E6    #Symbol Rate (1MHz or 2MHz)
Coded   = False  #False - no coding, True - S=8
N       = 50    #statistics count
SweTime = 10e-3  #Capture Time in sec 

CorrThresh = 80  #correlation threshold for the pattern search (Sync Word)
###############################################################################

import numpy as np
import os

np.seterr(invalid='ignore')

###############################################################################

# Function definitions

#----------------------------------------
def Calc_df1avg(ov,x,Estart,ELen,FreqAbs,mode):
    #compute df1avg 
    
    Nseq = int(ELen / 8)
    
    if (mode == 0):
        #11110000 pattern (uncoded)
        offset = 4
        Nend = Nseq-1
    else:
        #00001111 pattern (coded)
        #use full eval range
        offset = 31
        Nend = Nseq+1           

    #bits 2,3,6 and 7
    bitpos = [1,2,5,6]

    idx = np.zeros(2,dtype=int)
    jdx = np.zeros(2,dtype=int)
    f1ccf = np.zeros(Nend)
    df1max = np.zeros([Nend,len(bitpos)])
    
    for i in range(0,Nend): #don't use th last 4 bits
        #start f1ccf at the fifth bit (oversampled samples are left and  right of the symbol)
        tmp = (x == Estart + offset + 8*i - 0.5)
        idx[0] = np.argwhere(tmp == True)
        tmp = (x == Estart + offset + 8*(i+1)-1 + 0.5)
        idx[1] = np.argwhere(tmp == True)+1
        
        #fig = flib.pltframe.figure
        #fig.clear()
        #ax = fig.add_subplot(211, autoscale_on = False, xlim = (x[idx[0]], x[idx[1]]), ylim = (-300000,300000))
        #ax.plot(x[idx[0]:idx[1]], FreqAbs[idx[0]:idx[1]])
        #flib.pltframe.Show()
         
        #step 5
        f1ccf[i] = np.mean(FreqAbs[idx[0]:idx[1]])
        
        #step 6
        for j in range(0,len(bitpos)):
            jdx[0] = idx[0] + ov * bitpos[j]
            jdx[1] = idx[0] + ov * (bitpos[j]+1)
            df1max[i,j]  = np.abs(np.mean(FreqAbs[jdx[0]:jdx[1]]) - f1ccf[i]);
            
            #ax = fig.add_subplot(212, autoscale_on = False, xlim = (x[idx[0]], x[idx[1]]), ylim = (-300000,300000))
            #ax.plot(x[jdx[0]:jdx[1]], FreqAbs[jdx[0]:jdx[1]])
            #flib.pltframe.Show()
            
    df1avg = np.mean(np.mean(df1max))
    df1max_out = df1max.flatten()
            
    #step 7: mean of df1max per packet
    return df1avg, df1max_out

#------------------------------------------------------------------------------

def Calc_df2avg(ov,x,Estart,ELen,FreqAbs):

    #compute df2avg - use 10101010 pattern

    Nseq = int(ELen / 8)
    
    idx = np.zeros(2,dtype=int)
    jdx = np.zeros(2,dtype=int)
    f2ccf  = np.zeros(Nseq-1)
    df2max = np.zeros((Nseq-1,8))
    
    for i in range(0,Nseq-1): #don't use the last 4 bits
        #start f2ccf at the fifth bit (oversampled samples are left and  right of the symbol)
        tmp = (x == Estart + 4 + 8*i - 0.5)
        idx[0] = np.argwhere(tmp == True)
        tmp = (x == Estart + 4 + 8*(i+1)-1 + 0.5)
        idx[1] = np.argwhere(tmp == True) + 1 # +1 because a:b only goes to b-1
        
        #step 10
        f2ccf[i] = np.mean(FreqAbs[idx[0]:idx[1]])
        
        #step 11
        for j in range(0,8):
            jdx[0] = idx[0] + ov * j
            jdx[1] = idx[0] + ov * (j+1) + 1
            df2max[i,j]  = np.max(np.abs(FreqAbs[jdx[0]:jdx[1]] - f2ccf[i]))

    #step 12
    df2avg = np.mean(np.mean(df2max))
    df2max_out = df2max.flatten()
    
    return df2avg, df2max_out

#------------------------------------------------------------------------------

def Calc_fn(ov,x,Estart,ELen,FreqAbs):

    #differentiate betwen LE1 and LE2
    FSV.write('DDEM:SRAT?')
    SymRate = float(FSV.read())
    if (SymRate == 1e6):
        Nbit = [8, 10]; #LE1
    elif (SymRate == 2e6):
        Nbit = [16, 20]; #LE2

    Nseq = int(np.floor(ELen / Nbit[1]))

    #K70 removes carrier freq offset from FreqAbs trace.
    #use K70 result Carrier Frequency Error and add it back to the trace
    #this is relevant for FreqOff result
    FSV.write('CALC2:MARK:FUNC:DDEM:STAT:CFER?')
    CFE = float(FSV.read())
    FreqAbs = FreqAbs + CFE

    #step 4
    f0 = np.mean(FreqAbs[0:(Nbit[0]*ov)])

    idx = np.zeros(2,dtype=int)
    jdx = np.zeros(2,dtype=int)
    fn  = np.zeros(Nseq)
    dfn = np.zeros(Nseq-5)
    
    for i in range(0,Nseq):
        #start f1 at the 2nd bit (oversampled samples are left and  right of the symbol)       
        tmp = (x == Estart + 1 + Nbit[1]*i - 0.5)
        idx[0] = np.argwhere(tmp == True)         
        tmp = (x == Estart + 1 + Nbit[1]*(i+1)-1 + 0.5)
        idx[1] = np.argwhere(tmp == True) + 1
        
        #step 6
        fn[i] = np.mean(FreqAbs[idx[0]:idx[1]])

    #CMW: Freq. Offset
    n = np.argmax(np.abs(fn))
    FreqOff = fn[n]

    #CMW: Freq. Drift
    #|f0-fn| <= 50kHz
    df0n = f0 - fn[1:]
    n = np.argmax(np.abs(df0n))
    FreqDrift = df0n[n]

    #CMW: Initial freq. drift
    #|f1-f0| <= 23kHz
    InitFreqDrift = fn[0] - f0

    #CMW: Max. drift rate
    #|fn-fn-5| <=20kHz
    for i in range(5,Nseq):
        dfn[i-5] = fn[i] - fn[i-5]
    
    n = np.argmax(np.abs(dfn))
    MaxDriftRate = dfn[n]

    #no result for frequency accuracy

    return FreqOff, FreqDrift, InitFreqDrift, MaxDriftRate


#------------------------------------------------------------------------------

def Calc_fn_Coded(ov,x,Estart,ELen,FreqAbs):

    Nbit = [16, 16]
    offset = [2, 27]
        
    Nseq = int(np.floor(ELen / Nbit[1]))

    #K70 removes carrier freq offset from FreqAbs trace.
    #use K70 result Carrier Frequency Error and add it back to the trace
    #this is relevant for FreqOff result
    FSV.write('CALC2:MARK:FUNC:DDEM:STAT:CFER?')
    CFE = float(FSV.read())
    FreqAbs = FreqAbs + CFE

    #step 4
    f0123 = np.zeros(4)
    for i in range(0,4):
        tmp = (x == offset[0] + Nbit[0]*i - 0.5)
        istart = np.argwhere(tmp == True)[0,0] 
        tmp = (x == offset[0] + Nbit[0]*(i+1)-1 + 0.5)
        istop = np.argwhere(tmp == True)[0,0] + 1
        f0123[i] = np.mean(FreqAbs[istart:istop])

    idx = np.zeros(2,dtype=int)
    jdx = np.zeros(2,dtype=int)
    fn  = np.zeros(Nseq)
    dfn = np.zeros(Nseq-3)
    
    for i in range(0,Nseq):
        #start f1 at the 27th bit (oversampled samples are left and  right of the symbol)       
        tmp = (x == Estart + offset[1] + Nbit[1]*i - 0.5)
        idx[0] = np.argwhere(tmp == True)         
        tmp = (x == Estart + offset[1] + Nbit[1]*(i+1)-1 + 0.5)
        idx[1] = np.argwhere(tmp == True) + 1
        
        #step 6
        fn[i] = np.mean(FreqAbs[idx[0]:idx[1]])

    #CMW: Freq. Offset
    n = np.argmax(np.abs(fn))
    FreqOff = fn[n]

    #CMW: Freq. Drift
    #|f0-fn| <= 50kHz
    df0n = f0123[0] - fn[1:]
    n = np.argmax(np.abs(df0n))
    FreqDrift = df0n[n]

    #CMW: Initial freq. drift
    #|f3-f0| <= 19.2kHz
    InitFreqDrift = f0123[3] - f0123[0]

    #CMW: Max. drift rate
    #|fn-fn-3| <=19.2kHz
    for i in range(3,Nseq):
        dfn[i-3] = fn[i] - fn[i-3]
    
    n = np.argmax(np.abs(dfn))
    MaxDriftRate = dfn[n]

    #no result for frequency accuracy

    return FreqOff, FreqDrift, InitFreqDrift, MaxDriftRate
#------------------------------------------------------------------------------

def getTraceBinary(SCPIstr):
    
    #store format setting
    form = FSV.query('FORM?')
    
    #change to binary format
    FSV.write('FORM REAL,32')

    #send SCPI command
    FSV.write(SCPIstr)
    #get binary data as string
    tmp = FSV.read(report_mode=ReportMode.Off)

    #number of length digits
    NDigits = int(tmp[1])
    
    #number of bytes
    NBytes =  int(tmp[2:2+NDigits])
        
    #convert to float32
    data = np.frombuffer(tmp[2+NDigits:], dtype="float32")
    
    FSV.write('FORM ' + form)
    
    return data
#------------------------------------------------------------------------------

def filter_iq_data(SymRate):
     
    if (SymRate == 1E6):
        #define 1.3MHz channel filter
        #CMW SW filter with fs = 32MHz
        b = np.array([3.2716e-16,-3.6052e-16,3.8977e-16,-4.2659e-16,4.6825e-16,-5.1267e-16,5.6521e-16,-6.2022e-16,6.8693e-16,-7.5901e-16,8.4155e-16,-9.3389e-16,1.0416e-15,-1.1611e-15,1.3015e-15,-1.46e-15,1.6449e-15,-1.8566e-15,2.1038e-15,-2.3906e-15,2.7254e-15,-3.1192e-15,3.582e-15,-4.13e-15,4.7793e-15,-5.5573e-15,6.4859e-15,-7.6085e-15,8.963e-15,-1.0616e-14,1.2635e-14,-1.5125e-14,1.8206e-14,-2.2052e-14,2.6879e-14,-3.2989e-14,4.0778e-14,-5.0798e-14,6.3803e-14,-8.0855e-14,1.0345e-13,-1.3376e-13,1.7495e-13,-2.3175e-13,3.1139e-13,-4.2518e-13,5.9136e-13,-8.4036e-13,1.2252e-12,-1.8431e-12,2.8854e-12,-4.7631e-12,8.4794e-12,-1.6978e-11,4.1757e-11,-1.5587e-10,1.7587e-09,2.2665e-08,1.1376e-07,3.852e-07,1.0611e-06,2.6193e-06,5.9765e-06,1.278e-05,2.5594e-05,4.7684e-05,8.2562e-05,0.0001331,0.00020073,0.0002853,0.00038535,0.00049902,0.00062499,0.00076247,0.00091086,0.0010684,0.0012308,0.0013915,0.001541,0.0016689,0.0017649,0.0018193,0.0018232,0.0017685,0.0016466,0.0014495,0.0011698,0.00080218,0.00034447,-0.00020149,-0.00082958,-0.0015297,-0.0022879,-0.0030869,-0.0039056,-0.0047193,-0.0054993,-0.0062134,-0.0068266,-0.0073025,-0.0076044,-0.0076967,-0.0075457,-0.0071208,-0.0063953,-0.0053473,-0.0039608,-0.002227,-0.00014537,0.0022759,0.0050193,0.0080584,0.011358,0.014876,0.018561,0.022356,0.026199,0.030022,0.033757,0.037333,0.04068,0.043734,0.046434,0.048723,0.050556,0.051894,0.052708,0.052982,0.052708,0.051894,0.050556,0.048723,0.046434,0.043734,0.04068,0.037333,0.033757,0.030022,0.026199,0.022356,0.018561,0.014876,0.011358,0.0080584,0.0050193,0.0022759,-0.00014537,-0.002227,-0.0039608,-0.0053473,-0.0063953,-0.0071208,-0.0075457,-0.0076967,-0.0076044,-0.0073025,-0.0068266,-0.0062134,-0.0054993,-0.0047193,-0.0039056,-0.0030869,-0.0022879,-0.0015297,-0.00082958,-0.00020149,0.00034447,0.00080218,0.0011698,0.0014495,0.0016466,0.0017685,0.0018232,0.0018193,0.0017649,0.0016689,0.001541,0.0013915,0.0012308,0.0010684,0.00091086,0.00076247,0.00062499,0.00049902,0.00038535,0.0002853,0.00020073,0.0001331,8.2562e-05,4.7684e-05,2.5594e-05,1.278e-05,5.9765e-06,2.6193e-06,1.0611e-06,3.852e-07,1.1376e-07,2.2665e-08,1.7587e-09,-1.5587e-10,4.1757e-11,-1.6978e-11,8.4794e-12,-4.7631e-12,2.8854e-12,-1.8431e-12,1.2252e-12,-8.4036e-13,5.9136e-13,-4.2518e-13,3.1139e-13,-2.3175e-13,1.7495e-13,-1.3376e-13,1.0345e-13,-8.0854e-14,6.3803e-14,-5.0797e-14,4.0778e-14,-3.2988e-14,2.688e-14,-2.2052e-14,1.8207e-14,-1.5124e-14,1.2636e-14,-1.0616e-14,8.9636e-15,-7.6074e-15,6.4871e-15,-5.5565e-15,4.7803e-15,-4.1292e-15,3.5826e-15,-3.1188e-15,2.7257e-15,-2.3902e-15,2.1035e-15,-1.857e-15,1.6446e-15,-1.4601e-15,1.3006e-15,-1.1628e-15,1.04e-15,-9.3522e-16,8.3962e-16,-7.5993e-16,6.8415e-16,-6.2315e-16,5.619e-16,-5.1518e-16,4.6531e-16,-4.2976e-16,3.8787e-16,-3.6058e-16,3.26e-16])
    else:
        #define 2.6MHz channel filter
        b = np.array([-1.949e-16,1.9734e-16,-2.2191e-16,2.2852e-16,-2.5195e-16,2.632e-16,-2.8722e-16,3.0412e-16,-3.2753e-16,3.5009e-16,-3.7548e-16,4.0315e-16,-4.3132e-16,4.6652e-16,-4.9607e-16,5.3733e-16,-5.7551e-16,6.2321e-16,-6.6712e-16,7.2251e-16,-7.766e-16,8.4384e-16,-9.0824e-16,9.8486e-16,-1.0645e-15,1.1561e-15,-1.2523e-15,1.3612e-15,-1.4776e-15,1.61e-15,-1.7528e-15,1.9098e-15,-2.0848e-15,2.2771e-15,-2.4912e-15,2.726e-15,-2.9912e-15,3.2807e-15,-3.6091e-15,3.9714e-15,-4.3787e-15,4.8326e-15,-5.3439e-15,5.9188e-15,-6.5646e-15,7.2957e-15,-8.1203e-15,9.0603e-15,-1.0124e-14,1.1343e-14,-1.2732e-14,1.433e-14,-1.6165e-14,1.8289e-14,-2.0747e-14,2.361e-14,-2.6955e-14,3.0881e-14,-3.5512e-14,4.0996e-14,-4.7538e-14,5.5377e-14,-6.4834e-14,7.6316e-14,-9.0362e-14,1.0767e-13,-1.2917e-13,1.5613e-13,-1.9024e-13,2.3384e-13,-2.9018e-13,3.6382e-13,-4.6126e-13,5.9191e-13,-7.696e-13,1.0149e-12,-1.3593e-12,1.8511e-12,-2.5673e-12,3.6319e-12,-5.2514e-12,7.778e-12,-1.1834e-11,1.8557e-11,-3.0125e-11,5.0914e-11,-9.0283e-11,1.698e-10,-3.4416e-10,7.7065e-10,-1.988e-09,6.3929e-09,-3.0611e-08,4.0522e-07,5.3908e-06,2.7201e-05,9.6963e-05,0.00026641,0.00057125,0.0010004,0.0015272,0.0021378,0.0027837,0.0033381,0.003638,0.0035356,0.0028969,0.001601,-0.00040686,-0.0030627,-0.006177,-0.0094415,-0.012428,-0.014605,-0.015391,-0.014239,-0.01069,-0.0044484,0.0045578,0.016122,0.029756,0.044714,0.060045,0.074663,0.087465,0.097442,0.10378,0.10596,0.10378,0.097442,0.087465,0.074663,0.060045,0.044714,0.029756,0.016122,0.0045578,-0.0044484,-0.01069,-0.014239,-0.015391,-0.014605,-0.012428,-0.0094415,-0.006177,-0.0030627,-0.00040686,0.001601,0.0028969,0.0035356,0.003638,0.0033381,0.0027837,0.0021378,0.0015272,0.0010004,0.00057125,0.00026641,9.6963e-05,2.7201e-05,5.3908e-06,4.0522e-07,-3.0611e-08,6.3929e-09,-1.988e-09,7.7065e-10,-3.4416e-10,1.698e-10,-9.0283e-11,5.0914e-11,-3.0125e-11,1.8557e-11,-1.1834e-11,7.778e-12,-5.2514e-12,3.6319e-12,-2.5673e-12,1.8511e-12,-1.3593e-12,1.0149e-12,-7.696e-13,5.9191e-13,-4.6126e-13,3.6382e-13,-2.9018e-13,2.3384e-13,-1.9024e-13,1.5613e-13,-1.2917e-13,1.0767e-13,-9.0362e-14,7.6316e-14,-6.4839e-14,5.5379e-14,-4.7539e-14,4.0998e-14,-3.5509e-14,3.0881e-14,-2.6953e-14,2.3611e-14,-2.0747e-14,1.829e-14,-1.6165e-14,1.4331e-14,-1.2732e-14,1.1343e-14,-1.0125e-14,9.0604e-15,-8.1202e-15,7.2961e-15,-6.5638e-15,5.9186e-15,-5.3441e-15,4.8329e-15,-4.3792e-15,3.9714e-15,-3.6086e-15,3.2824e-15,-2.9903e-15,2.7281e-15,-2.4905e-15,2.2764e-15,-2.0839e-15,1.9086e-15,-1.7521e-15,1.6084e-15,-1.4785e-15,1.3591e-15,-1.2522e-15,1.155e-15,-1.065e-15,9.8458e-16,-9.0832e-16,8.4278e-16,-7.777e-16,7.2306e-16,-6.6746e-16,6.2322e-16,-5.7497e-16,5.3662e-16,-4.9848e-16,4.6489e-16,-4.3248e-16,4.0208e-16,-3.7633e-16,3.4917e-16,-3.2816e-16,3.0474e-16,-2.8517e-16,2.6533e-16,-2.4785e-16,2.3282e-16,-2.1493e-16,2.0704e-16,-1.8872e-16])
        
    #get reference level
    RefLev = FSV.ask('DISP:TRAC:Y:RLEV?');
    #convert from dBm to Volts
    RefLev = np.sqrt(50.0/1000) * 10**(float(RefLev[0])/20.0)
    
    #get IQ data from instrument
    FSV.query(':CALC3:FORM RIM; *OPC?')
    tmpR = np.array(getTraceBinary((':TRAC3:DATA? TRACE1R')))
    tmpI = np.array(getTraceBinary((':TRAC3:DATA? TRACE1I')))
    
    iq = tmpR + tmpI * 1j
    
    #scale iq data to Volts
    iq = RefLev * iq
    
    #filter
    iq_len = iq.size
    iq_filt = np.concatenate((iq,iq,iq),axis=0)
    iq_filt = np.convolve(b,iq_filt,'same')
    iq_filt = iq_filt[iq_len:2*iq_len]
    
    iqfile = 'tmp.iq.tar'  
    
    #save iq.tar file
    WriteIqTar( iq_filt, 32e6, iqfile)
    
    #upload iq file to instrument
    upload_file(iqfile,'C:\\Temp\\' + iqfile)
    
    #delete iq.tar file
    os.remove(iqfile)
    
    #load iq file into K70
    if (instr.find('FSV') >= 0):
        FSV.write('MMEM:LOAD:IQ:STAT 1, \'C:\\Temp\\' + iqfile + '\'')
        FSV.ask('*OPC?')
    else:
        FSV.write(':INP:FILE:PATH \'C:\\Temp\\' + iqfile + '\'')
        FSV.ask(':INP:SEL FIQ; *OPC?')
        
        #run single to update results
        FSV.ask(':INIT:IMM; *OPC?')

        
#------------------------------------------------------------------------------

def upload_file(file_source,file_destintion):
    #function to upload a file stored on the remote PC to a spectrum analyzer

    filesize = os.path.getsize(file_source)
    
    BlockSize = FSV.instrument.chunk_size - 1024;
    
    FSV.instrument.send_end = False

    header = '#' + str(len(str(filesize))) + str(filesize);
    FSV.write('MMEM:DATA \'' + file_destintion + '\',' + header)

    iFile = open(file_source,mode='rb');
    
    #the number of bytes remaining is initially set to the total number of bytes to send
    Blocks = 0;
    remaining = filesize;

    while remaining > 0:
        nextblock = remaining;
        if(nextblock > BlockSize):
            nextblock = BlockSize;
    
        # set the new number of bytes remaining
        remaining = remaining - nextblock;
        Blocks = Blocks + 1;
    
        rawdata = iFile.read(nextblock);
    
        # last block to send
        if (remaining==0):
            FSV.instrument.send_end = True
        
        #send binary data
        FSV.write(rawdata,report_mode=ReportMode.Off)
        
    iFile.close()

#------------------------------------------------------------------------------


def WriteIqTar( iqData, fs, FileName):
    #Writes an iq.tar file. Complex iqData values are interpreted as Volts.
    #iqData can be a list of complex or list of floats (iqiqiq format).
    
    import tarfile
    import os
    import re
    
    path,filename = os.path.split(FileName)

    #Create binary file
    binaryfile = re.sub( "iq.tar", "complex.1ch.float32", filename, flags=re.IGNORECASE)
    NumberOfSamples = WriteIqw( iqData, os.path.join(path, binaryfile))
    if NumberOfSamples == 0:
        return 0
    
    xmlfilename = re.sub( "iq.tar", "xml", filename, flags=re.IGNORECASE)
    __WriteXml( fs, NumberOfSamples, binaryfile, os.path.join(path, xmlfilename))
    
    try:
        tar = tarfile.open( FileName, "w")
        tar.add( os.path.join(path, binaryfile), arcname=binaryfile)
        tar.add( os.path.join(path, xmlfilename), arcname=xmlfilename)
        tar.close()
        os.remove( os.path.join(path, binaryfile))
        os.remove( os.path.join(path, xmlfilename))
    except:
        print("IqTar (" + FileName +") write error!" )
        return 0
    
    return NumberOfSamples

#------------------------------------------------------------------------------


def WriteIqw( iqData, FileName):
    #writes an IQW file (file of binary floats).
    
    import struct

    #check if iqData is complex
    if isinstance(iqData[0], complex):
        iqData = Complex2Iqiq( iqData)
       
    NumberOfSamples = len(iqData) // 2
        
    try:
        file = open( FileName, "wb")
        file.write( struct.pack("f"*len(iqData),*iqData))
        file.close
    except:
        print("File (" + FileName +") write error!" )
        return 0
    
    return NumberOfSamples

#------------------------------------------------------------------------------ 

def __WriteXml( fs, NumberOfSamples, filenameiqw, filenamexml):
    #Function to write the xml part of the iq.tar
    
    from datetime import datetime
    
    xmlfile = open ( filenamexml, "w")
    
    xmlfile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    xmlfile.write("<?xml-stylesheet type=\"text/xsl\" href=\"open_IqTar_xml_file_in_web_browser.xslt\"?>\n")
    xmlfile.write("<RS_IQ_TAR_FileFormat fileFormatVersion=\"2\" xsi:noNamespaceSchemaLocation=\"http://www.rohde-schwarz.com/file/RsIqTar.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n")
    xmlfile.write("<Name>Python iq.tar Writer (iqdata.py)</Name>\n")
    xmlfile.write("<DateTime>"+ datetime.now(None).isoformat() +"</DateTime>\n")
    xmlfile.write("<Samples>" + str(NumberOfSamples) + "</Samples>\n")
    xmlfile.write("<Clock unit=\"Hz\">" + str(fs) + "</Clock>\n")
    xmlfile.write("<Format>complex</Format>\n")
    xmlfile.write("<DataType>float32</DataType>\n")
    xmlfile.write("<ScalingFactor unit=\"V\">1</ScalingFactor>\n")
    xmlfile.write("<NumberOfChannels>1</NumberOfChannels>\n")
    xmlfile.write("<DataFilename>" + filenameiqw+ "</DataFilename>\n")
    xmlfile.write("</RS_IQ_TAR_FileFormat>\n")
    xmlfile.close()

    
    return

#------------------------------------------------------------------------------

def Complex2Iqiq( complexList):
    #Returns a list of I/Q samples from a complex list.
    
    import numpy as np
    
    iqiqiqList = np.empty(2*len(complexList))
    iqiqiqList[0::2] = complexList.real
    iqiqiqList[1::2] = complexList.imag  

    return iqiqiqList

#------------------------------------------------------------------------------

    
###############################################################################

#Main program

#change directory to script location
os.chdir(os.path.dirname(__file__))
        
print("----------------------------------------------------------------------------------")
print("1MA282 - Bluetooth for IoT\n")
print("Measuring Modulation Characteristics\n")
print("RF-PHY/TRM-LE/CA/BV-05-C [Modulation Characteristics, uncoded data at 1 Ms/s]")
print("RF-PHY/TRM-LE/CA/BV-06-C [Carrier frequency offset and drift, uncoded data at 1 Ms/s]")
print("RF-PHY/TRM-LE/CA/BV-09-C [Stable Modulation Characteristics, uncoded data at 1 Ms/s]")
print("RF-PHY/TRM-LE/CA/BV-10-C [Modulation Characteristics at 2 Ms/s]")
print("RF-PHY/TRM-LE/CA/BV-11-C [Stable Modulation Characteristics at 2 Ms/s]")
print("RF-PHY/TRM-LE/CA/BV-12-C [Carrier frequency offset and drift at 2 Ms/s]")
print("RF-PHY/TRM-LE/CA/BV-13-C [Modulation Characteristics, LE Coded (S=8)]")
print("RF-PHY/TRM-LE/CA/BV-14-C [Carrier frequency offset and drift, LE Coded (S=8)]")
print("----------------------------------------------------------------------------------\n")

FSV.connect()
FSV.instrument.send_end = True

#instrument indentification
idn = FSV.ask('*IDN?')

if (idn[1].find('FSV') >= 0):
    instr = 'FSV';
else:
    instr = 'FSW_FPS';

#basic settings    
FSV.ask('*RST; *CLS; *OPC?')
FSV.write('FREQ:CENT '  + str(Freq))

#open VSA
if (instr.find('FSV') >= 0):
    FSV.ask('INST:SEL DDEM; *OPC?')
else:
    FSV.ask('INST:CRE DDEM,\'VSA\'; *OPC?')
  
#load preset    
if (Coded == False): #uncoded
    if (SymRate == 1E6):
        FSV.ask('DDEM:PRES \'Bluetooth_5_LE1M\'; *OPC?')
    else:
        FSV.ask('DDEM:PRES \'Bluetooth_5_LE2M\'; *OPC?')
else: #coded
    SymRate = 1E6 #LE Coded always uses 1MHz
    FSV.ask('DDEM:PRES \'Bluetooth_5_LE_S8\'; *OPC?')

#increase oversampling rate to 32 to have 32MHz sample rate (CMW setting)
if (SymRate == 1E6):
    FSV.ask(":SENS:DDEM:PRAT 32; *OPC?")
else: #2MHz
    FSV.ask(":SENS:DDEM:PRAT 16; *OPC?")
    
FSV.ask('ADJ:LEV; *OPC?')
FSV.write('INIT:CONT OFF')
FSV.write('DDEM:RLEN ' + str(SweTime))
FSV.write('SWE:COUN 1')

FSV.write('DDEMod:SEARch:SYNC:IQCThreshold ' + str(CorrThresh))

#turn off freq drift compensation. BT spec assumes it to be in the FreqAbs trace
FSV.write('DDEM:NORM:CFDR OFF')

#32 oversampling of frequency deviation trace
ov = 32;
FSV.write('DISP:WIND1:PRAT:AUTO MAN')
FSV.write('DISP:WIND1:PRAT ' +str(ov));

#meas only if burst was found
FSV.write('DDEM:SEAR:BURS:MODE BURS')

#capture oversampling
FSV.write('DDEM:PRAT?')
PRAT = float(FSV.read())
    
#Eval Range
FSV.write('CALC:ELIN1:VAL?')
Estart = float(FSV.read())
FSV.write('CALC:ELIN2:VAL?')
Estop  = float(FSV.read())
ELen   = Estop - Estart + 1/PRAT;
Estart = int(Estart)
    
#FreqAbs trace start/stop in symbols
FSV.write('DISP:TRAC1:X:STAR?')
Tstart = float(FSV.read())
FSV.write('DDEM:TIME?')
Tstop = float(FSV.read())

Nseq = int(ELen / 8)

cnt = [0, 0, 0, 0, 0]
df1avg = np.empty(N) * np.nan
df2avg = np.empty(N) * np.nan
df1max = np.empty(((Nseq+1)*4,N)) * np.nan
df2max = np.empty(((Nseq-1)*8,N)) * np.nan
FreqOff       = np.empty(N) * np.nan
FreqDrift     = np.empty(N) * np.nan
InitFreqDrift = np.empty(N) * np.nan
MaxDriftRate  = np.empty(N) * np.nan

for i in range(0,N): 
            
    print("Statistic Count = {}/{}".format(i+1,N))
    
    FSV.ask('INIT:IMM; *OPC?')    
   
    SYNC = FSV.query('STATus:QUEStionable:SYNC:CONDition?')
    
    #FSW return binary value as string, 
    if (instr == 'FSW' or instr == 'FPS'):
        SYNC = flib.bin2dec(SYNC)
    else:
        #FSV returns decimal as string
        SYNC = int(SYNC)
        
    #check for proper signal sync
    if (SYNC > 0):
        #no burst found or pattern not found
        cnt[4] = cnt[4]+1
        
    else:
        
        #apply channel filter to raw IQ data
        filter_iq_data(SymRate)
        
        #get Frequency Deviation trace
        tmp = np.array(getTraceBinary(('TRAC1? TRACE1')))
        FreqAbs = tmp.astype(np.float)
        Ntrace = FreqAbs.size
        
        #time vector
        step = (Tstart + Tstop/Ntrace)
        x = np.arange(Tstart, (Tstop - Tstop/Ntrace) + step, step)
        
        #get demodulated bitstream
        tmp = np.array(FSV.ask('TRAC4? TRACE1', report_mode=ReportMode.Off))
        BitStream = tmp.astype(np.int)
        
        #find payload pattern
        #bitstream starts at symbol 0
        if Coded == False:
            if (np.sum(np.abs(BitStream[(Estart):(Estart+8)] - [1,1,1,1,0,0,0,0])) == 0):
                mode = 0
            elif (np.sum(np.abs(BitStream[Estart:(Estart+8)] - [1,0,1,0,1,0,1,0])) == 0):
                mode = 1
            else: 
                mode = 3 #unknown pattern
        else: #S=8 payload
            if (np.sum(np.abs(BitStream[(31+Estart):(31+Estart+8)] - [0,0,0,0,1,1,1,1])) == 0):
                mode = 2
            else:
                mode = 3   
                
        cnt[mode] = cnt[mode]+1

        # RF-PHY/TRM-LE/CA/BV-05-C [Modulation Characteristics, uncoded data at 1 Ms/s]
        # RF-PHY/TRM-LE/CA/BV-09-C [Stable Modulation Characteristics, uncoded data at 1 Ms/s]
        # RF-PHY/TRM-LE/CA/BV-10-C [Modulation Characteristics at 2 Ms/s]
        # RF-PHY/TRM-LE/CA/BV-11-C [Stable Modulation Characteristics at 2 Ms/s]        
        # calculations apply to 1 and 2MHz mode
        
        if (mode == 0):
           #11110000 pattern (uncoded)           
                         
           #compute df1avg - use 11110000 pattern
           df1avg[i], _ = Calc_df1avg(ov,x,Estart,ELen,FreqAbs,mode)                  
        elif (mode == 1):
            #10101010 pattern
                   
            df2avg[i], df2max[:,i] = Calc_df2avg(ov,x,Estart,ELen,FreqAbs)
                        
            # RF-PHY/TRM-LE/CA/BV-06-C [Carrier frequency offset and drift, uncoded data at 1 Ms/s]
            # RF-PHY/TRM-LE/CA/BV-12-C [Carrier frequency offset and drift at 2 Ms/s]        
            FreqOff[i], FreqDrift[i], InitFreqDrift[i], MaxDriftRate[i] = Calc_fn(ov,x,Estart,ELen,FreqAbs)        
        elif (mode == 2):
            #00001111 pattern (coded S8)
            
            # RF-PHY/TRM-LE/CA/BV-13-C [Modulation Characteristics, LE Coded (S=8)]
            # RF-PHY/TRM-LE/CA/BV-14-C [Carrier frequency offset and drift, LE Coded (S=8)]
            df1avg[i], df1max[:,i] = Calc_df1avg(ov,x,Estart,ELen,FreqAbs,mode)
            FreqOff[i], FreqDrift[i], InitFreqDrift[i], MaxDriftRate[i] = Calc_fn_Coded(ov,x,Estart,ELen,FreqAbs) 
        
        inp = FSV.ask(':INP:SEL?')
        if inp[0] == 'FIQ':
            FSV.ask(':INP:SEL RF; *OPC?')
            
FSV.close()

#CMW: Freq Dev delta f1_avg
#     Freq Dev delta f1_max/min not required by spec
#stats of df1_avg over all packets
if (cnt[0] > 0) | (cnt[2] > 0):
    df1_avg_out = [np.nanmean(df1avg), np.nanmax(df1avg), np.nanmin(df1avg)]
else:
    df1_avg_out = np.tile(np.nan,3)

#CMW: Freq Dev delta f2_avg
#      Freq Dev delta f2_max/min not required by spec
#stats of df2_avg over all packets
if (cnt[1] > 0):
    df2_avg_out = [np.nanmean(df2avg), np.nanmax(df2avg), np.nanmin(df2avg)]
else:
    df2_avg_out = np.tile(np.nan,3)

#CMW: delta f1 99.9%
if (cnt[2] > 0):    
    df1max_999 = np.nanpercentile(df1max.flatten(),99.9)
else:
    df1max_999 = np.nan

#CMW: delta f2 99.9%
if (cnt[1] > 0):    
    #df2max_999 = np.nanpercentile(df2max.flatten(),99.9)
    df2max_999 = np.nanpercentile(df2max.flatten(),0.1,interpolation='lower')
else:
    df2max_999 = np.nan

#CMW: min(delta f2_avg) / max(delta f1_avg)
dfavg_coef = df2_avg_out[2] / df1_avg_out[1]

#CMW: Freq Offset
if (cnt[1] > 0) | (cnt[2] > 0):
    n_max = np.nanargmax(np.abs(FreqOff))
    FreqOff_out = [np.nanmean(FreqOff), FreqOff[n_max], np.nanstd(FreqOff)]
else:
    FreqOff_out = np.tile(np.nan,3)
    
#CMW: Freq Drift
if (cnt[1] > 0) | (cnt[2] > 0):
    n_max = np.nanargmax(np.abs(FreqDrift))
    FreqDrift_out = [np.nanmean(FreqDrift), FreqDrift[n_max], np.nanstd(FreqDrift)]
else:
    FreqDrift_out = np.tile(np.nan,3)
    
#CMW: Initial Freq Drift
if (cnt[1] > 0) | (cnt[2] > 0):
    n_max = np.nanargmax(np.abs(InitFreqDrift))
    InitFreqDrift_out = [np.nanmean(InitFreqDrift), InitFreqDrift[n_max], np.nanstd(InitFreqDrift)]
else:
    InitFreqDrift_out = np.tile(np.nan,3)
    
#CMW: Max Drift Rate
if (cnt[1] > 0) | (cnt[2] > 0):
    n_max = np.nanargmax(np.abs(MaxDriftRate))
    MaxDriftRate_out  = [np.nanmean(MaxDriftRate), MaxDriftRate[n_max], np.nanstd(MaxDriftRate)]
else:
    MaxDriftRate_out = np.tile(np.nan,3)
    

###############################################################################

# Generate Output String:

greek_delta = unichr(916).encode('utf-8')
print("\n----------------------------------------------------------------------------------")
print("Frequency                 {}MHz".format(Freq/1E6))
print("Symbol Rate               {}MHz".format(SymRate/1e6))
print("Coded (S=8)               {}".format(Coded))
print("Statistic Count           {}\n".format(N))

if (Coded == False):
    print("Pattern                   11110000  10101010 Other no_sync")
    print("                          {:>}{:>10}{:>8}{:>9}\n".format(cnt[0], cnt[1], cnt[3], cnt[4]))
    print("{}f2 99.9% [kHz]           {:>.1f}".format(greek_delta, df2max_999/1e3))
    print("{}f2avg/{}f1avg             {:>.2f}\n".format(greek_delta, greek_delta, dfavg_coef))
    print("                          Average   Max      Std")
    print("Freq Offset [kHz]         {:>.1f} {:>9.1f} {:>8.1f}".format(FreqOff_out[0]/1e3, FreqOff_out[1]/1e3, FreqOff_out[2]/1e3))
    print("Freq Drift [kHz]          {:>.1f} {:>9.1f} {:>7.1f}".format(FreqDrift_out[0]/1e3, FreqDrift_out[1]/1e3, FreqDrift_out[2]/1e3))
    print("Initial Freq Drift [kHz]  {:>.1f} {:>9.1f} {:>8.1f}".format(InitFreqDrift_out[0]/1e3, InitFreqDrift_out[1]/1e3, InitFreqDrift_out[2]/1e3))
    print("Max Drift Rate [kHz/50us] {:>.1f} {:>9.1f} {:>8.1f}\n".format(MaxDriftRate_out[0]/1e3, MaxDriftRate_out[1]/1e3, MaxDriftRate_out[2]/1e3))
    print("                          Average   Max      Min")
    print("Freq Dev {}f1avg [kHz]     {:>.1f} {:>9.1f} {:>7.1f}".format(greek_delta, df1_avg_out[0]/1e3, df1_avg_out[1]/1e3, df1_avg_out[2]/1e3))
    print("Freq Dev {}f2avg [kHz]     {:>.1f} {:>9.1f} {:>7.1f}".format(greek_delta, df2_avg_out[0]/1e3, df2_avg_out[1]/1e3, df2_avg_out[2]/1e3))
else:
    print("Pattern                   11111111 Other no_sync")
    print("                          {:>}{:>10}{:>8}\n".format(cnt[2], cnt[3], cnt[4]))
    print("{}f1 99.9% [kHz]           {:>.1f}".format(greek_delta, df1max_999/1e3))    
    print("                          Average   Max      Std")
    print("Freq Offset [kHz]         {:>.1f} {:>9.1f} {:>8.1f}".format(FreqOff_out[0]/1e3, FreqOff_out[1]/1e3, FreqOff_out[2]/1e3))
    print("Freq Drift [kHz]          {:>.1f} {:>9.1f} {:>7.1f}".format(FreqDrift_out[0]/1e3, FreqDrift_out[1]/1e3, FreqDrift_out[2]/1e3))
    print("Initial Freq Drift [kHz]  {:>.1f} {:>9.1f} {:>8.1f}".format(InitFreqDrift_out[0]/1e3, InitFreqDrift_out[1]/1e3, InitFreqDrift_out[2]/1e3))
    print("Max Drift Rate [kHz/50us] {:>.1f} {:>9.1f} {:>8.1f}\n".format(MaxDriftRate_out[0]/1e3, MaxDriftRate_out[1]/1e3, MaxDriftRate_out[2]/1e3))
    print("                          Average   Max      Min")
    print("Freq Dev {}f1avg [kHz]     {:>.1f} {:>9.1f} {:>7.1f}".format(greek_delta, df1_avg_out[0]/1e3, df1_avg_out[1]/1e3, df1_avg_out[2]/1e3))
print("----------------------------------------------------------------------------------")