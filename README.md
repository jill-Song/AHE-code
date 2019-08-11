# AHE-code
This code is about the neural mass model AHE-CM, which is applied to simulate the generalized periodic discharges, triphasic waves, in acute hepatic encephalopathy (AHE). 
Outline:

   Folder structure
   Usage
   
Folder structure:
   mmc6.m: the computational model AHE-CM, which is a set of differential equations. The four parameters of AHE-CM are variables and other parameters are set to be their nominal values.  
   % Inputs(tau_rec_e:  recovery time for excitatory synapses, tau_rec_i:  recovery time for inhibitory synapses, ltp_factor: amplification factor, propofol:   control parameter relating to the decay time of the IPSP)
   % Ouputs:  simulated EEG signal
   
   GetPara.m: function for obtaining different model outputs when the four parameters are set to be different values in their own physiological ranges. 
   
   
We also refer the code of Ruijter2017_model. 
% Ruijter, Barry Johannes, et al. "Synaptic damage underlies EEG abnormalities in postanoxic encephalopathy: A computational study." Clinical neurophysiology 128.9 (2017): 1682-1695.
