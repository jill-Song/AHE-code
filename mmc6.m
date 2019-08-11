%% Function to generate Generalized Periodic Discharges in Acute Hepatic Encephalopathy(AHE)
% This code is about a new computational modl of acute hepatic encephalopathy. 
% It is constructed to study the Generalized Periodic Discharges in Acute Hepatic Encephalopathy(AHE)
% the inpus are four paramters in equations corresponding to three AHE mechanisms 
%
% Inputs:   tau_rec_e:  recovery time for excitatory synapses
%           tau_rec_i:  recovery time for inhibitory synapses
%           ltp_factor: amplification factor 
%           propofol:   control parameter relating to the decay time of the IPSP
%
% Ouputs:  simulated EEG signal
% Reference: code Ruijter2017_model
% Ruijter, Barry Johannes, et al. 
% "Synaptic damage underlies EEG abnormalities in postanoxic encephalopathy: A computational study." 
% Clinical neurophysiology 128.9 (2017): 1682-1695.

%%
function EEG=mmc6(parset_tau_rec_e,parset_tau_rec_i,parset_ltp_factor,parset_propofol)

		
		tau_rec_e	=	parset_tau_rec_e;
		tau_rec_i	=	parset_tau_rec_i;
		ltp_factor	=	parset_ltp_factor;
		propofol	=	parset_propofol;

  
    % 0. General settings
    sim_time=60;           % simulation time (sec)
    steps=8960;            % simulation steps per second / 'sample frequency' (Hz)
    white_noise=1;         % 0= no noise input, % 1= white noise input

    % 1. Define model parameters
	
		% part a: basic Liley model parameters
		p.S_e_max 		= 	0.5; 	% maximum spike-rate of e-neurons (1/ms) (0.05-0.5 / ms)
		p.S_i_max 		= 	0.5; 	% maximum spike-rate of i-neurons (1/ms) (0.05-0.5 / ms)
		p.h_e_r			= 	-70; 	% resting potential of e-neurons (mV) (-80:-60 mV)
		p.h_i_r			= 	-70; 	% resting potential of i-neurons (mV) (-80:-60 mV)
		p.mu_e 			= 	-50; 	% spike-threshold for e-neurons (mV) (-55:-40)
		p.mu_i 			= 	-50; 	% spike-threshold for i-neurons (mV) (-55:-40)
		p.sigma_e 		= 	5; 		% standarddeviation of spike-thresholds for e-neurons (mV) (2-7)
		p.sigma_i 		= 	5; 		% standarddeviation of spike-thresholds for i-neurons (mV) (2-7)
		p.tau_e 		= 	94; 	% membrane time-constant of e-neurons (ms) (5-150)
		p.tau_i 		= 	42; 	% membrane time-constant of i-neurons (ms) (5-150)
		p.h_ee_eq 		= 	45; 	% reversal potential of excitatory synaptic currents (mV) (-20:10) !!! In artikel: 40!!!
		p.h_ei_eq 		= 	45; 	% reversal potential of excitatory synaptic currents (mV) (-20:10) !!! In artikel: 40!!!
		p.h_ie_eq 		= 	-90; 	% reversal potential of inhibitory synaptic currents (mV)(-90:h_k_r-5)
		p.h_ii_eq 		= 	-90; 	% reversal potential of inhibitory synaptic currents (mV)(-90:h_k_r-5)
		p.Gamma_ee 		= 	0.71; 	% peak amplitude of EPSP (mV)(0.1-2.0)
		p.Gamma_ei 		= 	0.71; 	% peak amplitude of EPSP (mV)(0.1-2.0)
		p.Gamma_ie 		= 	0.71; 	% peak amplitude of IPSP (mV) (0.1-2.0)
		p.Gamma_ii 		= 	0.71; 	% peak amplitude of IPSP (mV) (0.1-2.0)
		p.gamma_ee 		= 	0.3; 	% EPSP rate constant (1/ms) (0.1-1.0)
		p.gamma_ei 		= 	0.3; 	% EPSP rate constant (1/ms) (0.1-1.0)
		p.gamma_ie 		= 	0.065; 	% 65*10^(-3); %IPSP rate constant (1/ms) (0.01-0.5)
		p.gamma_ii 		= 	0.065; 	% 65*10^(-3); %IPSP rate constant (1/ms) (0.01-0.5)
		p.p_ee 			= 	3.460; 	% nonspecific firing-rate to e-neurons (1/ms) (0-10)
		p.p_ee_sd 		= 	1.000; 	% standard-deviation of fluctuations in p
		p.p_ei 			= 	5.070; 	% nonspectific firing-rate to i-neurons (1/ms) (0-10)
		p.p_ei_sd 		= 	0; 		% standard-deviation of fluctuations in p
		p.p_ie			=	0; 		% nonspecific firing-rate from i to e-neurons (1/ms) (0-10)
		p.p_ii			=	0; 		% nonspecific firing-rate from i to e-neurons (1/ms) (0-10)
		p.N_ei_b 		= 	3000; 	% number of synaptic contacts from e-neurons to i-neurons (2000-5000)
		p.N_ee_b 		= 	3000; 	% number of synaptic contacts from e-neurons to e-neurons (2000-5000)
		p.N_ie_b 		= 	500; 	% number of synaptic contacts from i-neurons to e-neurons (100-1000)
		p.N_ii_b 		= 	500; 	% number of synaptic contacts from i-neurons to i-neurons (100-1000)
        
        p.N_ee_a		=	0; 		% (not used!)  Number of excitatory cortico-cortical synapses (1000-5000)
        p.N_ei_a		=	0;		% (not used!)  Number of excitatory cortico-cortical synapses (1000-5000)
        p.Lambda_ee		=	0; 		% (not used!)  Decay rate cortico-cortical connectivity (1/mm) (0.01-0.1)
        p.Lambda_ei		=	0; 		% (not used!)  Decay rate cortico-cortical connectivity (1/mm) (0.01-0.1)
        p.nu_ee			=	0; 		% (not used!)  Axonal conduction velocity (0.1-1 mm ms^-1)
        p.nu_ei			=	0; 		% (not used!)  Axonal conduction velocity (0.1-1 mm ms^-1)
		
		% part b: additional parameters to model activity-dependent synaptic depression
		p.tau_rec_e		=	tau_rec_e;  % recovery time constant for e-synapses
		p.tau_rec_i		=	tau_rec_i;  % recovery time constant for i-synapses 
		p.rho_e			=	0.003;      % synaptic depletion factor for e-synapses 
		p.rho_i			=	0.003;      % synaptic depletion factor for i-synapses 

		% part c: parameter adaptations to model anoxic Long Term Potentiation (LTP)
		p.ltp_factor	=	ltp_factor; % long-term-potentiation factor
		p.Gamma_ee		=	(1+p.ltp_factor)*p.Gamma_ee;
		p.Gamma_ei	=	(1+p.ltp_factor)*p.Gamma_ei;
        
        p.Gamma_ee_rest	=	p.Gamma_ee; % resting value of maximum PSP (e --> e)
		p.Gamma_ei_rest	=	p.Gamma_ei; % resting value of maximum PSP (e --> i)
		p.Gamma_ie_rest	=	p.Gamma_ie; % resting value of maximum PSP (i --> e)
		p.Gamma_ii_rest	=	p.Gamma_ii; % resting value of maximum PSP (i --> i)
	
		% part c: parameter adaptations to model effects of anesthesia 
		p.propofol		=	propofol;	% propofol concentration (no units)
		[p]				=	liley_anesthetics(p);

    % 4. Calculate equillibrium potentials
    [p.v_e_equil,p.v_i_equil]=liley_calculate_equilibrium(p);

    % 5. Define initial values, start in equillibrium state
    h 		= 	1000/steps;         % time step (msec)
    T 		= 	sim_time*10^3;      % observation time (ms)
    N 		= 	T/h-1;            % number of time steps
   
    
    X 		= 	zeros(18,N);        % initialization of state vector
    X(1,1) 	= 	p.v_e_equil;     
    X(2,1) 	= 	p.v_i_equil;     
    X(3,1)	= 	exp(p.gamma_ee/p.gamma_ee_0)/p.gamma_ee*p.Gamma_ee*(p.N_ee_b*S_e(p.v_e_equil) + p.N_ee_a*S_e(p.v_e_equil)+p.p_ee);
    X(5,1)	= 	exp(p.gamma_ei/p.gamma_ei_0)/p.gamma_ei*p.Gamma_ei*(p.N_ei_b*S_e(p.v_e_equil) + p.N_ei_a*S_e(p.v_e_equil)+p.p_ei);
    X(7,1)	= 	exp(p.gamma_ie/p.gamma_ie_0)/p.gamma_ie*p.Gamma_ie*(p.N_ie_b*S_i(p.v_i_equil));
    X(9,1)	= 	exp(p.gamma_ii/p.gamma_ii_0)/p.gamma_ii*p.Gamma_ii*(p.N_ii_b*S_i(p.v_i_equil));
    X(11,1)	=	p.N_ee_a*S_e(p.v_e_equil); 
    X(13,1)	=	p.N_ei_a*S_e(p.v_e_equil); 
    X(15,1) = 	p.Gamma_ee/(1+p.tau_rec_e*p.rho_e*S_e(p.v_e_equil));       
    X(16,1) = 	p.Gamma_ei/(1+p.tau_rec_e*p.rho_e*S_e(p.v_e_equil));   
    X(17,1) = 	p.Gamma_ie/(1+p.tau_rec_i*p.rho_i*S_i(p.v_i_equil));        
    X(18,1) = 	p.Gamma_ii/(1+p.tau_rec_i*p.rho_i*S_i(p.v_i_equil));     

    % 6. Start simulation
%       noise = importdata('noise.mat');
%     noise = zeros(18,N); 

    for n=1:N

        noise = zeros(18,1); 
        if white_noise==1
         noise(4)= p.gamma_ee_tilde*exp(p.gamma_ee/p.gamma_ee_0)*X(15,n)*p.p_ee_sd*randn(1,1);
           
        end

        X(:,n+1) = X(:,n) + h.*dynamics(X(:,n))+sqrt(h).*noise;  
    end
 
    
    % 7. downsample and filter data
    EEG=-X(1,:); 				% minus sign since negative=upward

    %% I. Function definitions 
	
	function dX = dynamics(X)

        dX = zeros(18,1);

        % Calculate synaptic reversal potentials
        psi_ee=(p.h_ee_eq-X(1))/abs(p.h_ee_eq-p.h_e_r);
        psi_ie=(p.h_ie_eq-X(1))/abs(p.h_ie_eq-p.h_e_r);
        psi_ei=(p.h_ei_eq-X(2))/abs(p.h_ei_eq-p.h_i_r);
        psi_ii=(p.h_ii_eq-X(2))/abs(p.h_ii_eq-p.h_i_r);

        % Calculate synaptic inputs A_jk 
        A_ee=p.N_ee_b*S_e(X(1))+X(11)+p.p_ee;
        A_ei=p.N_ei_b*S_e(X(1))+X(13)+p.p_ei;
        A_ie=p.N_ie_b*S_i(X(2));
        A_ii=p.N_ii_b*S_i(X(2));   

        % Calculate state vector
        dX(1) = (1/p.tau_e)*(p.h_e_r-X(1)+psi_ee*X(3)+psi_ie*X(7)); % V_e
        dX(2) = (1/p.tau_i)*(p.h_i_r-X(2)+psi_ei*X(5)+psi_ii*X(9)); % V_i
        dX(3) = X(4); % I_ee
        dX(4) = -(p.gamma_ee+p.gamma_ee_tilde)*X(4)     -p.gamma_ee*p.gamma_ee_tilde*X(3)   +p.gamma_ee_tilde*exp(p.gamma_ee/p.gamma_ee_0)*X(15)*A_ee; % J_ee
        dX(5) = X(6); % I_ei
        dX(6) = -(p.gamma_ei+p.gamma_ei_tilde)*X(6)     -p.gamma_ei*p.gamma_ei_tilde*X(5)   +p.gamma_ei_tilde*exp(p.gamma_ei/p.gamma_ei_0)*X(16)*A_ei; % J_ei
        dX(7) = X(8); % I_ie
        dX(8) = -(p.gamma_ie+p.gamma_ie_tilde)*X(8)     -p.gamma_ie*p.gamma_ie_tilde*X(7)   +p.gamma_ie_tilde*exp(p.gamma_ie/p.gamma_ie_0)*X(17)*A_ie; % J_ie
        dX(9) = X(10); % I_ii
        dX(10) = -(p.gamma_ii+p.gamma_ii_tilde)*X(10)   -p.gamma_ii*p.gamma_ii_tilde*X(9)   +p.gamma_ii_tilde*exp(p.gamma_ii/p.gamma_ii_0)*X(18)*A_ii; % J_ii
        dX(11) = X(12); %theta_ee
        dX(12) = -2*p.nu_ee*p.Lambda_ee*X(12)           -p.nu_ee^2*p.Lambda_ee^2*X(11)      +p.nu_ee^2*p.Lambda_ee^2*p.N_ee_a*S_e(X(1)); %d_theta_ee
        dX(13) = X(14); %theta_ei
        dX(14) = -2*p.nu_ei*p.Lambda_ei*X(14)           -p.nu_ei^2*p.Lambda_ei^2*X(13)      +p.nu_ei^2*p.Lambda_ei^2*p.N_ei_a*S_e(X(1)); %d_theta_ei
        dX(15) = (p.Gamma_ee_rest-X(15))/p.tau_rec_e-p.rho_e*S_e(X(1))*X(15); % Gamma_ee
        dX(16) = (p.Gamma_ei_rest-X(16))/p.tau_rec_e-p.rho_e*S_e(X(1))*X(16); % Gamma_ei   
        dX(17) = (p.Gamma_ie_rest-X(17))/p.tau_rec_i-p.rho_i*S_i(X(2))*X(17); % Gamma_ie 
        dX(18) = (p.Gamma_ii_rest-X(18))/p.tau_rec_i-p.rho_i*S_i(X(2))*X(18); % Gamma_ii    
    end

    function spikerate = S_e(v)    
        spikerate = p.S_e_max./(1 + exp(-sqrt(2)*(v - p.mu_e)./p.sigma_e));
    end

    function spikerate = S_i(v)
        spikerate = p.S_i_max./(1 + exp(-sqrt(2)*(v - p.mu_i)./p.sigma_i));
    end
	
	function [p2]=liley_anesthetics(p1)

		amplitude_effect=1; % Effect on maximum amplitude PSPs on (1) or off (0)
		decaytime_effect=1; % Effect on decay time iPSPs on (1) or off (0)

		c=p1.propofol; % c=propofol concentration

		% I. calculate PSP maximum amplitudes
		if amplitude_effect==1
			factor_e=0.707^2.22/(0.707^2.22+c^2.22);
			factor_i=(0.79^2.6+0.56*c^2.6)/(0.79^2.6+c^2.6);
		else
			factor_e=1;
			factor_i=1;
		end

		p2.Gamma_ee_0=p1.Gamma_ee;
		p2.Gamma_ei_0=p1.Gamma_ei;
		p2.Gamma_ie_0=p1.Gamma_ie;
		p2.Gamma_ii_0=p1.Gamma_ii;

		p2.Gamma_ee=factor_e*p1.Gamma_ee;
		p2.Gamma_ei=factor_e*p1.Gamma_ei;
		p2.Gamma_ie=factor_i*p1.Gamma_ie;
		p2.Gamma_ii=factor_i*p1.Gamma_ii;

		% II. Calculate PSP decay times (for inhibitory synapses only)
		p2.gamma_ie_0=      p1.gamma_ie;
		p2.gamma_ii_0=      p1.gamma_ii;

		if decaytime_effect==1 && c~=0
			kappa_i=            (0.32^2.7+4.7*c^2.7)/(0.32^2.7+c^2.7);
			epsilon_i=          exp(2.5466-1.3394*kappa_i)*sqrt(kappa_i-1)+...
								(exp(-1.2699*(kappa_i-1))-1)*(1/kappa_i^2+....
								lambertw(-1,exp(1)^(-0.23630/kappa_i^2)/(1-3.1462*kappa_i)));
			p2.epsilon_ie=      epsilon_i;
			p2.epsilon_ii=      epsilon_i;
			p2.gamma_ie=        p2.epsilon_ie*p2.gamma_ie_0/(exp(p2.epsilon_ie)-1);
			p2.gamma_ie_tilde=  p2.gamma_ie*exp(p2.epsilon_ie);
			p2.gamma_ii=        p2.epsilon_ii*p2.gamma_ii_0/(exp(p2.epsilon_ii)-1);
			p2.gamma_ii_tilde=  p2.gamma_ii*exp(p2.epsilon_ii);
		else
			p2.gamma_ie=        p1.gamma_ie;
			p2.gamma_ie_tilde=	p1.gamma_ie;
			p2.gamma_ii=        p1.gamma_ii;
			p2.gamma_ii_tilde=  p1.gamma_ii;
		end

		p2.gamma_ee_0=      p1.gamma_ee;
		p2.gamma_ei_0=      p1.gamma_ei;
		p2.gamma_ee_tilde=      p1.gamma_ee;
		p2.gamma_ei_tilde=      p1.gamma_ei;

		% III. Update parameter set
		p1=rmfield(p1,{'Gamma_ee','Gamma_ei','Gamma_ie','Gamma_ii','gamma_ie','gamma_ii'});
		names = [fieldnames(p1); fieldnames(p2)]; 
		p2 = cell2struct([struct2cell(p1); struct2cell(p2)], names, 1);

	end
	
	function [v_e_equil,v_i_equil]=liley_calculate_equilibrium(p,guess_ve)

		if nargin==1
			guess_ve=p.h_e_r;
		end

		r_abs=0;  

		if isfield(p,'gamma_ee_0')
			delta_ee=1/p.gamma_ee_0;
			delta_ei=1/p.gamma_ei_0;
			delta_ie=1/p.gamma_ie_0;
			delta_ii=1/p.gamma_ii_0;
		else
			delta_ee=1/p.gamma_ee;
			delta_ei=1/p.gamma_ei;
			delta_ie=1/p.gamma_ie;
			delta_ii=1/p.gamma_ii;    
		end

		y_ee=p.h_ee_eq/p.h_e_r-1;
		y_ei=p.h_ei_eq/p.h_i_r-1;
		y_ie=p.h_ie_eq/p.h_e_r-1;
		y_ii=p.h_ii_eq/p.h_i_r-1;
		u_e=p.mu_e/p.h_e_r-1;
		u_i=p.mu_i/p.h_i_r-1;
		r_e=1-r_abs*p.S_e_max;
		r_i=1-r_abs*p.S_i_max;
		v_e=p.sigma_e/(sqrt(2)*p.h_e_r);
		v_i=p.sigma_i/(sqrt(2)*p.h_i_r);

		s_ee=exp(p.gamma_ee*delta_ee)*p.Gamma_ee/(p.gamma_ee*p.h_e_r)*(p.N_ee_a+p.N_ee_b)*p.S_e_max;
		s_ei=exp(p.gamma_ei*delta_ei)*p.Gamma_ei/(p.gamma_ei*p.h_i_r)*(p.N_ei_a+p.N_ei_b)*p.S_e_max;
		s_ie=exp(p.gamma_ie*delta_ie)*p.Gamma_ie/(p.gamma_ie*p.h_e_r)*p.N_ie_b*p.S_i_max;
		s_ii=exp(p.gamma_ii*delta_ii)*p.Gamma_ii/(p.gamma_ii*p.h_i_r)*p.N_ii_b*p.S_i_max;

		p_ee=exp(p.gamma_ee*delta_ee)*p.Gamma_ee/(p.gamma_ee*p.h_e_r)*p.p_ee;
		p_ei=exp(p.gamma_ei*delta_ei)*p.Gamma_ei/(p.gamma_ei*p.h_i_r)*p.p_ei;
		p_ie=exp(p.gamma_ie*delta_ie)*p.Gamma_ie/(p.gamma_ie*p.h_e_r)*p.p_ie;
		p_ii=exp(p.gamma_ii*delta_ii)*p.Gamma_ii/(p.gamma_ii*p.h_i_r)*p.p_ii;

		syms x_e
		syms x_i

		x_i=u_i-v_i*log(1/r_i*(s_ie*(-1/(p_ie+abs(y_ie)/(y_ie-x_e)*(x_e+(y_ee-x_e)/abs(y_ee)*(s_ee/(1+r_e*exp(-(x_e-u_e)/v_e))+p_ee))))-1));
		s=vpasolve(x_i+(y_ei-x_i)/abs(y_ei)*(s_ei/(1+r_e*exp(-(x_e-u_e)/v_e))+p_ei)+(y_ii-x_i)/abs(y_ii)*(s_ii/(1+r_i*exp(-(x_i-u_i)/v_i))+p_ii)==0,x_e,-70/guess_ve-1);

		if isempty(s)
            s=vpasolve(x_i+(y_ei-x_i)/abs(y_ei)*(s_ei/(1+r_e*exp(-(x_e-u_e)/v_e))+p_ei)+(y_ii-x_i)/abs(y_ii)*(s_ii/(1+r_i*exp(-(x_i-u_i)/v_i))+p_ii)==0,x_e,[-90/p.h_e_r-1 0/p.h_e_r-1] );
		end

		s=double(s);
		clear x_i x_e

        try
			x_e=s;
			x_i=u_i-v_i*log(1/r_i*(s_ie*(-1/(p_ie+abs(y_ie)/(y_ie-x_e)*(x_e+(y_ee-x_e)/abs(y_ee)*(s_ee/(1+r_e*exp(-(x_e-u_e)/v_e))+p_ee))))-1));
			v_e_equil=(1+x_e)*p.h_e_r;
			v_i_equil=(1+x_i)*p.h_i_r;
		catch
			disp('Failed to calculate equilibrium potentials, using resting membrane potentials') 
			v_e_equil=p.h_e_r;
			v_i_equil=p.h_i_r;
        end

    end



end
