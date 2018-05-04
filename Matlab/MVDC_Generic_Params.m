 %--------------------------------------------------------------------------
% MVDC_Generic_Params.m
% Author: Nathan Ainsworth, DRS Naval Power Systems
% Description: Matlab script that defines all of the parameters for the
% MVDC fundamental subsystem model in the Matlab workspace. 

Ts=12.5e-6;

% Define the script configuration
verbose = 0;
fignum = 1;

% Test Toolboxes
Control_Systems_Toolbox_Installed = license('test', 'Control_Toolbox');
if(Control_Systems_Toolbox_Installed)
    s = tf('s');
end

% Simulation Paramters
f_DSP_Clock = 100e6;
T_DSP_Clock = 10*1/f_DSP_Clock;
t_enable = 50e-3;
t_clear = 1e-3;

% Variant Control - Scalar Variable
SWITCHING = 0;

%--------------------------------------------------------------------------
% System Rated Values

% MVDC Bus
SystemRatings.MVDC.P_Rated = 10e6;
SystemRatings.MVDC.Vdc_Rated = 12000;
SystemRatings.MVDC.Idc_Rated = SystemRatings.MVDC.P_Rated/SystemRatings.MVDC.Vdc_Rated;
SystemRatings.MVDC.Zdc_Rated = SystemRatings.MVDC.Vdc_Rated/SystemRatings.MVDC.Idc_Rated;

% MVAC Bus
% SystemRatings.MVAC.Vac_Rated = 6.9e3/sqrt(3);
SystemRatings.MVAC.P_Rated = 10e6;
SystemRatings.MVAC.Vac_Rated =  2000/sqrt(3); %4160/sqrt(3);
SystemRatings.MVAC.Iac_Rated = SystemRatings.MVAC.P_Rated/SystemRatings.MVAC.Vac_Rated/3;
SystemRatings.MVAC.Zac_Rated = SystemRatings.MVAC.Vac_Rated/SystemRatings.MVAC.Iac_Rated;
SystemRatings.MVAC.f_Rated = 60; %was 400 before for some reason.

% LVDC Bus
SystemRatings.LVDC.P_Rated = 3e6;
SystemRatings.LVDC.Vdc_Rated = 1000;
SystemRatings.LVDC.Idc_Rated = SystemRatings.LVDC.P_Rated/SystemRatings.LVDC.Vdc_Rated;
SystemRatings.LVDC.Zdc_Rated = SystemRatings.LVDC.Vdc_Rated/SystemRatings.LVDC.Idc_Rated;

% LVAC Bus
SystemRatings.LVAC.P_Rated = 12e3;
SystemRatings.LVAC.Vac_Rated = 450/sqrt(3);
SystemRatings.LVAC.Iac_Rated = SystemRatings.LVAC.P_Rated/SystemRatings.LVAC.Vac_Rated/3;
SystemRatings.LVAC.Zac_Rated = SystemRatings.LVAC.Vac_Rated/SystemRatings.LVAC.Iac_Rated;
SystemRatings.LVAC.f_Rated = 60;

%--------------------------------------------------------------------------
% PGM Synchronous Machine Parameters
PGM.SM.S_Rated = SystemRatings.MVAC.P_Rated;
PGM.SM.Iac_Rated = SystemRatings.MVAC.Iac_Rated; 
PGM.SM.Vac_Rated = SystemRatings.MVAC.Vac_Rated; 
PGM.SM.Ifield_Rated = 1000;
PGM.SM.N_PolePairs = 4;


PGM.SM.fe_Rated = SystemRatings.MVAC.f_Rated;
PGM.SM.we_Rated = 2*pi*PGM.SM.fe_Rated;
PGM.SM.fm_Rated = PGM.SM.fe_Rated/PGM.SM.N_PolePairs;
PGM.SM.wm_Rated = 2*pi*PGM.SM.fm_Rated;
PGM.SM.Tm_Rated = PGM.SM.S_Rated/PGM.SM.wm_Rated;
PGM.SM.Tm_Max = 2*PGM.SM.Tm_Rated; 
PGM.SM.rpm_Rated = PGM.SM.fm_Rated*60;

PGM.SM.Xd_Synchronous = 1.81; 
PGM.SM.Xq_Synchronous = 1.76;
PGM.SM.Xd_Transient = 0.3; 
PGM.SM.Xd_Subtransient = 0.23; 
PGM.SM.Xq_Subtransient = 0.25;
PGM.SM.Td0_Transient = 1;
PGM.SM.Td0_Subtransient = 0.07;
PGM.SM.Tq0_Subtransient = 0.03; 

PGM.SM.H = 3; % sec
PGM.SM.D = 0.01; % PU
PGM.SM.E_Rated = PGM.SM.S_Rated*PGM.SM.H;
PGM.SM.I = 2*PGM.SM.E_Rated/PGM.SM.we_Rated^2;
PGM.SM.F = PGM.SM.Tm_Rated*0.01/PGM.SM.we_Rated;
PGM.SM.Tf = 0;

PGM.SM.Ls = 835e-6;
PGM.SM.Ll = 0.1*PGM.SM.Ls;
PGM.SM.Lm = PGM.SM.Ls - PGM.SM.Ll;
Saliency = 1.76/1.81;
PGM.SM.Lmd = PGM.SM.Ls;
PGM.SM.Lmq = Saliency*PGM.SM.Ls;
PGM.SM.Rs = 0.18; 
PGM.SM.Vm_Const = PGM.SM.Vac_Rated*sqrt(3)*sqrt(2)/(PGM.SM.rpm_Rated/1000);

%--------------------------------------------------------------------------
% PGM Rectifier Paramters, Nathan's version
PGM.Rectifier.P_Rated = SystemRatings.MVAC.P_Rated;
PGM.Rectifier.Vac_Rated = PGM.SM.Vac_Rated;
PGM.Rectifier.Vdc_Rated = SystemRatings.MVDC.Vdc_Rated;
PGM.Rectifier.Vdc_Natural = 3*sqrt(3)*sqrt(2)*PGM.Rectifier.Vac_Rated/pi;

PGM.Rectifier.Vf = 0.8;
PGM.Rectifier.Ron = 0.01;
PGM.Rectifier.Roff = 100e6;
PGM.Rectifier.Rsnubber = 500;
PGM.Rectifier.Csnubber = 250e-6;
PGM.Rectifier.Rshunt = 100e6;

PGM.Rectifier.Lfilter = 100e-6;
PGM.Rectifier.Rfilter = 5e-6;
PGM.Rectifier.Cfilter = 1e-3;
PGM.Rectifier.ESRfilter = 10e-3;
PGM.Rectifier.wc = 1/sqrt(PGM.Rectifier.Lfilter*PGM.Rectifier.Cfilter);
PGM.Rectifier.Vdc_Init = PGM.Rectifier.Vdc_Rated;

%from Willy's model
% PGM.Rectifier.feRated = 60;  %PGM
% PGM.Rectifier.VdcRated = 4000;
% PGM.Rectifier.Vcmd = 4000;
% 
% PGM.Rectifier.VARated = 10e6;
% PGM.Rectifier.VLLRated = 2000;
% 
% PGM.Rectifier.Sb = PGM.Rectifier.VARated; 
% PGM.Rectifier.Vb =PGM.Rectifier.VdcRated %sqrt(2/3)*VLLRated;
% PGM.Rectifier.Ib =PGM.Rectifier.Sb/PGM.Rectifier.Vb; %VARated/sqrt(3)/PGM.Rectifier.VLLRated*sqrt(2);
% PGM.Rectifier.Zb = PGM.Rectifier.Vb/PGM.Rectifier.Ib;
% PGM.Rectifier.Zfpu = 0.02;
% PGM.Rectifier.XoverR = 20;
% PGM.Rectifier.Zf = PGM.Rectifier.Zfpu*PGM.Rectifier.Zb/2;
% PGM.Rectifier.Rf = PGM.Rectifier.Zf/sqrt(1 + PGM.Rectifier.XoverR^2);
% PGM.Rectifier.Xf = sqrt(PGM.Rectifier.Zf^2 - PGM.Rectifier.Rf^2);
% PGM.Rectifier.Lf = PGM.Rectifier.Xf/2/pi/PGM.Rectifier.feRated;
% PGM.Rectifier.Rinrush = 10;
% 
% PGM.Rectifier.Prated = PGM.Rectifier.VARated;
% 
% %%%
% PGM.Rectifier.Droop = 3;
% PGM.Rectifier.KDroop_VSR = PGM.Rectifier.Droop/100*320/PGM.Rectifier.Prated;
% PGM.Rectifier.fsw = 20000; %39960;
% PGM.Rectifier.Ldm1 = 4*62.5e-6;
% PGM.Rectifier.Ldm2 = 20*62.5e-6;
% PGM.Rectifier.Cdm = 10e-6;
% PGM.Rectifier.Cd = 3*PGM.Rectifier.Cdm;
% PGM.Rectifier.Rd = sqrt(PGM.Rectifier.Ldm1/3/PGM.Rectifier.Cdm);
% PGM.Rectifier.Lcm1 = 1e-3;
% PGM.Rectifier.Lcm2 = 1e-3;
% PGM.Rectifier.Ccm = 1e-6;

%end Willy's model


%--------------------------------------------------------------------------
% PGM Controller Parameters

PGM_Controller.T_s = 100e-6;
PGM_Controller.we_Min = PGM.SM.we_Rated *2/3;

PGM_Controller.Governor.wc = 2*pi*1*6;
PGM_Controller.Governor.kp = 1/0.05;
PGM_Controller.Governor.ki = PGM_Controller.Governor.kp*PGM_Controller.Governor.wc;

PGM_Controller.Vdc_Regulator.Kp = 1*pi/3/sqrt(3)/sqrt(2)*PGM.Rectifier.Vdc_Rated/PGM.Rectifier.Vac_Rated;
PGM_Controller.Vdc_Regulator.Ki = PGM_Controller.Vdc_Regulator.Kp*2*pi/PGM.SM.Td0_Transient;

%--------------------------------------------------------------------------
% PEBB 6000 Parameters
PEBB6000.Cdc = 1e-3;
PEBB6000.ESRdc = 10e-3;
PEBB6000.Ron = 1e-3;
PEBB6000.Roff = 100e6;
PEBB6000.Vth = 3;
PEBB6000.Vf = 1.5;
PEBB6000.Rsnubber = 60e6;
PEBB6000.Csnubber = inf;
PEBB6000.Vdc_Rated = 6e3;
PEBB6000.Vdc_Init = PEBB6000.Vdc_Rated;

% PEBB 1000 Parameters
PEBB1000.Cdc = 1e-3;
PEBB1000.ESRdc = 10e-3;
PEBB1000.Ron = 1e-3;
PEBB1000.Roff = 100e6;
PEBB1000.Vth = 2;
PEBB1000.Vf = 1.2;
PEBB1000.Rsnubber = 60e6;
PEBB1000.Csnubber = inf;
PEBB1000.Vdc_Rated = 1e3;
PEBB1000.Vdc_Init = PEBB1000.Vdc_Rated;

%--------------------------------------------------------------------------
% PEBB Control Parameters
PEBB6000Controller.f_PWM = 10e3;
PEBB6000Controller.T_PWM = 1/PEBB6000Controller.f_PWM;
PEBB6000Controller.T_PWM_Update = T_DSP_Clock;
PEBB6000Controller.t_DT = 2e-6;
PEBB6000Controller.DutyCycle_Max = 1 - PEBB6000Controller.t_DT/PEBB6000Controller.T_PWM;
PEBB6000Controller.DutyCycle_Min = PEBB6000Controller.t_DT/PEBB6000Controller.T_PWM;

PEBB1000Controller.f_PWM = 10e3;
PEBB1000Controller.T_PWM = 1/PEBB1000Controller.f_PWM;
PEBB1000Controller.T_PWM_Update = T_DSP_Clock;
PEBB1000Controller.t_DT = 2e-6;
PEBB1000Controller.DutyCycle_Max = 1 - PEBB1000Controller.t_DT/PEBB1000Controller.T_PWM;
PEBB1000Controller.DutyCycle_Min = PEBB1000Controller.t_DT/PEBB1000Controller.T_PWM;

%--------------------------------------------------------------------------
% JFET Parameters
JFET_UJN1205K.Vbi = 1;
JFET_UJN1205K.Vth = -8;
JFET_UJN1205K.Go = 30;
JFET_UJN1205K.Cgs = 1020e-12;
JFET_UJN1205K.Rgs = 30e-9/JFET_UJN1205K.Cgs/3;
JFET_UJN1205K.Cds = 154e-12;
JFET_UJN1205K.Rds = 0.045;

%--------------------------------------------------------------------------
% MMC Parameters
MMC.Vdc_Rated = SystemRatings.MVDC.Vdc_Rated;
MMC.Vac_Rated = SystemRatings.MVAC.Vac_Rated;
MMC.S_Rated = SystemRatings.MVAC.P_Rated;
MMC.N_Modules_per_Arm = 2;

% MMC.Larm = 300e-6;
MMC.Larm = 100e-6;
MMC.Rarm = 10e-3;
% MMC.Lfilter = 600e-6;
MMC.Lfilter = 200e-6;
MMC.Rfilter = 20e-3;
MMC.Lac_Eff = MMC.Lfilter+MMC.Larm/2;
MMC.Rac_Eff = MMC.Rfilter+MMC.Rarm/2;
% MMC.Cfilter = 33e-6;
MMC.Cfilter = 100e-6;
% MMC.ESRfilter = 400e-3;
MMC.ESRfilter = 150e-3;
if(Control_Systems_Toolbox_Installed)
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        bode_pretty((1 + MMC.ESRfilter*MMC.Cfilter*s)/(MMC.Lac_Eff*MMC.Cfilter*s^2 + (MMC.Rac_Eff + MMC.ESRfilter)*MMC.Cfilter*s + 1), 1, 10e3)
        title('MMC AC Filter Open-Loop Transfer Function');
    end
end

MMC.Cdc = 1e-3;
MMC.ESRdc = 5e-3;

MMC.Vdc_Init = MMC.Vdc_Rated;
MMC.T_Breaker_Open = 1e-3;

%--------------------------------------------------------------------------
% MMC Controller Parameters
MMC_Controller.T_s = PEBB6000Controller.T_PWM;
MMC_Controller.f_PWM = PEBB6000Controller.f_PWM;

if(Control_Systems_Toolbox_Installed)
    z = tf('z', MMC_Controller.T_s);
end

% Idq0 Regulator Parameters
% Design the Idq0 regulator so that Ki compensated the zero introduced by
% Rac, and Kp places the gain-crossover frequency at 1000 Hz. 
MMC_Controller.Idq_Regulator.wc = 2*pi*1000;
MMC_Controller.Idq_Regulator.wz = 10*MMC.Rac_Eff/MMC.Lac_Eff;
MMC_Controller.Idq_Regulator.Kp = MMC.Lac_Eff*MMC_Controller.Idq_Regulator.wc;
MMC_Controller.Idq_Regulator.Ki = MMC_Controller.Idq_Regulator.Kp*MMC_Controller.Idq_Regulator.wz;
if(Control_Systems_Toolbox_Installed)
    MMC_Controller.Idq_Regulator.G_OpenLoop_TF = (MMC_Controller.Idq_Regulator.Kp+MMC_Controller.Idq_Regulator.Ki/s)*1/(MMC.Lac_Eff*s+MMC.Rac_Eff)*exp(-MMC_Controller.T_s*s);
    MMC_Controller.Idq_Regulator.G_ClosedLoop_TF = feedback((MMC_Controller.Idq_Regulator.Kp+MMC_Controller.Idq_Regulator.Ki/s)*1/(MMC.Lac_Eff*s+MMC.Rac_Eff), exp(-MMC_Controller.T_s*s));
    [MMC_Controller.Idq_Regulator.Gm, MMC_Controller.Idq_Regulator.Pm, MMC_Controller.Idq_Regulator.wgm, MMC_Controller.Idq_Regulator.wpm] = margin(MMC_Controller.Idq_Regulator.G_OpenLoop_TF);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin(MMC_Controller.Idq_Regulator.G_OpenLoop_TF);
        disp(['MMC Idq Regulator: Gm = ' num2str(MMC_Controller.Idq_Regulator.Gm) ' @ ' num2str(MMC_Controller.Idq_Regulator.wgm/2/pi) ' Hz, Pm = ' num2str(MMC_Controller.Idq_Regulator.Pm) ' deg @ ' num2str(MMC_Controller.Idq_Regulator.wpm/2/pi) ' Hz']);
        title('MMC DQ0 Current Regulator Open-Loop Transfer Function');
    end
end

% Vdq0 Regulator Parameters
% Design the Vdq regulator so that Kp places the gain-crossover frequency
% at 100 Hz. 
MMC_Controller.Vdq_Regulator.wc = 2*pi*500;
MMC_Controller.Vdq_Regulator.wz = 2*pi*20;
MMC_Controller.Vdq_Regulator.Kp = MMC.Cfilter*MMC_Controller.Vdq_Regulator.wc;
MMC_Controller.Vdq_Regulator.Ki = MMC_Controller.Vdq_Regulator.Kp*MMC_Controller.Vdq_Regulator.wz;
if(Control_Systems_Toolbox_Installed)
    MMC_Controller.Vdq_Regulator.G_OpenLoop_TF = (MMC_Controller.Vdq_Regulator.Kp+MMC_Controller.Vdq_Regulator.Ki/s)*MMC_Controller.Idq_Regulator.G_ClosedLoop_TF*(1/s/MMC.Cfilter + MMC.ESRfilter)*exp(-MMC_Controller.T_s*s);
    [MMC_Controller.Vdq_Regulator.Gm, MMC_Controller.Vdq_Regulator.Pm, MMC_Controller.Vdq_Regulator.wgm, MMC_Controller.Vdq_Regulator.wpm] = margin(MMC_Controller.Vdq_Regulator.G_OpenLoop_TF);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin(MMC_Controller.Vdq_Regulator.G_OpenLoop_TF);
        disp(['MMC Vdq Regulator: Gm = ' num2str(MMC_Controller.Vdq_Regulator.Gm) ' @ ' num2str(MMC_Controller.Vdq_Regulator.wgm/2/pi) ' Hz, Pm = ' num2str(MMC_Controller.Vdq_Regulator.Pm) ' deg @ ' num2str(MMC_Controller.Vdq_Regulator.wpm/2/pi) ' Hz']);
        title('MMC DQ0 Voltage Regulator Open-Loop Transfer Function');
    end
end

% Vdc Regulator Parameters
MMC_Controller.Vdc_Regulator.wc = 2*pi*500;
MMC_Controller.Vdc_Regulator.wz = 0; %2*pi*100;
% MMC_Controller.Vdc_Regulator.Kp = 1;
MMC_Controller.Vdc_Regulator.Kp = MMC.Cdc*MMC_Controller.Vdc_Regulator.wc^2/sqrt(MMC_Controller.Vdc_Regulator.wc^2+MMC_Controller.Vdc_Regulator.wz^2);
MMC_Controller.Vdc_Regulator.Ki = MMC_Controller.Vdc_Regulator.Kp*MMC_Controller.Vdc_Regulator.wz;
if(Control_Systems_Toolbox_Installed)
    MMC_Controller.Vdc_Regulator.G_OpenLoop_TF = (MMC_Controller.Vdc_Regulator.Kp+MMC_Controller.Vdc_Regulator.Ki/s)*(1/(MMC.Cdc*s)+MMC.ESRdc)*MMC_Controller.Idq_Regulator.G_ClosedLoop_TF*exp(-MMC_Controller.T_s*s);
    MMC_Controller.Vdc_Regulator.G_ClosedLoop_TF = feedback((MMC_Controller.Vdc_Regulator.Kp+MMC_Controller.Vdc_Regulator.Ki/s)*(1/(MMC.Cdc*s)+MMC.ESRdc)*MMC_Controller.Idq_Regulator.G_ClosedLoop_TF, exp(-MMC_Controller.T_s*s));
    [MMC_Controller.Vdc_Regulator.Gm, MMC_Controller.Vdc_Regulator.Pm, MMC_Controller.Vdc_Regulator.wgm, MMC_Controller.Vdc_Regulator.wpm] = margin(MMC_Controller.Vdc_Regulator.G_OpenLoop_TF);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin(MMC_Controller.Vdc_Regulator.G_OpenLoop_TF);
        disp(['MMC Vdc Regulator: Gm = ' num2str(MMC_Controller.Vdc_Regulator.Gm) ' @ ' num2str(MMC_Controller.Vdc_Regulator.wgm/2/pi) ' Hz, Pm = ' num2str(MMC_Controller.Vdc_Regulator.Pm) ' deg @ ' num2str(MMC_Controller.Vdc_Regulator.wpm/2/pi) ' Hz']);
        title('MMC DC Voltage Regulator Open-Loop Transfer Function');
    end
end

% PLL Parameters
MMC_Controller.Kp_PLL = 300;
MMC_Controller.Ki_PLL = MMC_Controller.Kp_PLL*2*pi*100;
MMC_Controller.Vq_PLL_lock_threshold = 0.02;

% Ileg_DC Regulator Parameters
MMC_Controller.omega_gc_Ileg = 2*pi*100*6;
% MMC_Controller.Kp_Ileg_DC = 1;
MMC_Controller.Kp_Ileg_DC = MMC_Controller.omega_gc_Ileg*MMC.Larm; 
MMC_Controller.Ki_Ileg_DC = MMC_Controller.Kp_Ileg_DC*MMC.Rarm/MMC.Larm;
if(Control_Systems_Toolbox_Installed)
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin((MMC_Controller.Kp_Ileg_DC+MMC_Controller.Ki_Ileg_DC/s)*1/(s*MMC.Larm+MMC.Rarm)*exp(-PEBB6000Controller.T_PWM*s))
        title('MMC Ileg DC Regulator Open-Loop Bode Plot');
    end
end

% Circulating Current Elimination Repetitive Controller
Circulating_Current_Repetitive.f_Rated = SystemRatings.MVAC.f_Rated;
Circulating_Current_Repetitive.T_Rated = 1/Circulating_Current_Repetitive.f_Rated;
Circulating_Current_Repetitive.T_s = MMC_Controller.T_s;
Circulating_Current_Repetitive.fs = 1/Circulating_Current_Repetitive.T_s;
Circulating_Current_Repetitive.Ns = floor(Circulating_Current_Repetitive.T_Rated/Circulating_Current_Repetitive.T_s);
% if(Circulating_Current_Repetitive.Ns ~= Circulating_Current_Repetitive.T_Rated/Circulating_Current_Repetitive.T_s)
%     error(['MVDC_Generic_Params: Circulating current repetitive control requires that T_Rated/T_s is an integer']);
% end    %I commented this out because I checked that T_Rated/T_s is an
% integer and it is, but Simulink is still generating an error. -Anna
% Brinck
Circulating_Current_Repetitive.L = 1;
Circulating_Current_Repetitive.Kr = 0.5;

Circulating_Current_Repetitive.Moving_Average_Filter.Kp = 0.99;
Circulating_Current_Repetitive.Moving_Average_Filter.N = 3;
Circulating_Current_Repetitive.Moving_Average_Filter.Ns_Lead = (Circulating_Current_Repetitive.Moving_Average_Filter.N-1)/2;
Circulating_Current_Repetitive.Moving_Average_Filter.Num = Circulating_Current_Repetitive.Moving_Average_Filter.Kp/Circulating_Current_Repetitive.Moving_Average_Filter.N*ones(1, Circulating_Current_Repetitive.Moving_Average_Filter.N);
Circulating_Current_Repetitive.Moving_Average_Filter.Den = [1 zeros(1, Circulating_Current_Repetitive.Moving_Average_Filter.N-1)];
if(Control_Systems_Toolbox_Installed)
    Circulating_Current_Repetitive.Moving_Average_Filter.DTF = tf(Circulating_Current_Repetitive.Moving_Average_Filter.Num, Circulating_Current_Repetitive.Moving_Average_Filter.Den, Circulating_Current_Repetitive.T_s);
    % bode(Circulating_Current_Repetitive.Moving_Average_Filter.DTF);
end

Circulating_Current_Repetitive.N_Fwd_Delay = Circulating_Current_Repetitive.Ns - Circulating_Current_Repetitive.Moving_Average_Filter.Ns_Lead - Circulating_Current_Repetitive.L;
Circulating_Current_Repetitive.N_Rev_Delay = Circulating_Current_Repetitive.L;
if(Control_Systems_Toolbox_Installed)
    Circulating_Current_Repetitive.DTF = Circulating_Current_Repetitive.Kr*z^(-Circulating_Current_Repetitive.Ns+Circulating_Current_Repetitive.L)/(1 - Circulating_Current_Repetitive.Moving_Average_Filter.DTF*z^(Circulating_Current_Repetitive.Ns+Circulating_Current_Repetitive.L));
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        bode_pretty(Circulating_Current_Repetitive.DTF, 1, Circulating_Current_Repetitive.fs/2);
        title('MMC Ileg Repetitive Controller Bode Plot');
    end
end

% Vc Average Regulator
MMC_Controller.omega_gc_Vc = 2*pi*10;
MMC_Controller.Kp_Vc = 0.1; %PEBB6000.Cdc*MMC_Controller.omega_gc_Vc/MMC_Controller.omega_gc_Ileg;
MMC_Controller.Ki_Vc = MMC_Controller.Kp_Vc*MMC_Controller.omega_gc_Ileg; 

% Vc Upper/Lower Balance Regulator
MMC_Controller.omega_gc_Vc = 2*pi*10;
MMC_Controller.Kp_Vc_UL_Balance = 0.1/SystemRatings.MVDC.Idc_Rated; %PEBB6000.Cdc*MMC_Controller.omega_gc_Vc/MMC_Controller.omega_gc_Ileg;
MMC_Controller.Ki_Vc_UL_Balance = MMC_Controller.Kp_Vc_UL_Balance*2*pi*10; 

%--------------------------------------------------------------------------
% DAB Parameters
DAB.N_Series_MVDC_PEBBs = 1;
DAB.N_Parallel_LVDC_PEBBs = 3;

% Define the ratings of the DAB
DAB.P_Rated = SystemRatings.LVDC.P_Rated/2;
DAB.V_MVDC_Rated = SystemRatings.MVDC.Vdc_Rated/2;
DAB.I_MVDC_Rated = DAB.P_Rated/DAB.V_MVDC_Rated;
DAB.Z_MVDC_Rated = DAB.V_MVDC_Rated/DAB.I_MVDC_Rated;
DAB.V_LVDC_Rated = SystemRatings.LVDC.Vdc_Rated;
DAB.I_LVDC_Rated = DAB.P_Rated/DAB.V_LVDC_Rated;
DAB.Z_LVDC_Rated = DAB.V_LVDC_Rated/DAB.I_LVDC_Rated;

DAB.f_Rated = 10e3;
DAB.w_Rated = DAB.f_Rated*2*pi;

% Define the DAB transformer characteristics
DAB.Pmax = DAB.P_Rated*2;
DAB.Vpri = DAB.V_MVDC_Rated*2*sqrt(2)/pi;
DAB.Vsec = DAB.V_LVDC_Rated*2*sqrt(2)/pi;
DAB.N_Transformer = DAB.V_MVDC_Rated/DAB.V_LVDC_Rated;
DAB.N_Turns_Primary = DAB.N_Transformer;
DAB.N_Turns_Secondary = 1; 
DAB.Lpri = DAB.V_MVDC_Rated^2/DAB.w_Rated/DAB.Pmax/2;
DAB.Lsec = DAB.V_LVDC_Rated^2/DAB.w_Rated/(DAB.Pmax/DAB.N_Parallel_LVDC_PEBBs)/2;
DAB.Rpri = 0.1*DAB.Lpri*DAB.w_Rated; 
DAB.Rsec = 0.1*DAB.Lsec*DAB.w_Rated; 
DAB.Pmag = DAB.Pmax*0.005;
DAB.Qmag = DAB.Pmax*0.02;
DAB.Rmag = DAB.Vpri^2/DAB.Pmag;
DAB.Lmag = DAB.Vpri^2/DAB.Qmag/DAB.w_Rated;

% Calculate the the equivalent impedances as seen on the primary and
% secondary sides of the DAB
DAB.Cdc_LVDC_Eff = DAB.N_Parallel_LVDC_PEBBs*PEBB1000.Cdc;
DAB.Cdc_MVDC_Eff = 2*PEBB6000.Cdc/DAB.N_Series_MVDC_PEBBs;
DAB.ESR_MVDC_Eff = PEBB6000.ESRdc*DAB.N_Series_MVDC_PEBBs; 
DAB.ESR_LVDC_Eff = PEBB1000.ESRdc/DAB.N_Parallel_LVDC_PEBBs; 

DAB.Leff_MVDC = DAB.Lpri+DAB.Lsec*DAB.N_Transformer^2;
DAB.Reff_MVDC = DAB.Rpri+DAB.Rsec*DAB.N_Transformer^2+2*PEBB6000.Ron+(2*PEBB1000.Ron*DAB.N_Transformer^2);
DAB.Zeff_MVDC = DAB.Reff_MVDC + j*DAB.w_Rated*DAB.Leff_MVDC; 
DAB.Yeff_MVDC = 1/DAB.Zeff_MVDC; 
DAB.Geff_MVDC = real(DAB.Yeff_MVDC);
DAB.Beff_MVDC = imag(DAB.Yeff_MVDC);

DAB.Leff_LVDC = DAB.Leff_MVDC/DAB.N_Transformer^2;
DAB.Reff_LVDC = DAB.Reff_MVDC/DAB.N_Transformer^2;
DAB.Zeff_LVDC = DAB.Reff_LVDC + j*DAB.w_Rated*DAB.Leff_LVDC; 
DAB.Yeff_LVDC = 1/DAB.Zeff_LVDC; 
DAB.Geff_LVDC = real(DAB.Yeff_LVDC);
DAB.Beff_LVDC = imag(DAB.Yeff_LVDC);

DAB.Loss_Angle_Txfr = atan2(DAB.Reff_LVDC, DAB.w_Rated*DAB.Leff_LVDC);

%--------------------------------------------------------------------------
% DAB Controller Parameters

% Assign the basic timing parameters of the DAB controller
DAB_Controller.fs = DAB.f_Rated;
DAB_Controller.ws = 2*pi*DAB.f_Rated; 
DAB_Controller.T_s = 1/DAB_Controller.fs;
DAB_Controller.T_DSP_Clock = T_DSP_Clock;
DAB_Controller.N_Fullbridge_PWMClock_2_Delay = floor(PEBB1000Controller.T_PWM/PEBB1000Controller.T_PWM_Update/2);

% Define the controller current limits
DAB_Controller.Iin_Max = 1.25*DAB.I_MVDC_Rated;
DAB_Controller.Iin_Min = -1.25*DAB.I_MVDC_Rated; 
DAB_Controller.Io_Total_Max = 1.1*DAB.I_LVDC_Rated;
DAB_Controller.Io_Total_Min = -1.1*DAB.I_LVDC_Rated;

% Design a feedback Anti-Aliasing filter to capture the average values of
% the V and I DC signals
DAB_Controller.VI_LowPass_Filter.fc = 3000;
DAB_Controller.VI_LowPass_Filter.wc = 2*pi*DAB_Controller.VI_LowPass_Filter.fc;
DAB_Controller.VI_LowPass_Filter.N = 2;
% [DAB_Controller.VI_LowPass_Filter.Num, DAB_Controller.VI_LowPass_Filter.Den] = butter(DAB_Controller.VI_LowPass_Filter.N, DAB_Controller.VI_LowPass_Filter.fc*2*pi, 's');
DAB_Controller.VI_LowPass_Filter.Num = DAB_Controller.VI_LowPass_Filter.wc^2;
DAB_Controller.VI_LowPass_Filter.Den = [1 sqrt(2)*DAB_Controller.VI_LowPass_Filter.wc DAB_Controller.VI_LowPass_Filter.wc^2];
if(Control_Systems_Toolbox_Installed)
    DAB_Controller.VI_LowPass_Filter.TF = tf(DAB_Controller.VI_LowPass_Filter.Num, DAB_Controller.VI_LowPass_Filter.Den);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        bode_pretty(DAB_Controller.VI_LowPass_Filter.TF, 1, DAB_Controller.fs);
        title('DAB LowPass Filter Bode Plot');
    end
end

% Design a feedback Notch filter to eliminate the switching frequency from
% DC measurements
DAB_Controller.VI_DC_Notch_Filter.wo = DAB_Controller.ws;
DAB_Controller.VI_DC_Notch_Filter.Lf = 10e-6;
DAB_Controller.VI_DC_Notch_Filter.Cf = 1/DAB_Controller.VI_DC_Notch_Filter.Lf/DAB_Controller.VI_DC_Notch_Filter.wo^2;
DAB_Controller.VI_DC_Notch_Filter.gamma = 1;
DAB_Controller.VI_DC_Notch_Filter.Rf = DAB_Controller.VI_DC_Notch_Filter.gamma/DAB_Controller.VI_DC_Notch_Filter.Cf/DAB_Controller.VI_DC_Notch_Filter.wo; 
DAB_Controller.VI_DC_Notch_Filter.Num = [DAB_Controller.VI_DC_Notch_Filter.Lf*DAB_Controller.VI_DC_Notch_Filter.Cf 0 1];
DAB_Controller.VI_DC_Notch_Filter.Den = [DAB_Controller.VI_DC_Notch_Filter.Lf*DAB_Controller.VI_DC_Notch_Filter.Cf DAB_Controller.VI_DC_Notch_Filter.Rf*DAB_Controller.VI_DC_Notch_Filter.Cf 1];
if(Control_Systems_Toolbox_Installed)
%     DAB_Controller.VI_DC_Notch_Filter.TF = (DAB_Controller.VI_DC_Notch_Filter.Lf*DAB_Controller.VI_DC_Notch_Filter.Cf*s^2 + 1)/(DAB_Controller.VI_DC_Notch_Filter.Lf*DAB_Controller.VI_DC_Notch_Filter.Cf*s^2 + DAB_Controller.VI_DC_Notch_Filter.Rf*DAB_Controller.VI_DC_Notch_Filter.Cf*s + 1);
%     [DAB_Controller.VI_DC_Notch_Filter.Num, DAB_Controller.VI_DC_Notch_Filter.Den] = tfdata(DAB_Controller.VI_DC_Notch_Filter.TF, 'v');
    DAB_Controller.VI_DC_Notch_Filter.TF = tf(DAB_Controller.VI_DC_Notch_Filter.Num, DAB_Controller.VI_DC_Notch_Filter.Den);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        bode_pretty(DAB_Controller.VI_DC_Notch_Filter.TF, 1, DAB_Controller.fs);
        title('DAB Notch Filter Bode Plot');
    end
end

% Design a regulator to control the output current Io
DAB_Controller.SPS_Regulator_Idc.Po_max = 8/pi^2*DAB.V_MVDC_Rated/DAB.N_Transformer*DAB.V_LVDC_Rated*abs(DAB.Yeff_LVDC); %/sqrt(DAB.Reff_LVDC^2+DAB_Controller.ws^2*DAB.Leff_LVDC^2);
DAB_Controller.SPS_Regulator_Idc.Io_Max = 8/pi^2*DAB.V_MVDC_Rated/DAB.N_Transformer*abs(DAB.Yeff_LVDC);
DAB_Controller.SPS_Regulator_Idc.Io_Min = -DAB_Controller.SPS_Regulator_Idc.Io_Max; 
DAB_Controller.SPS_Regulator_Idc.Io_Offset = 8/pi^2*DAB.Geff_LVDC*DAB.V_LVDC_Rated;

DAB_Controller.SPS_Regulator_Idc.G_FF_Vo =  8/pi^2*DAB.Geff_LVDC*DAB.N_Transformer/DAB.V_MVDC_Rated/abs(DAB.Yeff_LVDC); 
DAB_Controller.SPS_Regulator_Idc.G_FF_Io_Ref =  0; 
DAB_Controller.SPS_Regulator_Idc.Kp = 0.1/DAB_Controller.SPS_Regulator_Idc.Io_Max;
DAB_Controller.SPS_Regulator_Idc.Ki = DAB_Controller.SPS_Regulator_Idc.Kp*2*pi*DAB_Controller.VI_LowPass_Filter.fc; 
if(Control_Systems_Toolbox_Installed)
    DAB_Controller.SPS_Regulator_Idc.G_OpenLoop_TF = (DAB_Controller.SPS_Regulator_Idc.Kp + DAB_Controller.SPS_Regulator_Idc.Ki/s)*DAB_Controller.SPS_Regulator_Idc.Io_Max*DAB_Controller.VI_LowPass_Filter.TF*DAB_Controller.VI_DC_Notch_Filter.TF*exp(-DAB_Controller.T_s*s);
    DAB_Controller.SPS_Regulator_Idc.G_ClosedLoop_TF = feedback((DAB_Controller.SPS_Regulator_Idc.Kp + DAB_Controller.SPS_Regulator_Idc.Ki/s)*DAB_Controller.SPS_Regulator_Idc.Io_Max, DAB_Controller.VI_LowPass_Filter.TF*DAB_Controller.VI_DC_Notch_Filter.TF*exp(-DAB_Controller.T_s*s));
    [DAB_Controller.SPS_Regulator_Idc.Gm, DAB_Controller.SPS_Regulator_Idc.Pm, DAB_Controller.SPS_Regulator_Idc.wgm, DAB_Controller.SPS_Regulator_Idc.wpm] = margin(DAB_Controller.SPS_Regulator_Idc.G_OpenLoop_TF);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin(DAB_Controller.SPS_Regulator_Idc.G_OpenLoop_TF);
        disp(['DAB Idc Regulator: Gm = ' num2str(DAB_Controller.SPS_Regulator_Idc.Gm) ' @ ' num2str(DAB_Controller.SPS_Regulator_Idc.wgm/2/pi) ' Hz, Pm = ' num2str(DAB_Controller.SPS_Regulator_Idc.Pm) ' deg @ ' num2str(DAB_Controller.SPS_Regulator_Idc.wpm) ' Hz']);
        title('DAB Idc Regulator Open Loop Bode Plot');
    end
end

% Design a regulator to control the output DC voltage
DAB_Controller.V_LVDC_Regulator.omega_c = 2*pi*200;
% DAB_Controller.V_LVDC_Regulator.omega_z = 2*pi*50;
DAB_Controller.V_LVDC_Regulator.Kp = DAB.Cdc_LVDC_Eff*DAB_Controller.V_LVDC_Regulator.omega_c;
DAB_Controller.V_LVDC_Regulator.Ki = 0; %DAB_Controller.V_LVDC_Regulator.Kp*DAB_Controller.V_LVDC_Regulator.omega_c;
if(Control_Systems_Toolbox_Installed)
    DAB_Controller.V_LVDC_Regulator.G_OpenLoop_TF = DAB_Controller.V_LVDC_Regulator.Kp*(1/s/DAB.Cdc_LVDC_Eff + DAB.ESR_LVDC_Eff)*DAB_Controller.SPS_Regulator_Idc.G_ClosedLoop_TF*DAB_Controller.VI_LowPass_Filter.TF*DAB_Controller.VI_DC_Notch_Filter.TF*exp(-DAB_Controller.T_s*s);
    [DAB_Controller.V_LVDC_Regulator.Gm, DAB_Controller.V_LVDC_Regulator.Pm, DAB_Controller.V_LVDC_Regulator.wgm, DAB_Controller.V_LVDC_Regulator.wpm] = margin(DAB_Controller.V_LVDC_Regulator.G_OpenLoop_TF);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin(DAB_Controller.V_LVDC_Regulator.G_OpenLoop_TF);
        disp(['DAB Vdc Regulator: Gm = ' num2str(DAB_Controller.V_LVDC_Regulator.Gm) ' @ ' num2str(DAB_Controller.V_LVDC_Regulator.wgm/2/pi) ' Hz, Pm = ' num2str(DAB_Controller.V_LVDC_Regulator.Pm) ' deg @ ' num2str(DAB_Controller.V_LVDC_Regulator.wpm/2/pi) ' Hz']);
        title('DAB Vdc Regulator Open Loop Bode Plot');
    end
end

%--------------------------------------------------------------------------
% IPNC Parameters
IPNC.Output_Module1.Rfilter = 10e-3;
IPNC.Output_Module1.Lfilter = 500e-6;
IPNC.Output_Module1.Cfilter = 50e-6;
IPNC.Output_Module1.ESRfilter = 1/(2*pi*10e3*IPNC.Output_Module1.Cfilter); %Set Wz >= 10kHz
if(Control_Systems_Toolbox_Installed)
    IPNC.Output_Module1.Filter_TF = (1 + s*IPNC.Output_Module1.Cfilter*IPNC.Output_Module1.ESRfilter)/(IPNC.Output_Module1.Lfilter*IPNC.Output_Module1.Cfilter*s^2+(IPNC.Output_Module1.Rfilter + IPNC.Output_Module1.ESRfilter)*IPNC.Output_Module1.Cfilter*s + 1);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        bode_pretty(IPNC.Output_Module1.Filter_TF, 10, 10e3);
        title('IPNC Output Filter Bode Plot');
    end
end

%--------------------------------------------------------------------------
% IPNC Control Parameters
IPNC_Controller.T_s = 100e-6;
IPNC_Controller.VI_LowPass_Filter = DAB_Controller.VI_LowPass_Filter; 
IPNC_Controller.VI_Notch_Filter = DAB_Controller.VI_DC_Notch_Filter;
IPNC_Controller.Vdq_Regulator.Kp = 1e-3;
IPNC_Controller.Vdq_Regulator.Ki = 0;
IPNC_Controller.Idq_Regulator.wc = 2*pi*300;
IPNC_Controller.Idq_Regulator.wz = (IPNC.Output_Module1.Rfilter+PEBB1000.Ron)/IPNC.Output_Module1.Lfilter;
IPNC_Controller.Idq_Regulator.Kp = IPNC.Output_Module1.Lfilter*IPNC_Controller.Idq_Regulator.wc;
IPNC_Controller.Idq_Regulator.Ki = IPNC_Controller.Idq_Regulator.Kp*IPNC_Controller.Idq_Regulator.wz;
if(Control_Systems_Toolbox_Installed)
    IPNC_Controller.Idq_Regulator.G_OpenLoop_TF = (IPNC_Controller.Idq_Regulator.Kp+IPNC_Controller.Idq_Regulator.Ki/s)*1/(IPNC.Output_Module1.Lfilter*s+IPNC.Output_Module1.Rfilter+PEBB1000.Ron)*IPNC_Controller.VI_LowPass_Filter.TF*IPNC_Controller.VI_Notch_Filter.TF*exp(-IPNC_Controller.T_s*s);
    [IPNC_Controller.Idq_Regulator.Gm, IPNC_Controller.Idq_Regulator.Pm, IPNC_Controller.Idq_Regulator.wgm, IPNC_Controller.Idq_Regulator.wpm] = margin(IPNC_Controller.Idq_Regulator.G_OpenLoop_TF);
    if(verbose)
        figure(fignum)
        fignum = fignum+1;
        margin(IPNC_Controller.Idq_Regulator.G_OpenLoop_TF);
        disp(['IPNC Idq Regulator: Gm = ' num2str(IPNC_Controller.Idq_Regulator.Gm) ' @ ' num2str(IPNC_Controller.Idq_Regulator.wgm/2/pi) ' Hz, Pm = ' num2str(IPNC_Controller.Idq_Regulator.Pm) ' deg @ ' num2str(IPNC_Controller.Idq_Regulator.wpm/2/pi) ' Hz']);
        title('IPNC Current Regulator Open-Loop Transfer Function');
    end
end

