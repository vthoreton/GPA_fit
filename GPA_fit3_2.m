% GPA_fit v3.2 15/04/2021
% by  Vincent Thoréton (email address for questions, comments and requests: vincent.thoreton.smn.uio.no (UiO) or  work.vincent@thoreton.net (permanent) )
% This script opens raw data (ascii files exported by QUADSTAR 32-Bit), plots and processes the data in a semi-automatic way.
% It has been written for a research purpose but also for an educational
% purpose, as the code is explained to a great extent. You may copy it
% or take inspiration from it to make your own routines. Suggestions and 
% requests are very welcome!
% Acknowledgements are also welcome if you are using this script for
% publication :-)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   List of associated functions: %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   External functions:
%
%   -import_GPA_QS32 (imports raw data files from QUADSTAR 32-Bit)
%
%   -reconstruct_data_file.m (homogeneise the number of columns of the data file)
%
%   -savloarch_log.m (to load, save and archive log files)
%
%   -analysis_log.m (to export vectors in columns)
%
%   -exponential_decay_extract_tau_and_delay.m (fitting tau, delay or both)
%
%   -Den_Otter_two_steps.m (fitting to the 2 steps model; 1 M. W. den Otter, B. A. Boukamp and H. J. M. Bouwmeester, Theory of oxygen isotope exchange, Solid State Ionics, 2001, 139, 89–94.
%
%   -Boukamp_k.m (fitting of the dissociation coefficient k_dis according
%   to Boukamp's theory in: 1 B. A. Boukamp, H. J. M. Bouwmeester and
%   A. J. Burggraaf, The surface oxygen process in oxygen ion conducting
%   materials, Proc. 2nd Intl. Symp. Ion. Mix. Conduct. Oxides, 1994, 141–149.
%
%   -Mizusaki_kD.m (fitting of k and D (as well as a delay and the particle
%   radius) in Mizusaki model (1 T. Ishigaki, S. Yamauchi, J. Mizusaki, K.
%   Fueki and H. Tamura, Tracer diffusion coefficient of oxide ions in
%   LaCoO3 single crystal, J. Solid State Chem., 1984, 54, 100–107.) based
%   on Edwards's equations (1 J. T. B. H.S. Edwards, A.F. Rosenberg, Report NO. ASD-TDR-63-635, 1963.))
%
%   -Mizusaki_D.m (same as previous but fitting D only (and possibly a delay and the particle radius))
%
%   -Mizusaki_D_roots.m (calculates the roots of tan(x)=3*x/(3+alpha*x^2),
%   used for Mizusaki_D.m model)
%
%   -Mizusaki_kD_roots.m (calculates the roots of tan(x)=(3M*x-alpha*x^3)/(3M+(M-1)*alpha*x^2),
%   used for Mizusaki_kD.m model)
%
%	-bg_corr.m (experimental: trying to refine the background of masses 32,
%	34 and 36 to have a nicer exponential decay fit)

%   Internal functions (routines at the end of this file):
%
%   -restrain_time_range (quick choice of a time range, giving some preset options)
%
%   -basic_restrain_time_range (same but no preset option)
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   LIST OF POTENTIAL IMPROVEMENTS: %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	*   close intermediary useless figures ** OK**
%   *   chose between using particle size or SSA in the calculation **OK**
%   *   choice  between using the initial mass (no labelling yet OR the final mass (with 18O enrichment) )  **OK**
%   *   define the initial content of 18O in the powder **OK**
%   *   import the experimental conditions by recognising the content, independently of the number and order of the parameters. **OK**
%   *   clean the code **in progress**
%   *   export the data as csv **OK**
%   *   export the fitted parameters as csv **OK, improved**
%   *   implement the sphericity in the calculations (MEH!)
%   *   make the main database great again *** sort of OK
%   *   merge exponential_decay_extract_tau.m and exponential_decay_extract_tau_and_delay.m **OK**
%   *   improve the way figures are saved (with a list and eval) ** sort of ready but not implemented **
%   *   improve the way the logged values are exported (make it more elegant/universal) **OK**
%   *   include Klier fitting
%   *   ULTIMATE: make a GUI
%   *   ...
%
% clear the workspace, figures etc:
clear all;
clearvars;
close all force;
% declare the global variables that may be shared with subroutines and functions:
global frac_18_O2_inf norm_36_t0 p corr_time frac_18_O2_t0 q r_s A B S r_dis_guess norm_36 delta_f18 P V R T nu m M norm_norm_36 norm_36_inf possible_mass_num possible_mass_str raw_time raw_data_headers raw_data recorded_mass_n_columns experiment_name experiment_date experiment_time frac_18_O frac_18_O2 MS_pressure extract_tau_data two_steps_model_data experimental_conditions_log rho SSA time_origin t_max experimental_conditions_log_path background_n_yields_database timestamp raw_data_path raw_data_filename raw_data_root_folder experimental_conditions raw_time_min sample raw_time_min f18_sample_t0 M_16O M_18O t_max_plot_scale c_o fitting_data nroots fitting_param D_k_delay_and_radius_calc guess radius bg_corr_data custom_indexes custom_message SSA_log radius_log sphericity m_i m_f low_index high_index log_m_i_instead m_f_log m_i_log corr_time_min norm_16 norm_17 norm_18 norm_19 norm_20 norm_21 norm_22 norm_23 norm_28 norm_30 norm_32 norm_34 norm_36 norm_40 norm_44 norm_46 norm_48 frac_18_O2 frac_18_O frac_18_O2_calc_Mizu_D frac_18_O2_calc_Mizu_kD frac_18_O2_calc_DenOtter f_32 f_34 f_36 f_32_bis f_34_bis f_36_bis norm_32_calc norm_34_calc norm_36_calc y regLC regLCfull m_f_calc m_i_calc SSA_calc radius_calc sphericity D_Mizu_D delay_Mizu_D k_Mizu_kD D_Mizu_kD L_C_Mizu_um delay_Mizu_kD radius_Mizu_kD ks_DO kdis_DO delay_DO p_DO ks_DO_free kdis_DO_free delay_DO_free p1_DO_free p2_DO_free ks kdis r_dis_fit nu_calc ;






%say hi
disp('Hello World!');
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ********************************************************************* %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ***************** %% Choice of parameters and options %% ************ %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ********************************************************************* %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%
%% Data handling parameters:
discard_pressure=0; % if the pressure needs to be modified...
skip_plots=0;% skip plotting multiple figures at the beginning and go straight to fitting
override_18O_fractions=0;%skip the starting and final fraction of 18O indicated in the experimental_condition_log.csv file
override_time_window=0;%skip the time window indicated in the experimental_condition_log.csv file
override_volume=0;% to override the volume (not so useful, but gives the option to set a new volume without editing the log file)
force_initial_mass=0;% By default, the final mass is used (including 18O enrichment)for usual GPA performed at chemical equilibrium. Using the initial mass instead can be useful for titration or if the final mass is missing.
force_delta_to_zero=0; % 
calculate_SSA_from_radius=0;% does as the variable says. The values from the log are kept.
calculate_radius_from_SSA=0;% does as the variable says; obviously, these two cannot be both set to 1...
plot_frac_18_O=0; % plotting the 18O fraction determined from mass 16 and 18 (additionally to the one determined from masses 32, 34 and 36).
plot_F=0; % plotting function F (see Boukamp's paper for details)
t_min_plot_scale=-5;% draws plots with additional time before the time origin (minutes. It is useful to visualise the initial isotopologues' fractions.
diff_step=4;% used in the drawing of differentials: interval from the
% central datapoint to each of the two datapoints on either sides used to
% calculate the local differentials. Should be minimised if possible but
% increased for noisy signals.% default was 10
%
%% Fitting parameters:
options_fit=optimset('Display','notify','TolX',1e-11,'TolFun',1e-11,'MaxIter',1000);% default optimisation of the fminsearch function used for fitting

% Parameters to refine for Mizusaki's (Edwards) fitting of D only (diffusion limitation):
% !! no choice of parameters !! (it might be implemented later) D and a delay are fitted.
% However the script can possibly:
% 1*fit D only
% 2*fit D and a delay (this is the chosen mode)
% 3*fit D, a delay and the particle size

% Parameters to refine for Mizusaki's (Edwards) fitting of k and D (mixed regime):
fitD=1;% 1 by default
fitk=1;% 1 by default
fitdelay=1;% 1 by default
fitradius=0;% 0 by default
nroots=24;% number of roots used in Mizusaki models. 16 was found enough.

%Options for Den Otter's approach:
fit_tau=1;% should always be 1
delay_correction=1;% Fitting of the time at which the sample starts
% exchanging. May be different from the time at which the sample
%reaches the exchange temperature or at least the time at which the
%information reaches the MS. It is critical in the case of fast exchange
%where a relatively large delay appears at the begining. This option allows
%fitting that delay, thus getting a better fit. It recommended to keep it
%set to 1. Set it to 0 e.g. if you get a negative delay despite a meaningful
%choice of the time origin. You should then be able to set the delay to a
%value more suitable with the situation.
plot_tau=0;%optional figure of the normalised 18O fraction
advanced_DO_fitting=1; % for further fitting with the Den Otter's 2 steps model
% additional parameters for further fitting with the Den Otter's 2 steps model:
fit_p1=1; % if fit_p1 is set to 1 and fit_p2 to 0, p1=p2 ; if fit_p1 is set to 0 and fit_p2 to 1, only one is fitted (NB: p1 and p2 are mathematically equivalent); if both are set to 1, both are fitted independently; if both are 0, none is fitted
fit_p2=1;
fit_delay_DO=0; % high probablility of failure if set to 1... not sure why
fit_tau_DO=1;

% Options for Boukamp's approach:
k_as_a_function_of_time=0; % plots k as a function of time from the slope of the curve; this is for OBSERVATION ONLY because this is wrong to do it (variable k violates the model!)
diff_step_Boukamp=2; % used fot the determination of k as a function of time (see diff_step)
fiftyfiftymix_option=0;% gives the choice of setting f18(g) and f18(s) to 0.5 in Boukamp's model (default is 0). Set it to 1 if you want to do some "old school" 50/50 GPA.
quartz_background_fit=0;% gives the choice of fiting 18O2 fraction (mass 36) as an exponential decay; can be used when fitting the quartz background (~deprecated)
%
%% Exporting the data:
%Local export of vectors and figures:
save_norm_data=1;%This option allows to save the raw data (after background correction and normalisation)
save_all_in_one=0;% The vectors related to the Mizusaki, Den Otther or Boukamp fitting are exported separatelly for convenience. This option allows to save all the vectors in one single big file.
fig_res=300;% resolution of exported figures in dpi
%Export to the main database:
option_log_to_database=0;% Offers to log to the "main database" ### BROKEN! ###
saving_folder=pwd;%working directory of matlab
% or somewhere else, e.g. : saving_folder = 'C:\Users\username\Documents';
GPA_log_db_path = strcat(saving_folder,'\GPA_log_db.csv');%location of the external main database
%
%% Experimental features:
option_update_weigth_and_nu=0;% offers the option to replace the molar mass and oxygen stoichiometry with the one calculated (not recommended)
mass_transport_correction=0;% mass tranport correction (experimental, not recommended)
background_correction=0;% background intensity refinement of masses 32 34 and 36 (experimental, can work sometimes (ONLY USED IN Den Otter fitting) if the signal is very noisy)... might correct for data shift when the leak valve is too much closed during the measurement (experimental)
carbotruc=0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ********************************************************************* %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% * %%

%% * %%

%% * %%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Program starts %%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a few constants %
%%%%%%%%%%%%%%%%%%%%%%%%%%
R=8.314462;% ideal gas constant
M_16O=15.9949146196; % atomic weight of oxygen 16
M_18O=17.99915961;% atomic weight of oxygen 18
frac_18_O_nat=0.00205;% Natural abundance of oxygen 18 ; not to get mixed up with the initial fraction of 18O in the sample.



%%%%%%%%%%%%%%%%%%%%
%% Default values %%
%%%%%%%%%%%%%%%%%%%%
default_frac_18=0.964;%depends on the labelled gas
default_frac_36=0.93; %depends on the labelled gas
%% background correction and ionisation yields database for auto-picked masses:
background_n_yields_database=[14 16 17 18 19 20 21 22 23 24 28 30 32 34 36 40 44 46 48;
    0 1.34E-12 4.328E-13 7.8E-13 3.42E-13 3.38E-13 0 0 0 0 3.459E-47 0 4.37E-13 3.23E-13 1.90841E-13 3.09e-13 3.33e-13 0 0;%% must fill the missing values
    1 1 1 1 1 1 1 1 1 1 1.3315 1 1 1 1 1.28 1.12 1.12 1.12];%% must fill the missing values if relevant (CO2 and water in particular)


%%%%%%%%%%%%%%%%%%%%
%% Importing data %%
%%%%%%%%%%%%%%%%%%%%
timestamp=datestr(now,'mmmm-dd-yyyy HHMMSS'); %timestamp used for naming log folders
%choose the raw data file (directly exported from the MS software)
[raw_data_filename,raw_data_root_folder] = uigetfile('*.*');
raw_data_path = strcat(raw_data_root_folder,raw_data_filename);
message=msgbox('Importing data, please wait...','Notification');
recorded_mass_n_columns=import_GPA_QS32(raw_data_path);%call import_GPA_QS32 function
%checking the consistency of 2 different mass lists
if prod(background_n_yields_database(1,:)==possible_mass_num)~=1
    msgbox('Inconsistent mass list; might cause a problem! Script aborted.','Warning!')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Split raw data to vectors, correct for background and normalise to O2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_time_min=raw_time/60;
for i = 1:size(recorded_mass_n_columns,2)
    eval(sprintf('raw_intensity%g=raw_data(1:end,recorded_mass_n_columns(2,i));',recorded_mass_n_columns(1,i)));
    col=find(background_n_yields_database(1,:)==recorded_mass_n_columns(1,i));
    eval(sprintf('corr_intensity%g=(raw_intensity%g-background_n_yields_database(2,col))/background_n_yields_database(3,col);',recorded_mass_n_columns(1,i),recorded_mass_n_columns(1,i)));
end
sum_intensityO2=corr_intensity32+corr_intensity34+corr_intensity36;
for i = 1:size(recorded_mass_n_columns,2)
    eval(sprintf('norm_%g=corr_intensity%g./sum_intensityO2;',recorded_mass_n_columns(1,i),recorded_mass_n_columns(1,i)));
end
delete(message); %after data got sucessfully imported


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve or assign experimental conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine the log and figures folder
logging_folder=sprintf('%s',raw_data_filename(1:end-4));% removes the .ext
saving_path=strcat(raw_data_root_folder,logging_folder,'\');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% locate and reload the experimental conditions... %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

old_parameters_list={'sample' 'raw_data_filename'; % for compatibility with former log version       
        'm_f' 'm_f_log';
        'm_i' 'm_i_log';
        'M' 'M';
        'nu' 'nu';
        'rho' 'rho';
        'SSA'  'SSA_log';
        'V' 'V';
        'P' 'P';
        'T' 'T';
        't0'   'time_origin' ;
        'tmax' 't_max';
        't_max_scale' 't_max_plot_scale';
        'f18_t0' 'frac_18_O2_t0';
        'f18_inf'    'frac_18_O2_inf';
        'f36_t0'  'norm_36_t0';
        'f36_inf'  'norm_36_inf';
        'particle_radius'  'radius_log';
        'sphericity' 'sphericity';
        'spher'  'sphericity';
        'f18_sample_t0' 'f18_sample_t0' ;
        'm' 'm' };
    
    new_parameters_list={'sample'  'raw_data_filename'; % current version (more readable with 1 parameters per row)
        'm_f [g]' 'm_f_log';
        'm_i [g]' 'm_i_log';
        'M [g/mol]' 'M';
        'nu' 'nu';
        'rho [g/m³]' 'rho';
        'SSA [m²/g]'  'SSA_log';
        'V [m³]' 'V';
        'P [Pa]' 'P';
        'T [K]' 'T';
        't_0 [minutes]'   'time_origin' ;
        't_max [minutes]' 't_max';
        't_max_scale [minutes]' 't_max_plot_scale';
        'f18_t0' 'frac_18_O2_t0';
        'f18_inf'    'frac_18_O2_inf';
        'f36_t0'  'norm_36_t0';
        'f36_inf'  'norm_36_inf';
        'particle_radius [m]'  'radius_log';
        'sphericity' 'sphericity';
        'f18_sample_t0' 'f18_sample_t0'};

experimental_conditions_log_path = strcat(saving_path,'experimental_conditions_log.csv');
if exist(experimental_conditions_log_path,'file') == 2
    experimental_conditions=readtable(experimental_conditions_log_path);
    experimental_conditions_content=experimental_conditions.Properties.VariableNames;%get the names of variables from the log
    experimental_conditions=table2cell(experimental_conditions);%converts the table to array before extracting the values
    
if size(experimental_conditions,1)==1
%% old format in 2 lines!
    parameters_to_load=old_parameters_list;    
    elseif size(experimental_conditions,2)==2    
        parameters_to_load=new_parameters_list;
    end
    parameters_to_save=new_parameters_list;    
 savloarch_log(parameters_to_load,experimental_conditions_log_path,'load'); 
    
    % for compatibility with former log versions:
    if isa(T,'double')~=1
        T=301;
    end

    if isa(sphericity,'char')==1
        sphericity=1;
    end
    
    if isa(f18_sample_t0,'char')==1
        f18_sample_t0=frac_18_O_nat;
    end
    
    if isempty(m)==0
        m_f_log=m; % as m was usually the final mass
    end
    
    if isa(m_f,'double')==1 % it works because m_f is either a double or it is defined as 'NA'
        m_f=m_f_log;
    end
    
    if isa(m_i,'double')==1
        m_i=m_i_log;
    end
    
    if isa(f18_sample_t0,'char')==1 %% for compatibility with previous files that do not include f18_sample_t0
        f18_sample_t0=frac_18_O_nat;
    end
    
    if strcmp(radius_log,'NA')==1%for setting automatically the radius in former files where only the SSA was indicated
        radius_calc=3*sphericity/SSA_log/rho;
        radius=radius_calc;
        radius_log=radius;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OR input parameters in a dialog %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %defaults parameters
    m=0.020;%sample mass (g) preferably AFTER the exchange (except if indicated in the parameters). (we assume that the sample hase been degased and that the sample contains a significant amount of 18O) 0.0358 0.0401
    M=136.74;% M is the molar mass of the material with natural content of 18O with ideal oxygen (g/mol). EVEN IF the sample is initially labelled, this must not be taken into account (the script takes it into account).
    SSA_log=0.01;%sample surface area (m^2/g)1.99
    %SSA=SSA_log;
    radius_um=50;%spherical particle radius (um)
    sphericity=1;
    V=103*1E-6;%volume of the tube (m³) 63*1E-6
    nu=3;%number of oxygen per molecular unit
    P=4000;%total pressure in the chanber (Pa) 4000
    rho=4000000;%theoretical density of the phase (g/m³)
    T=302;%Temperature at filling (K)
    V_ml=V*1E6;
    T_C=round(T-273.15,0);
    m_mg=m*1000;
    P_mbar=P/100;
    rho_gcc=rho/1E6;
    time_origin='NA';
    t_max='NA';
    t_max_plot_scale='NA';
    frac_18_O2_t0='NA';
    frac_18_O2_inf='NA';
    norm_36_t0='NA';
    norm_36_inf='NA';
    radius_log=radius_um/1E6;
    f18_sample_t0=frac_18_O_nat;
    
    % input of the experimental parameters
    default_input =[m_mg;M;nu;rho_gcc;f18_sample_t0;SSA_log;radius_um;sphericity;P_mbar;V_ml;T_C];
    experimental_parameters = inputdlg({'m (mg) [it will be chosen at next step if this is the initial of final mass]','M (g/mol) [with ideal oxygen stoichiometry and NO 18O enrichment!]','ideal oxygen stoichiometry','theoretical density (g/cm³)','initial 18O content in the sample [useful fo reused sample; make sure to measure it at the end of the pre-anneal and that it is equilibrated!]','specific surface area (m²/g) [used in Den Otter and Boukamp''s models]','particle radius (um) [used in Mizusaki/Edwards''s models]','sphericity [used in Mizusaki/Edwards''s models but not implemented yet]','Pressure of the reactor AT FILLING (mbar)','Volume of the reactor (mL)','Temperature of the reactor AT FILLING (°C)'},...
        'Input of the experimental parameters',...
        1,string(default_input));
    noideawhy=size(experimental_parameters);
    if noideawhy(1)~=0 %actually it is because the size is 0 if canceled is pressed
        experimental_parameters=str2double(experimental_parameters);
        m_mg=experimental_parameters(1);
        m=m_mg/1E3;
        M=experimental_parameters(2);
        nu=experimental_parameters(3);
        rho_gcc=experimental_parameters(4);
        rho=rho_gcc*1E6;
        f18_sample_t0=experimental_parameters(5);
        SSA_log=experimental_parameters(6);
        radius_um=experimental_parameters(7);
        radius_log=radius_um/1E6;
        sphericity=experimental_parameters(8);
        P_mbar=experimental_parameters(9);
        P=P_mbar*100;%Pa
        V_ml=experimental_parameters(10);
        V=V_ml/1E6;
        T_C=experimental_parameters(11);
        T=T_C+273.15;
    elseif size(experimental_parameters)==0  % if cancel is pressed
        choice = questdlg('Continue with default values?',...
            'Warning', ...
            'yes','no','no');
        switch choice
            case 'no'
                return
            case 'yes'
        end
    end
    
    question = questdlg('Is the mass:', ...
        'Choice of mass', ...
        'the final mass (default)','the initial mass','the final mass (default)');
    switch question
        case 'the final mass (default)'
            force_initial_mass=0;
        case 'the initial mass'
            force_initial_mass=1;
    end
    
    if force_initial_mass==1
        m_i=m;
        m_mg_i=m_i*1000;
        m_i_log=m_i;
        m_f_log='NA';
        %log_m_i_instead=1;
    else % default
        m_f=m;
        m_mg_f=m_f*1000;
        m_f_log=m_f;
        m_i_log='NA';
    end
    
    
    %save  the experimental conditions log if not there yet:
    mkdir(fullfile(raw_data_root_folder,logging_folder));%create the logging folder
    parameters_to_save=new_parameters_list;
    %save_conditions();%save
end % end of input of parameters (saved of dialog)

savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save


if override_volume==1;
    override_vol = questdlg('Override the volume?', ...
        'Volume override', ...
        'yes','no','yes');
    switch scale_limit
        case 'yes'
            input_new_volume = inputdlg({'new volume (mL)'},...
                '',...
                1);
            input_new_volume=str2double(input_new_volume);
            V=input_time_limit(1)/1E6;
        case 'no'
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choosing the area of interest and setting the time origine %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frac_18_O2=norm_36+1/2*norm_34;
if strcmp(time_origin,'NA')==1 || strcmp(t_max,'NA')==1 || override_time_window==1% OR
    % initial display (raw data)
    plot_raw323436=figure('Name','Raw signals','NumberTitle','off');
    hold on;
    plot(raw_time_min,corr_intensity32,'b-','Linewidth',2)
    plot(raw_time_min,corr_intensity34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
    plot(raw_time_min,corr_intensity36,'r-','Linewidth',2)
    legend('I_3_2','I_3_4','I_3_6','Location','east')
    xlabel('time (min)')
    ylabel('Current (A)')
    
    % display of the normalised data
    plot_norm323436=figure('Name','Normalised signals','NumberTitle','off');
    hold on;
    plot(raw_time_min,frac_18_O2,'k-','Linewidth',2)
    plot(raw_time_min,norm_32,'b-','Linewidth',2)
    plot(raw_time_min,norm_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
    plot(raw_time_min,norm_36,'r-','Linewidth',2)
    xlabel('time (min)')
    ylabel('Isotopologue fraction f_i')
    axis([0 inf 0 1])
    if ismember(28,recorded_mass_n_columns)==1
        choice = questdlg('Plot mass 28?', ...
            'Checking leakage (nitrogen)', ...
            'yes','no','no');
        switch choice
            case 'no'
                legend('f_1_8','f_3_2','f_3_4','f_3_6','Location','east')
            case 'yes'
                yyaxis right
                set(gca,'YColor',[0 1.0 0])
                semilogy(raw_time_min,norm_28,'color',[0 1.0 0],'Linewidth',2)
                ylabel('N_2/O_2 (OR CO/O_2...)')
                axis([0 inf 1e-3 1e-1])
                legend('f_1_8','f_3_2','f_3_4','f_3_6','f_2_8','Location','east')
        end
        
    end
    % choice of the zone of interest
    custom_message='Zoom in the range of interest (e.g. discard the pre-anneal if it was recorded in the same file)';
    [low_index high_index]=basic_restrain_time_range(raw_time_min);
    time_origin=raw_time_min(low_index);
    t_max=raw_time_min(high_index)-time_origin;
    % shrinks the dataset accordingly to the chosen range
    if exist('MS_pressure')==1 & size(MS_pressure)~=0
        MS_pressure=MS_pressure(low_index:high_index,1);
    end
    raw_time=raw_time(low_index:high_index,1);
    raw_time_min=raw_time_min(low_index:high_index,1);
    for i = 1:size(recorded_mass_n_columns,2)
        eval(sprintf('norm_%g=norm_%g(low_index:high_index,1);',recorded_mass_n_columns(1,i),recorded_mass_n_columns(1,i)));
    end
    corr_time_min=raw_time_min-time_origin;% center the window to the range of interest
    corr_time=corr_time_min*60;
    frac_18_O2=frac_18_O2(low_index:high_index,1);
    
    
    % Setting the time origin and fitting range
    %% the evolution of the pressure helps to set the time origin!
    plot_norm323436shift=figure('Name','Reduced range','NumberTitle','off');
    if exist('MS_pressure')==1 & numel(MS_pressure)~=0
        %% should handle that better, but it should work for most of the files
        yyaxis right
        semilogy(corr_time_min,MS_pressure,'color',[0 1.0 0],'Linewidth',2)
        ylabel('MS pressure')
        if max(MS_pressure)<=1e-8
            axis([t_min_plot_scale t_max 1e-9 1e-8])
        elseif    max(MS_pressure)<=1e-7
            axis([t_min_plot_scale t_max 1e-8 1e-7])
        elseif max(MS_pressure)<=1e-6
            axis([t_min_plot_scale t_max 1e-7 1e-6])
        elseif max(MS_pressure)<=1e-5
            axis([t_min_plot_scale t_max 1e-6 1e-5])
        elseif max(MS_pressure)<=1e-4
%             axis([t_min_plot_scale t_max 1e-5 1e-4])
%         elseif max(MS_pressure)<=1e-3
            axis([t_min_plot_scale t_max -Inf Inf])
        end
        title('Reduced range')
        set(gca,'YColor',[0 1.0 0])
    end
    yyaxis left
    hold on;
    plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
    plot(corr_time_min,norm_32,'b-','Linewidth',2)
    plot(corr_time_min,norm_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
    plot(corr_time_min,norm_36,'r-','Linewidth',2)
    xlabel('time (min)')
    ylabel('Isotopologue fraction f_i')
    set(gca,'YColor',[0 0 0])
    axis([t_min_plot_scale t_max 0 1])
    if exist('MS_pressure')==1 & size(MS_pressure) ~=0
        legend('f_1_8','f_3_2','f_3_4','f_3_6','MS pressure','Location','east')
    else
        legend('f_1_8','f_3_2','f_3_4','f_3_6','Location','east')
    end
    custom_message='Select the time bounds(the change of pressure can help to define t_0)';
    [low_index high_index]=basic_restrain_time_range(corr_time_min);
    time_origin=corr_time_min(low_index)+time_origin;
    t_max=corr_time_min(high_index)-corr_time_min(low_index);
    
    corr_time_min=corr_time_min-corr_time_min(low_index);
    corr_time=corr_time_min*60;%seconds
    average_MS_pressure=mean(MS_pressure);
    for i = 1:numel(recorded_mass_n_columns)/2
        eval(sprintf('norm_%g=norm_%g;',recorded_mass_n_columns(1,i),recorded_mass_n_columns(1,i)));
    end
    %save_conditions();
    savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save

else    % Setting automatically the origin time from the experimental conditions file
    s=size(raw_time,1);
    low_index=1;
    for a =1:s
        if raw_time_min(a)>time_origin+t_min_plot_scale
            break
        else
            low_index=a;
        end
    end
    low_index_close_cut=1;
    for a =1:s
        if raw_time_min(a)>time_origin
            break
        else
            low_index_close_cut=a;
        end
    end
    
    t_end=time_origin+t_max;
    for a = s:-1:1
        if raw_time_min(a)<t_end
            high_index=a;
            break
        end
    end
    corr_time_min=raw_time_min(low_index:high_index,1)-time_origin;
    corr_time=corr_time_min*60;%seconds
    if exist('MS_pressure')==1 & size(MS_pressure) ~=0
        MS_pressure=MS_pressure(low_index:high_index,1);
    end
    average_MS_pressure=mean(MS_pressure(low_index_close_cut-low_index:end,1));
    for i = 1:size(recorded_mass_n_columns,2)
        eval(sprintf('norm_%g=norm_%g(low_index:high_index,1);',recorded_mass_n_columns(1,i),recorded_mass_n_columns(1,i)));
    end
    for i = 1:size(recorded_mass_n_columns,2)
        eval(sprintf('raw_intensity%g=raw_intensity%g(low_index:high_index,1);',recorded_mass_n_columns(1,i),recorded_mass_n_columns(1,i)));
    end
end % end of choosing the zone of interest and time origin

% plot the restricted zone
% calculation of the 18O fraction from O2 or O (worth comparing!)
frac_18_O2=norm_36+1/2*norm_34;
if plot_frac_18_O ==1 && exist('norm_18')==1
    frac_18_O=norm_18./(norm_18+norm_16);
else
    plot_frac_18_O =0;
end
F=norm_34./(norm_32+norm_36);

plot_norm323436shift=figure('Name','Fitting range','NumberTitle','off');
hold on;
plot(corr_time_min,norm_32,'b-','Linewidth',2)
plot(corr_time_min,norm_34,'-','color',[0.5  0.0  0.5],'Linewidth',2)
plot(corr_time_min,norm_36,'r-','Linewidth',2)
plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
xlabel('time (min)')
ylabel('Isotopologue fraction f_i')
axis([t_min_plot_scale t_max 0 1])
legend('f_{32}','f_{34}','f_{36}','^{18}O fraction','Location','east')


%Set the starting and final fractions by hand or automatically (or keep the values from the file)
if strcmp(frac_18_O2_t0,'NA')==1 || strcmp(frac_18_O2_inf,'NA')==1 || strcmp(norm_36_t0,'NA')==1 || strcmp(norm_36_inf,'NA')==1 || override_18O_fractions==1
    choice = questdlg('Set the starting and final fractions by hand?', ...
        'Limit conditions', ...
        'yes','no (auto)','yes');
    switch choice
        case 'no (auto)'
            frac_18_O2_t0=max(frac_18_O2);
            frac_18_O2_inf=min(frac_18_O2);
            norm_36_t0=max(norm_36(:,1));
            norm_36_inf=min(norm_36(:,1));
        case 'yes'
            choice2 = questdlg('Setting initial 18O and 36O2 fractions:', ...
                'Starting  level', ...
                'Set visually','Type','Set visually');
            switch choice2
                case 'Set visually'
                    message = sprintf('Select the initial 18O fraction (black curve inset)');
                    uiwait(msgbox(message,'Notification'));
                    frac_18_O2_bounds=ginput(1);
                    frac_18_O2_t0=max(frac_18_O2_bounds(:,2));
                    message = sprintf('Select the initial 18O2 fraction (red curve inset)');
                    uiwait(msgbox(message,'Notification'));
                    fraction_36_bounds=ginput(1);
                    norm_36_t0=max(fraction_36_bounds(:,2));
                case 'Type'
                    frac_18_O2_t0=default_frac_18;
                    norm_36_t0=default_frac_36;
                    default_input =[frac_18_O2_t0;norm_36_t0];
                    experimental_parameters = inputdlg({'18O fraction','36 fraction'},...
                        '',...
                        1,string(default_input));
                    check_replied=size(experimental_parameters);
                    if check_replied(1)==2 %size is 0 if cancel is pressed
                        experimental_parameters=str2double(experimental_parameters);
                        frac_18_O2_t0=experimental_parameters(1);
                        norm_36_t0=experimental_parameters(2);
                    elseif size(experimental_parameters)==0  % if cancel is pressed
                        choice3 = questdlg('Continue with default values?',...
                            'Warning', ...
                            'yes','no','no');
                        switch choice3
                            case 'no'
                                return
                            case 'yes'
                        end
                    end
            end
            
            choice4 = questdlg('Setting final 18O and 36O2 fractions:', ...
                'Starting  level', ...
                'Set visually','Type','Calculate','Set visually');
            switch choice4
                case 'Calculate'
                    %% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %% This part is problematic and may be investigated properly... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    M_iO=M_18O*f18_sample_t0+M_16O*(1-f18_sample_t0);
                    A=2*P*V/R/T;
                    
                    if isa(m_i_log,'double')==1 && isempty(m_i_log)==0
                        m_f_approx=m_i_log*1.005; %% this is not right but "ok" for now
                        msgbox('Impossible to calculate properly the final fraction as the initial mass is used (calculation loop)!','Warning!') %Script aborted.
                        %return
                    elseif isa(m_f_log,'double')==1 && isempty(m_f_log)==0
                        m_f_approx=m_f_log;
                    end
                    
                    delta_hypo=0;
                    M_fO=0.5;
                    % converging M_fO and frac_18_O2_inf: (5 steps should be enough)
                    for i=1:10;
                        frac_18_O2_inf=(A*((M+nu*(f18_sample_t0-frac_18_O_nat))-delta_hypo*M_fO-(M_iO-M_fO)*nu)*frac_18_O2_t0-(delta_hypo-nu)*f18_sample_t0*m_f_approx)/(A*((M+nu*(f18_sample_t0-frac_18_O_nat))-delta_hypo*M_fO-(M_iO-M_fO)*nu)-(delta_hypo-nu)*m_f_approx);
                        norm_36_inf=frac_18_O2_inf^2;
                        M_fO=(1-frac_18_O2_inf)*M_16O+frac_18_O2_inf*M_18O;%
                    end
                    
                    
                case 'Set visually'
                    message = sprintf('Select the final 18O fraction (black curve end)');
                    uiwait(msgbox(message,'Notification'));
                    frac_18_O2_bounds=ginput(1);
                    frac_18_O2_inf=frac_18_O2_bounds(2);
                    message = sprintf('Select the final 18O2 fraction (red curve end)');
                    uiwait(msgbox(message,'Notification'));
                    fraction_36_bounds=ginput(1);
                    norm_36_inf=fraction_36_bounds(2);
                case 'Type'
                    frac_18_O2_inf=0.2;
                    norm_36_inf=0.15;
                    default_input =[frac_18_O2_inf;norm_36_inf];
                    experimental_parameters = inputdlg({'18O fraction','36 fraction'},...
                        '',...
                        1,string(default_input));
                    check_replied=size(experimental_parameters);
                    if check_replied(1)==2 %size is 0 if cancel is pressed
                        experimental_parameters=str2double(experimental_parameters);
                        frac_18_O2_inf=experimental_parameters(1);
                        norm_36_inf=experimental_parameters(2);
                        
                    elseif size(experimental_parameters)==0  % if cancel is pressed
                        choice4 = questdlg('Continue with default values?',...
                            'Warning', ...
                            'yes','no','no');
                        switch choice4
                            case 'no'
                                return
                            case 'yes'
                        end
                    end
                    
                    
            end
            %save_conditions();
            savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save

    end  %end of setting the limit condtions by hand or automatically
else
    waitfor(msgbox('Using the recorded values of the f18 and f36 bounds.','Notification'));
end

if strcmp(t_max_plot_scale,'NA')==1 || t_max_plot_scale==t_max
    scale_limit = questdlg('Constrain the time scale?', ...
        'Time scale limit', ...
        'yes','no','no');
    switch scale_limit
        case 'yes'
            input_time_limit = inputdlg({'time limit (min)'},...
                '',...
                1);
            input_time_limit=str2double(input_time_limit);
            t_max_plot_scale=input_time_limit(1);
        case 'no'
            t_max_plot_scale=t_max;
    end
end
close(plot_norm323436shift);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating useful values and analytical calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% message about the use of m_i if it is the case:
if force_initial_mass==1 || (isa(m_i_log,'double')==1 && isempty(m_i_log)~=1) %% if m_i is used
    uiwait(msgbox('Using the initial mass instead of the final mass. It is recommended to use the final mass (at equilibrium at the defined pO2) except if you are using the script for titration at room temperature; if so, ignore this message.','Warning!'));
    log_m_i_instead=1;
    m=m_i_log; %  WRONG, as m is actually m_f in all the rest of the script!!! very bad idea to use m as it is not well defined
elseif isa(m_f,'double')==1
    log_m_i_instead=0;
    m=m_f_log; % m should be equal to m_f in the whole script
end

A=2*P*V/R/T;%amount of oxygen in the gas phase (mol)
M_iO=M_18O*f18_sample_t0+M_16O*(1-f18_sample_t0); % M_iO=15.999 at natural 18 abundance %% OK
M_fO=(1-frac_18_O2_inf)*M_16O+frac_18_O2_inf*M_18O;% averaged atomic weight of oxygen after labelling %% OK

if force_initial_mass==1 || (isa(m_i,'double')==1 & isempty(m_i)~=1)==1  % m is the final mass by default BUT:
    delta=(A*(M+nu*(f18_sample_t0-frac_18_O_nat))*(frac_18_O2_t0-frac_18_O2_inf)-(frac_18_O2_inf-f18_sample_t0)*m_i*nu)/(A*(frac_18_O2_t0-frac_18_O2_inf)*M_iO-(frac_18_O2_inf-f18_sample_t0)*m_i);%OK
if force_delta_to_zero==1
delta=0;
end
else % default
    delta=(A*((M+nu*(f18_sample_t0-frac_18_O_nat))-(M_iO-M_fO)*nu)*(frac_18_O2_t0-frac_18_O2_inf)-(frac_18_O2_inf-f18_sample_t0)*m_f*nu)/(A*(frac_18_O2_t0-frac_18_O2_inf)*M_fO-(frac_18_O2_inf-f18_sample_t0)*m_f);%% OK %precise calculation of delta. Infinitely small matter but still meaningfull!
if force_delta_to_zero==1
delta=0;
end
end

c_o=nu*rho/M;%ideal concentration of oxygen in the solid (mol/m³) %%HAVE TO ASSUME THAT the real nu and M are closed to the ideal values considered here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% c_o=(nu-delta)*real_rho/real_M; at the beginning or at the end (you choose)

% M is the molar mass of the unlabeled material with ideal oxygen
% stoichiometry. NB: if the sample is already labelled at the beginning,
% the script does take it into account.
nu_calc=nu-delta;%% OK
M_calc_i=M+nu*(f18_sample_t0-frac_18_O_nat)-delta*M_iO;%molar mass at start, considering delta (and the initial content of 18O) OK
M_calc_f=M_calc_i+nu_calc*(M_fO-M_iO);% OK
m_increase_by_labelling_pc=100*(M_calc_f-M_calc_i)/M_calc_i; %OK

if calculate_SSA_from_radius==1 && calculate_radius_from_SSA==0 && isa(radius_log,'char')==0  % the case where the radius is not logged and calculate_SSA_from_radius set to 1 would be a stupid combination...
    SSA_calc=3*sphericity/radius_log/rho;
    SSA=SSA_calc;
    SSA_display='SSA (calc)';
elseif calculate_SSA_from_radius==0
    SSA=SSA_log;
    radius_calc=3*sphericity/SSA_log/rho;
    SSA_display='SSA';
end
if calculate_radius_from_SSA==1 && calculate_SSA_from_radius==0
    radius_calc=3*sphericity/SSA_log/rho;
    radius=radius_calc;
    radius_display='particle radius (calc)';
elseif calculate_radius_from_SSA==0
    radius=radius_log;
    SSA_calc=3*sphericity/radius_log/rho;
    radius_display='particle radius';
end
S=m*SSA;%sample surface (m^2)

%human friendly units:
V_ml=V*1E6;
P_mbar=P/100;
rho_gcc=rho*1E-6;
T_C=round(T-273.15,0);%Temperature at filling(*C) and round for easier display
radius_um=radius*1E6;


if isa(m_i,'double')==1 && isempty(m_i)==0%)==1 %force_initial_mass==1 ||
    m_f_calc=m_i*M_calc_f/M_calc_i;
    m=m_f_calc;
    m_mg_i=1E3*m_i;
    m_mg_f=1E3*m_f_calc;
        B=(nu-delta)*m_i/M_calc_i;
    
elseif isa(m_f,'double')==1 && isempty(m_f)==0%)==1
    m_i_calc=m_f*M_calc_i/M_calc_f;
    m_mg_i=1E3*m_i_calc;
    m_mg_f=1E3*m_f;
    

        B=(nu-delta)*m_f/M_calc_f;

end


%%B=nu*m/M;%amount of oxygen in the solid phase (mol) HAVE TO ASSUME THAT nu and M filled there are closed to real value AT THE END !!! the nature of m is not defined properly %% (m is m_f originally)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



alpha=A/B;%ratio of oxygen in the gas/in the solid




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% calculation of the ideal V and of the ideal final 18O fraction IF delta was equal to 0 etc %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
delta_hypo=0;
V_ideal=R*T/2/P*(-(delta_hypo-nu)*(frac_18_O2_inf-f18_sample_t0)*m)/((M-delta_hypo*M_fO-(M_iO-M_fO)*nu)*(frac_18_O2_t0-frac_18_O2_inf)); % m and M may not be right !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
V_ideal_cc=V_ideal*1E6;
P_ideal=R*T/2/V*(-(delta_hypo-nu)*(frac_18_O2_inf-f18_sample_t0)*m)/((M-delta_hypo*M_fO-(M_iO-M_fO)*nu)*(frac_18_O2_t0-frac_18_O2_inf));% m and M may not be right !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
P_ideal_mbar=P_ideal/100;
frac_18_O2_inf_ideal=(A*(M-delta_hypo*M_fO-(M_iO-M_fO)*nu)*frac_18_O2_t0-(delta_hypo-nu)*f18_sample_t0*m)/(A*(M-delta_hypo*M_fO-(M_iO-M_fO)*nu)-(delta_hypo-nu)*m);% m and M may not be right !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
m_f_ideal=(-A*(M-delta_hypo*M_fO-(M_iO-M_fO)*nu)*(frac_18_O2_t0-frac_18_O2_inf))/((delta_hypo-nu)*(frac_18_O2_inf-f18_sample_t0)); % M may not be right !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
m_f_ideal_mg=m_f_ideal*1E3;
T_ideal=(2*P*V)/R*1/((-(delta_hypo-nu)*(frac_18_O2_inf-f18_sample_t0)*m)/((M-delta_hypo*M_fO-(M_iO-M_fO)*nu)*(frac_18_O2_t0-frac_18_O2_inf))); %% M may not be right !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T_C_ideal=T_ideal-273.15;
alpha_ideal=(frac_18_O2_inf-f18_sample_t0)/(frac_18_O2_t0-frac_18_O2_inf);
error_estimate_V_pc=100*(V_ideal_cc-V_ml)/V_ml;
error_estimate_P_pc=100*(P_ideal_mbar-P_mbar)/P_mbar;
error_estimate_m_pc=100*(m_f_ideal_mg-m_mg_f)/m_mg_f;
error_estimate_f_pc=100*(frac_18_O2_inf_ideal-frac_18_O2_inf)/frac_18_O2_inf;
error_estimate_alpha_pc=100*(alpha_ideal-alpha)/alpha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the 18O fraction and the F function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if log_m_i_instead==1
    content_annotation={sprintf('Data file: %s',raw_data_filename),sprintf('Date & time : %s - %s',experiment_date,experiment_time),sprintf('F_i = %0.3g , F^{sample}_i = %0.3g , F_{\\infty} = %0.3g',frac_18_O2_t0,f18_sample_t0,frac_18_O2_inf),sprintf('P_i = %0.4g mbar , V_{reactor} = %0.4g mL , T_i = %0.4g %cC',P_mbar, V_ml,T_C,char(176)),sprintf('m_{i} = %0.4g mg , m^{(calc)}_{\\infty} = %0.4g mg , m_{gain} = %0.2g %%',m_mg_i,m_mg_f,m_increase_by_labelling_pc),sprintf('\\alpha = %0.3g , \\nu_{input} = %0.3g , \\nu_{calc} = %0.3g',alpha,nu,nu_calc),sprintf('M^{input}_i = %0.5g g/mol , M^{calc}_i = %0.5g g/mol',M,M_calc_i),sprintf('%s = %0.3g \\mum , %s = %0.3g m^2/g',radius_display,radius_um,SSA_display,SSA),sprintf(''),sprintf('Debugging:'),sprintf('P^{ideal} = %0.4g mbar , V^{ideal} = %0.4g mL , T^{ideal} = %0.1f %cC',P_ideal_mbar,V_ideal_cc,T_C_ideal,char(176)),sprintf('\\alpha^{ideal} = %0.3g , m_{\\infty}^{ideal} = %0.4g mg , F_{\\infty}^{ideal} = %0.3g',alpha_ideal,m_f_ideal_mg,frac_18_O2_inf_ideal),sprintf('P_{MS} = %0.3g mbar',average_MS_pressure)};
else
    content_annotation={sprintf('Data file: %s',raw_data_filename),sprintf('Date & time : %s - %s',experiment_date,experiment_time),sprintf('F_i = %0.3g , F^{sample}_i = %0.3g , F_{\\infty} = %0.3g',frac_18_O2_t0,f18_sample_t0,frac_18_O2_inf),sprintf('P_i = %0.4g mbar , V_{reactor} = %0.4g mL , T_i = %0.4g %cC',P_mbar, V_ml,T_C,char(176)),sprintf('m^{(calc)}_{i} = %0.4g mg , m_{\\infty} = %0.4g mg , m_{gain} = %0.2g %%',m_mg_i,m_mg_f,m_increase_by_labelling_pc),sprintf('\\alpha = %0.3g , \\nu_{input} = %0.3g , \\nu_{calc} = %0.3g',alpha,nu,nu_calc),sprintf('M^{input}_i = %0.5g g/mol , M^{calc}_i = %0.5g g/mol',M,M_calc_i),sprintf('%s = %0.3g \\mum , %s = %0.3g m^2/g',radius_display,radius_um,SSA_display,SSA),sprintf(''),sprintf('Debugging:'),sprintf('P^{ideal} = %0.4g mbar , V^{ideal} = %0.4g mL , T^{ideal} = %0.1f %cC',P_ideal_mbar,V_ideal_cc,T_C_ideal,char(176)),sprintf('\\alpha^{ideal} = %0.3g , m_{\\infty}^{ideal} = %0.4g mg , F_{\\infty}^{ideal} = %0.3g',alpha_ideal,m_f_ideal_mg,frac_18_O2_inf_ideal),sprintf('P_{MS} = %0.3g mbar',average_MS_pressure)};
end
plot_frac_18_O2=figure('Name','Summary of experimental conditions and titration','NumberTitle','off');
hold on;
set(gcf, 'Position',  [100, 100, 800, 600],'color','w')
plot(corr_time_min,frac_18_O2,'b-','Linewidth',2)
if plot_frac_18_O==1
    plot(corr_time_min,frac_18_O,'r-','Linewidth',2)
end
xlabel('time (min)')
ylabel('^1^8O fraction F')
axis([t_min_plot_scale t_max_plot_scale 0 1])
if plot_F==1
    yyaxis right
    set(gca,'YColor',[0 0 0])
    plot(corr_time_min,F,'k-','Linewidth',2)
    ylabel('F')
    axis([t_min_plot_scale t_max_plot_scale -inf inf])
    set(gcf, 'Position',  [100, 100, 800, 600],'color','w')
    if plot_frac_18_O==1
        legend('f_{18} (^*O_2 species)','f_{18} (^*O species)','F','Location','south','NumColumns',3)
    else
        legend('f_{18} (^*O_2 species)','F','Location','south','NumColumns',3)
    end
else
    if plot_frac_18_O==1
        legend('f_{18} (^*O_2 species)','f_{18} (^*O species)','F','Location','south','NumColumns',3)
    else
        legend('f_{18} (^*O_2 species)','Location','south','NumColumns',2)
    end
end
title('Summary of experimental conditions and titration')
annotation('textbox',[0.4 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')



if skip_plots==0
    
    choice5 = questdlg('Plot pressure and other species?', ...
        'Differentials', ...
        'yes','no','no');
    switch choice5
        case 'no'
        case 'yes'
            
            
            if exist('MS_pressure')==1
                plot_MS_pressure=figure('Name','Pressure','NumberTitle','off');
                semilogy(corr_time_min,MS_pressure,'color',[0 1.0 0],'Linewidth',2)
                hold on;
                ylabel('MS pressure')
                if max(MS_pressure)<=1e-8
                    axis([t_min_plot_scale t_max_plot_scale 1e-9 1e-8])
                elseif    max(MS_pressure)<=1e-7
                    axis([t_min_plot_scale t_max_plot_scale 1e-8 1e-7])
                elseif max(MS_pressure)<=1e-6
                    axis([t_min_plot_scale t_max_plot_scale 1e-8 1e-6])
                else
                    axis([t_min_plot_scale t_max_plot_scale 1e-7 1e-5])
                end
                legend('MS pressure','Location','east')
            end
            
            
            if exist('norm_17')+exist('norm_19')+exist('norm_20')+exist('norm_21')+exist('norm_22')~=0
                plot_frac_17_19_20_21_22=figure('Name','water & co','NumberTitle','off');
                hold on;
                dynamic_legend='';
                if exist('norm_17')==1 & isempty(norm_17)==0
                    dynamic_legend=append(dynamic_legend,'''I_{17}'',');
                    plot(corr_time_min,norm_17,'m-','Linewidth',2)
                end
                if exist('norm_19')==1 & isempty(norm_19)==0
                    dynamic_legend=append(dynamic_legend,'''I_{19}'',');
                    plot(corr_time_min,norm_19,'r-','Linewidth',2)
                end
                if exist('norm_20')==1 & isempty(norm_20)==0
                    dynamic_legend=append(dynamic_legend,'''I_{20}'',');
                    plot(corr_time_min,norm_20,'b-','Linewidth',2)
                end
                if exist('norm_21')==1 & isempty(norm_21)==0
                    dynamic_legend=append(dynamic_legend,'''I_{21}'',');
                    plot(corr_time_min,norm_21,'g-','Linewidth',2)
                end
                if exist('norm_22')==1 & isempty(norm_22)==0
                    dynamic_legend=append(dynamic_legend,'''I_{22}'',');
                    plot(corr_time_min,norm_22,'c-','Linewidth',2)
                end
                dynamic_legend = dynamic_legend(1:end-1);
                xlabel('time (min)')
                ylabel('I_i (normalised)')
                axis([t_min_plot_scale t_max_plot_scale 0 0.1])
                eval(sprintf('legend(%s,''Location'',''east'',''NumColumns'',2)',dynamic_legend));
            end
            
            %nitrogen? and argon
            if exist('norm_14')+ exist('norm_40')~=0
                plot_frac_14_40=figure('Name','nitrogen & co','NumberTitle','off');
                hold on;
                dynamic_legend='';
                if exist('norm_14')==1 & isempty(norm_14)==0
                    dynamic_legend=append(dynamic_legend,'''I_{14}'',');
                    plot(corr_time_min,norm_14,'g-','Linewidth',2)
                end
                if exist('norm_40')==1 & isempty(norm_40)==0
                    dynamic_legend=append(dynamic_legend,'''I_{40}'',');
                    plot(corr_time_min,norm_40,'k-','Linewidth',2)
                end
                dynamic_legend = dynamic_legend(1:end-1);
                xlabel('time (min)')
                ylabel('I_i (normalised)')
                axis([t_min_plot_scale t_max_plot_scale 0 inf])
                eval(sprintf('legend(%s,''Location'',''east'',''NumColumns'',2)',dynamic_legend));
            end
            
            
            
            %Carbonates
            if exist('norm_28')+exist('norm_44')+exist('norm_46')+exist('norm_48')~=0
                plot_frac_28_30_44_46_48=figure('Name','carbonates & co','NumberTitle','off');
                hold on;
                dynamic_legend='';
                if exist('norm_28')==1 & isempty(norm_28)==0
                    dynamic_legend=append(dynamic_legend,'''I_{28}'',');
                    plot(corr_time_min,norm_28,'k-','Linewidth',2)
                end
                if exist('norm_30')==1 & isempty(norm_30)==0
                    dynamic_legend=append(dynamic_legend,'''I_{30}'',');
                    plot(corr_time_min,norm_30,'c-','Linewidth',2)
                end
                if exist('norm_44')==1 & isempty(norm_44)==0
                    dynamic_legend=append(dynamic_legend,'''I_{44}'',');
                    plot(corr_time_min,norm_44,'b-','Linewidth',2)
                end
                if exist('norm_46')==1 & isempty(norm_46)==0
                    dynamic_legend=append(dynamic_legend,'''I_{46}'',');
                    plot(corr_time_min,norm_46,'m-','Linewidth',2)
                end
                if exist('norm_48')==1 & isempty(norm_48)==0
                    dynamic_legend=append(dynamic_legend,'''I_{48}'',');
                    plot(corr_time_min,norm_48,'r-','Linewidth',2)
                end
                dynamic_legend = dynamic_legend(1:end-1);
                xlabel('time (min)')
                ylabel('I_i (normalised)')
                axis([t_min_plot_scale t_max_plot_scale 0 0.1])
                eval(sprintf('legend(%s,''Location'',''east'',''NumColumns'',2)',dynamic_legend));
                
                
        
       if carbotruc==1 
                norm_CO2=norm_44+norm_46+norm_48;
                norm_norm_44=norm_44./norm_CO2;
                norm_norm_46=norm_46./norm_CO2;
                norm_norm_48=norm_48./norm_CO2;
                frac_18_CO2=norm_norm_48+norm_norm_46/2;                
                plot_norm_CO2=figure('Name','CO2 evolution','NumberTitle','off');
                hold on;                
                plot(corr_time_min,norm_norm_44,'b-','Linewidth',2)
                plot(corr_time_min,norm_norm_46,'m-','Linewidth',2)
                plot(corr_time_min,norm_norm_48,'r-','Linewidth',2)
                plot(corr_time_min,frac_18_CO2,'k-','Linewidth',2)
                xlabel('time (min)')
                ylabel('CO_2 isotopologues fraction')
                axis([t_min_plot_scale t_max_plot_scale 0 1])
                legend('f_{44}','f_{46}','f_{48}','f_{18} (^C*O2 species)','Location','east')
       end                
                
                
            end
            
    end
    
    
    choice6 = questdlg('Visualise how far is the gas phase to the thermodynamic equilibrium?', ...
        'Show statistical equilibrium', ...
        'yes','no','yes');
    switch choice6
        case 'no'
        case 'yes'
            % fraction of each oxygen Isotopologue at statistical equilibrium as a function of the 18O fraction
            th_36=frac_18_O2.^2;
            th_32=(1-frac_18_O2).^2;
            th_34=2*(1-frac_18_O2).*frac_18_O2;
            
            plot_thvsexp=figure('Name','Towards thermodynamic equilibrium','NumberTitle','off');
            hold on;
            plot(corr_time_min,norm_32,'b-','Linewidth',2)
            plot(corr_time_min,norm_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
            plot(corr_time_min,norm_36,'r-','Linewidth',2)
            plot(corr_time_min,th_32,'b--','Linewidth',0.5)
            plot(corr_time_min,th_34,'--','color',[1.0  0.0  1.0],'Linewidth',0.5)
            plot(corr_time_min,th_36,'r--','Linewidth',0.5)
            xlabel('time (min)')
            ylabel('Isotopologue fraction f_i')
            axis([0 t_max_plot_scale 0 1])
            legend('f_3_2','f_3_4','f_3_6','th_3_2','th_3_4','th_3_6','Location','east')
            
    end
    
    
    
    
    choice7 = questdlg('Plot the differrentials?', ...
        'Differentials', ...
        'yes','no','no');
    switch choice7
        case 'no'
        case 'yes'
            % difference between the statistical equilibrium and the real
            % fraction of Isotopologues (should be proportinnal to the driving
            % force?
            th_36=frac_18_O2.^2;
            th_32=(1-frac_18_O2).^2;
            th_34=2*(1-frac_18_O2).*frac_18_O2;
            delta_32=norm_32-th_32;
            delta_34=norm_34-th_34;
            delta_36=norm_36-th_36;
            
            
            
            %%calculation of differentials with a step
            for i=diff_step+1:size(corr_time_min)-diff_step
                diff_32(i,1)=(norm_32(i+diff_step)-norm_32(i-diff_step))/(corr_time_min(i+diff_step)-corr_time_min(i-diff_step));
                diff_34(i,1)=(norm_34(i+diff_step)-norm_34(i-diff_step))/(corr_time_min(i+diff_step)-corr_time_min(i-diff_step));
                diff_36(i,1)=(norm_36(i+diff_step)-norm_36(i-diff_step))/(corr_time_min(i+diff_step)-corr_time_min(i-diff_step));
            end
            diff_32=diff_32(diff_step+1:end-1,1);
            diff_34=diff_34(diff_step+1:end-1,1);
            diff_36=diff_36(diff_step+1:end-1,1);
            
            corr_time_diff_min=corr_time_min(diff_step+1:end-diff_step-1,1);
            
            
            plot_deltavsdiff=figure('Name','Differential vs gap between real-time and thq. equilibrium','NumberTitle','off');
            hold on;
            plot(corr_time_diff_min,diff_32,'b-','Linewidth',2)
            plot(corr_time_diff_min,diff_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
            plot(corr_time_diff_min,diff_36,'r-','Linewidth',2)
            xlabel('time (min)')
            ylabel('d(f_i) / dt (min^-^1)')
            axis([0 t_max_plot_scale -inf inf])
            yyaxis right
            set(gca,'YColor',[0 0 0])
            plot(corr_time_min,delta_32,'b-.','Linewidth',0.5)
            plot(corr_time_min,delta_34,'--','color',[1.0  0.0  1.0],'Linewidth',0.5)
            plot(corr_time_min,delta_36,'r--','Linewidth',0.5)
            ylabel('delta f_i')
            axis([0 t_max_plot_scale -inf inf])
            legend('d(f_3_2) / dt','d(f_3_4) / dt','d(f_3_6) / dt','delta_3_2','delta_3_4','delta_3_6','Location','east')
    end
    
end % end of multiple figures at the begining


%% Surface kinetics:
norm_frac_18_O2=(frac_18_O2-frac_18_O2_t0)/(frac_18_O2_inf-frac_18_O2_t0);
message_ks=questdlg('Extract the surface exchange kinetics?', ...
    'Next step', ...
    'Yes','No','Yes');
switch message_ks
    case 'No'
    case 'Yes'
                
        figure(plot_frac_18_O2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mass_transport_correction (experimental) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mass_transport_correction==1
            %get a delay
            extract_tau_data={corr_time,norm_frac_18_O2,1,size(corr_time)};
            f_mod=3;
            fitting_param=[0 0 f_mod];
            guess=[600 120];
            tau_and_delay_calc=fminsearch('exponential_decay_extract_tau_and_delay',guess,options_fit);
            tau_calc=tau_and_delay_calc(1);
            delay=tau_and_delay_calc(2);
            % delay=1.292696328178378e+02;;
            tau_mass=60;
            mass_transport_correction=1./(1-exp(-(corr_time-delay)/tau_mass));
            norm_frac_18_O2_corr=norm_frac_18_O2.*mass_transport_correction;
            %                  plot_test=figure('Name','test','NumberTitle','off');
            %                  hold on;
            %                          plot(corr_time,norm_frac_18_O2)
            %               plot(corr_time,norm_frac_18_O2_corr)
            extract_tau_data={corr_time,norm_frac_18_O2_corr,15,size(corr_time)};
            guess=[tau_calc delay+20];
            tau_and_delay_calc2=fminsearch('exponential_decay_extract_tau_and_delay',guess,options_fit);
            tau_calc2=tau_and_delay_calc2(1);
            tau_calc2_min=tau_calc2/60;
            r_s2=A*B/(tau_calc2*S*(A+B));%mol/m^2/s
            ks2=r_s2/c_o;
            delay2=tau_and_delay_calc2(2);
            norm_frac_18_O2_calc2=1-exp(-1/tau_calc2.*(corr_time-delay2));
            %replace by 1 for time values lower than the delay:
            for i=1:size(corr_time)
                if corr_time(i)<delay2
                    norm_frac_18_O2_calc2(i)=0;
                end
            end
            frac_18_O2_calc2=norm_frac_18_O2_calc2*(frac_18_O2_inf-frac_18_O2_t0)+frac_18_O2_t0;
            frac_18_O2_corr=norm_frac_18_O2_corr*(frac_18_O2_inf-frac_18_O2_t0)+frac_18_O2_t0;
            
            content_annotation={sprintf('delay correction = %0.3g s',delay2),sprintf('\\Re_0 = %0.2e mol/m^2/s',r_s2),sprintf('k_s = %0.3g m/s',ks2)};%sprintf('\\tau = %0.4g s',tau_calc),
            plot_fit_frac_18_O2_mass=figure('Name','18O fraction: fitting of surface exchange rate','NumberTitle','off');
            plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
            hold on;
            plot(corr_time_min,frac_18_O2_corr,'k.','Linewidth',2)
            plot(corr_time_min,frac_18_O2_calc2,'r-','Linewidth',2)
            xlabel('time (min)')
            ylabel('^1^8O fraction')
            axis([0 t_max_plot_scale 0 1])
            legend('f_{18} (^*O_2)','f_{18} (^*O_2)(corr)','best fit','Location','east')
            annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
            frac_18_O2=frac_18_O2_corr;
            norm_frac_18_O2=(frac_18_O2-frac_18_O2_t0)/(frac_18_O2_inf-frac_18_O2_t0);
            %norm_frac_18_O2=norm_frac_18_O2_corr;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% END mass_transport_correction (experimental) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% fitting a diffusion limited process (Mizusaki 1984) %%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        choice = questdlg('Fit a diffusion limited process (Mizusaki 1984)?', ...
            'Mizusaki fitting', ...
            'yes','no','yes');
        switch choice
            case 'no'
            case 'yes'
                q=Mizusaki_D_roots([alpha nroots]);
                custom_message='Choose the fitting window:';
                custom_indexes=[15 size(corr_time_min)];
                [low_index high_index]=restrain_time_range(corr_time_min);
                fitting_data={corr_time,norm_frac_18_O2,A,B,radius,low_index,high_index,q};
                % screening D over 30 orders of magnitude :)
                logDscr=1:30;
                D_screening=10.^(-logDscr);
                for i = 1:size(D_screening,2)
                    Mizu_scr(i,1)=Mizusaki_D([D_screening(i) 120]);
                end
                [Mizu,I] = min(Mizu_scr);
                D_guess=1*10^-I;
                guess=[D_guess 120];% 5e-5]; prefer to keep the radius constant for now
                
                %fitting
                D_and_delay_calc=fminsearch('Mizusaki_D',guess,options_fit);
                D_calc=D_and_delay_calc(1);
                delay=D_and_delay_calc(2);
                %%%%%%%%%%%%%%radius_calc=radius;%D_and_delay_calc(3); prefer to keep the radius constant for now
                %%%%%%%%%%%%% %radius_calc_um=radius_calc*1E6;
                norm_frac_18_O2_calc_Mizu_D=zeros(size(corr_time));
                for n=1:1:nroots
                    %norm_frac_18_O2_calc_Mizu_D=norm_frac_18_O2_calc_Mizu_D+(6*alpha*(alpha+1)*exp(-D_calc*q(n)^2.*(corr_time-delay)/radius_calc^2))/(9+9*alpha+q(n)^2*alpha^2);
                    norm_frac_18_O2_calc_Mizu_D=norm_frac_18_O2_calc_Mizu_D+(6*alpha*(alpha+1)*exp(-D_calc*q(n)^2.*(corr_time-delay)/radius^2))/(9+9*alpha+q(n)^2*alpha^2);
                end
                norm_frac_18_O2_calc_Mizu_D=1-norm_frac_18_O2_calc_Mizu_D;
                for i=1:size(corr_time)
                    if corr_time(i)<delay
                        norm_frac_18_O2_calc_Mizu_D(i)=0;
                    end
                end
                frac_18_O2_calc_Mizu_D=norm_frac_18_O2_calc_Mizu_D*(frac_18_O2_inf-frac_18_O2_t0)+frac_18_O2_t0;
                
                
                content_annotation={sprintf('delay correction = %0.3g s',delay),sprintf('D = %0.3g m^2/s',D_calc),sprintf('%s = %0.3g \\mum',radius_display,radius_um),sprintf('%s = %0.3g m^2/g',SSA_display,SSA),sprintf('\\alpha = %0.3g',alpha),sprintf('n. roots = %d ',nroots)};
                plot_fit_Mizusaki_D_limited=figure('Name','18O fraction: fitting of limiting D','NumberTitle','off');
                hold on;
                plot(corr_time_min,frac_18_O2,'k.','Linewidth',2)
                plot(corr_time_min(low_index:high_index),frac_18_O2_calc_Mizu_D(low_index:high_index),'r-','Linewidth',1)
                plot(corr_time_min,frac_18_O2_calc_Mizu_D,'r:','Linewidth',1)
                xlabel('time (min)')
                ylabel('^1^8O fraction')
                axis([0 t_max_plot_scale 0 1])
                title('Mizusaki''s(Edwards) approach (diffusion limitation)')
                
                legend('f_{18} (^*O_2)','best fit','Location','east')
                annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                
                % for export later:
                D_Mizu_D=D_calc;
                delay_Mizu_D=delay;
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% fitting D and k (Mizusaki 1984) %%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        choice = questdlg('Fit D and k (Mizusaki 1984)?', ...
            'Mizusaki fitting', ...
            'yes','no','yes');
        switch choice
            case 'no'
            case 'yes'
                
                q=Mizusaki_D_roots([alpha nroots]);
                custom_message='Choose the fitting window:';
                custom_indexes=[15 size(corr_time_min)];
                [low_index high_index]=restrain_time_range(corr_time_min);
                fitting_data={corr_time,norm_frac_18_O2,A,B,radius,low_index,high_index,q};
                
                
                %setting optimal starting values of delay_guess D_guess and k_guess
                if exist('delay')==1
                    delay_guess=delay;
                else
                    delay_guess=120;
                end
                if exist('radius')==1
                    radius_guess=radius;
                else
                    radius_guess=5e-5;
                end
                if exist('D_calc')==1
                    D_guess=D_calc;
                else
                    % screening D over 30 orders of magnitude :)
                    logDscr=1:30;
                    D_screening=10.^(-logDscr);
                    for i = 1:size(D_screening,2)
                        D_scr(i,1)=Mizusaki_D(D_screening(i));
                    end
                    [Mizu,I] = min(D_scr);
                    D_guess=1*10^-I;
                end
                % screening k over 30 orders of magnitude :)
                f_mod=4;
                fitting_param=[D_guess 0 delay_guess radius_guess f_mod];
                logkscr=1:30;
                k_screening=10.^(-logkscr);
                for i = 1:size(k_screening,2)
                    k_scr(i,1)=Mizusaki_kD([k_screening(i)]);
                end
                [Mizu,I] = min(k_scr);
                k_guess=1*10^-(I);
                % D_guess=1.18666E-13;%Benson650
                % D_guess=6.45747E-15;
                %k_guess=6.33186E-10;
                %D_guess=4.06448E-13;% Benson700
                %k_guess=1.45161E-8;% Benson700
                %D_guess=4.06448E-13;
                % k_guess=2.15e-11;
                %delay_guess=300;
%                 D_guess=1.1586E-18;
%                 k_guess=3.4962e-13;
%                 delay_guess=120;
                
                %setting fitting_param variable
                f_mod=0;
                
                if fitD==1
                    f_mod=bitset(f_mod,4);
                else
                      input_D = inputdlg({'D (m^2/s)'},...
                '',...
                1);
            input_D=str2double(input_D);
            D_guess=input_D(1);            
                end
                
                if fitk==1
                    f_mod=bitset(f_mod,3);
                else
                      input_k = inputdlg({'k (m/s)'},...
                '',...
                1);
            input_k=str2double(input_k);
            k_guess=input_k(1);
                end
                
                if fitdelay==1
                    f_mod=bitset(f_mod,2);
                else
                     choice = questdlg('Setting the delay:', ...
            'Input requested', ...
            'Set visually','Type','Set visually');
        switch choice
            case 'Set visually'
                    message = sprintf('Select the delay');
                    uiwait(msgbox(message,'Notification'));
                    delay_coord=ginput(1);                
             delay_guess=delay_coord(1);
            case 'Type'
                    input_delay = inputdlg({'Set the delay (s):'},...
                'Delay',...
                1);
            input_delay=str2double(input_delay);
            delay_guess=input_delay(1);      
        end     
%                       input_delay = inputdlg({'delay (s)'},...
%                 '',...
%                 1);
%             input_delay=str2double(input_delay);
%             delay_guess=input_D(1);
            end
                
                if fitradius==1
                    f_mod=bitset(f_mod,1);
                else
                choice_radius = questdlg('Set the radius from:', ...
        'Choice of radius', ...
        'the log','manual input','the log');
    switch choice_radius
        case 'the log'
         case 'manual input'
         input_radius = inputdlg({'particle radius (um)'},...
                '',...
                1);
            input_radius=str2double(input_radius);
            radius_guess=input_radius(1)/1E6;
    end      
                end
                
                fitting_param=[D_guess k_guess delay_guess radius_guess f_mod];
                
                % setting the guess for fminsearch
                n=1;
                clear('guess');
                if bitget(f_mod,4)==1
                    guess(n)=D_guess;
                    n=n+1;
                end
                if bitget(f_mod,3)==1
                    guess(n)=k_guess;
                    n=n+1;
                end
                if bitget(f_mod,2)==1
                    guess(n)=delay_guess;
                    n=n+1;
                end
                if bitget(f_mod,1)==1
                    guess(n)=radius_guess;
                end
                
                %fitting
                D_k_delay_and_radius_calc=fminsearch('Mizusaki_kD',guess,options_fit);
                n=1;
                if bitget(f_mod,4)==0
                    D_calc=D_guess;
                else
                    D_calc=D_k_delay_and_radius_calc(n);
                    n=n+1;
                end
                if bitget(f_mod,3)==0
                    k_calc=k_guess;
                else
                    k_calc=D_k_delay_and_radius_calc(n);
                    n=n+1;
                end
                if bitget(f_mod,2)==0
                    delay=delay_guess;
                else
                    delay=D_k_delay_and_radius_calc(n);
                    n=n+1;
                end
                if bitget(f_mod,1)==0
                    radius_calc=radius_guess;
                else
                    radius_calc=D_k_delay_and_radius_calc(n);
                end
                radius_calc_um=radius_calc*1E6;
                L_C_Mizu=D_calc/k_calc;
                L_C_Mizu_um=L_C_Mizu*1e6;
                
                kaD=k_calc*radius_calc/D_calc;
                q=Mizusaki_kD_roots([k_calc D_calc alpha radius_calc nroots]);
                norm_frac_18_O2_calc_Mizu_kD=zeros(size(corr_time));
                for n=1:1:size(q,1);%nroots
                    norm_frac_18_O2_calc_Mizu_kD=norm_frac_18_O2_calc_Mizu_kD+(6*alpha*kaD^2*(alpha+1)*exp(-D_calc*q(n)^2.*(corr_time-delay)/radius_calc^2))/(alpha^2*q(n)^4+alpha*kaD*(alpha*(kaD-1)-6)*q(n)^2+9*(1+alpha)*kaD^2);
                end
                norm_frac_18_O2_calc_Mizu_kD=1-norm_frac_18_O2_calc_Mizu_kD;
                
                for i=1:size(corr_time)
                    if corr_time(i)<delay
                        norm_frac_18_O2_calc_Mizu_kD(i)=0;
                    end
                end
                
                frac_18_O2_calc_Mizu_kD=norm_frac_18_O2_calc_Mizu_kD*(frac_18_O2_inf-frac_18_O2_t0)+frac_18_O2_t0;
                %%0.1f
                content_annotation={sprintf('delay correction = %0.3g s',delay),sprintf('k = %0.3g m/s',k_calc),sprintf('D = %0.3g m^2/s',D_calc),sprintf('L_c = %0.3g \\mum',L_C_Mizu_um),sprintf('%s = %0.3g \\mum',radius_display,radius_calc_um),sprintf('%s = %0.3g m^2/g',SSA_display,SSA),sprintf('\\alpha = %0.3g',alpha),sprintf('n. roots = %d ',nroots)};
                plot_fit_Mizusaki_D_an_k=figure('Name','18O fraction: fitting of D and k','NumberTitle','off');
                hold on;
                plot(corr_time_min,frac_18_O2,'k.','Linewidth',2)
                plot(corr_time_min(low_index:high_index),frac_18_O2_calc_Mizu_kD(low_index:high_index),'r-','Linewidth',1)
                plot(corr_time_min,frac_18_O2_calc_Mizu_kD,'r:','Linewidth',1)
                xlabel('time (min)')
                ylabel('^1^8O fraction')
                axis([0 t_max_plot_scale 0 1])
                title('Mizusaki''s(Edwards) approach (mixed regime)')
                legend('f_{18} (^*O_2)','best fit','Location','east')
                annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                
                % for export later:
                k_Mizu_kD=k_calc;
                D_Mizu_kD=D_calc;
                delay_Mizu_kD=delay;
                radius_Mizu_kD=radius_calc_um;
                
        end
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% fitting the exchange and dissociation rate (Den Otter 2001) %%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        choice = questdlg('Fit to the 2 steps model (Den Otter 2001)?', ...
            '2 steps model fitting', ...
            'yes','no','yes');
        switch choice
            case 'no'
            case 'yes'
                %Fit the 18O fraction as an exponential decay, get the
                %characteristic time tau and r_s, k_s
                norm_frac_18_O2=(frac_18_O2-frac_18_O2_t0)/(frac_18_O2_inf-frac_18_O2_t0);
               %could recall plot_norm323436shift instead but the time range might
                %not be good (after constraint):
                %figure(plot_norm323436shift)                
                plot_2steps_b4_fitting=figure('Name','fitting with the 2 steps model','NumberTitle','off');
                hold on;
                plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
                plot(corr_time_min,norm_32,'b-','Linewidth',2)
                plot(corr_time_min,norm_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
                plot(corr_time_min,norm_36,'r-','Linewidth',2)
                xlabel('time (min)')
                ylabel('Isotopologue fraction f_i')
                axis([0 t_max_plot_scale 0 1])
                title('Den Otter approach (2001)')
                legend('^{18}O fraction','f_{32}','f_{34}','f_{36}','Location','north','NumColumns',3)
                custom_message='Optimise the fitting window for the 18O fraction:';
                custom_indexes=[10 size(corr_time_min)];
                [low_index high_index]=restrain_time_range(corr_time_min);
                extract_tau_data={corr_time,norm_frac_18_O2,low_index,high_index};
                %if delay_correction==1
%                     guess=[600 120];
%                     tau_and_delay_calc=fminsearch('exponential_decay_extract_tau_and_delay',guess,options_fit);
%                     tau_calc=tau_and_delay_calc(1);
%                     delay=tau_and_delay_calc(2);
tau_guess=600;
delay_guess=120;
f_mod=0;
                if fit_tau==1
                    f_mod=bitset(f_mod,2);
                else
                      input_tau = inputdlg({'good luck with guessing tau! (s):'},...
                '',...
                1);
            input_tau=str2double(input_tau);
            tau_guess=input_tau(1);
                end
                if delay_correction==1
                    f_mod=bitset(f_mod,1);
                else
                    choice = questdlg('Setting the delay:', ...
            'Den Otter fitting', ...
            'Set visually','Type','Set visually');
        switch choice
            case 'Set visually'
                    message = sprintf('Select the delay');
                    uiwait(msgbox(message,'Notification'));
                    delay_coord=ginput(1);                
             delay_guess=delay_coord(1);
            case 'Type'
                    input_delay = inputdlg({'Set the delay (s):'},...
                'Delay',...
                1);
            input_delay=str2double(input_delay);
            delay_guess=input_delay(1);      
        end
                end
                fitting_param=[tau_guess delay_guess f_mod];
                n=1;
                clear('guess');
               if bitget(f_mod,2)==1
                    guess(n)=tau_guess;
                    n=n+1;
                end
                if bitget(f_mod,1)==1
                    guess(n)=delay_guess;
                end
                    tau_and_delay_calc=fminsearch('exponential_decay_extract_tau_and_delay',guess,options_fit);
                   % guess=[600 120];
                   % tau_calc=tau_and_delay_calc(1);
                   % delay=tau_and_delay_calc(2);
                    n=1;
              if bitget(f_mod,2)==0
                    tau_calc=tau_guess;
                else
                    tau_calc=tau_and_delay_calc(n);
                    n=n+1;
                end
                if bitget(f_mod,1)==0
                    delay_calc=delay_guess;
                else
                    delay_calc=tau_and_delay_calc(n);
                end
                %else
                %    tau_calc=fminsearch('exponential_decay_extract_tau',600,options_fit);
                %    delay=0;
                %end
                
                %%%%%%%%%%%%%%%%%%%%%%
                %%%% EXPERIMENTAL %%%%
                %%%%%%%%%%%%%%%%%%%%%%
                if  background_correction==1
                    bg_corr_data={corr_time,raw_intensity32,raw_intensity34,raw_intensity36,frac_18_O2_t0,frac_18_O2_inf,tau_calc,delay_calc,low_index,high_index};
                    bg_guess=[4.37E-13 3.23E-13 1.90841E-13];
                    bg_calc=fminsearch('bg_corr',bg_guess,options_fit);
                    bg_intensity32=bg_calc(1);
                    bg_intensity34=bg_calc(2);
                    bg_intensity36=bg_calc(3);
                    corr_intensity32=raw_intensity32-bg_intensity32;
                    corr_intensity34=raw_intensity34-bg_intensity34;
                    corr_intensity36=raw_intensity36-bg_intensity36;
                    sum_intensityO2=corr_intensity32+corr_intensity34+corr_intensity36;
                    norm_32=corr_intensity32./sum_intensityO2;
                    norm_34=corr_intensity34./sum_intensityO2;
                    norm_36=corr_intensity36./sum_intensityO2;
                    frac_18_O2=norm_36+1/2*norm_34;
                    norm_frac_18_O2=(frac_18_O2-frac_18_O2_t0)/(frac_18_O2_inf-frac_18_O2_t0);
                    norm_calc=1-exp(-1/tau_calc.*(corr_time-delay_calc));
                    
                    tau_calc_min=tau_calc/60;
%                     plot_tau=figure('Name','Characteristic time','NumberTitle','off');
%                     hold on;
%                     plot(corr_time_min,norm_frac_18_O2,'b-','Linewidth',2)
%                     plot(corr_time_min(low_index:high_index),norm_calc(low_index:high_index),'r-','Linewidth',2)
%                     plot(corr_time_min,norm_calc,'r--','Linewidth',2)
%                     xlabel('time (min)')
%                     ylabel('Normalised ^{18}O fraction')
%                     axis([0 t_max_plot_scale 0 1])
%                     title('Background correction')
%                     legend('normalised f_{18} (^*O_2)','best fit','Location','south')
%                     annotation('textbox',[0.35 0.35 0.3 0.3],'String',{sprintf('\\tau = %0.4g min',tau_calc_min),sprintf('delay correction = %0.3g s',delay)},'FitBoxToText','on')
                end
                %%%%%%%%%%%%%%%%%%%%%%
                %%% /EXPERIMENTAL/ %%%
                %%%%%%%%%%%%%%%%%%%%%%                
                
                tau_calc_min=tau_calc/60;
                delay_calc_min=delay_calc/60;
                r_s_DO=A*B/(tau_calc*S*(A+B));%mol/m^2/s
                ks_DO=r_s_DO/c_o;
                norm_calc=1-exp(-1/tau_calc.*(corr_time-delay_calc));
                
                if plot_tau==1
                plot_charac_time=figure('Name','Characteristic time','NumberTitle','off');
                hold on;
                plot(corr_time_min,norm_frac_18_O2,'b-','Linewidth',2)
                plot(corr_time_min(low_index:high_index),norm_calc(low_index:high_index),'r-','Linewidth',2)
                plot(corr_time_min,norm_calc,'r--','Linewidth',2)
                xlabel('time (min)')
                ylabel('Normalised ^{18}O fraction')
                axis([0 t_max_plot_scale 0 1])
                if background_correction==0
                title('Den Otter''s approach (surface exchange limitation)')
                elseif background_correction==1
                title({'Den Otter''s approach (surface exchange limitation)','WARNING! background correction is [on]'})   
                end
                legend('normalised f_{18} (^*O_2)','best fit','Location','south')
                annotation('textbox',[0.35 0.35 0.3 0.3],'String',{sprintf('\\tau = %0.4g min',tau_calc_min),sprintf('delay correction = %0.3g s',delay_calc)},'FitBoxToText','on')
                end
                
                norm_frac_18_O2_calc_DenOtter=1-exp(-1/tau_calc.*(corr_time-delay_calc)); %%% = norm_calc !!!!
                %replace by 1 for time values lower than the delay:
                for i=1:size(corr_time)
                    if corr_time(i)<delay_calc
                        norm_frac_18_O2_calc(i)=0;
                    end
                end
                frac_18_O2_calc_DenOtter=norm_frac_18_O2_calc_DenOtter*(frac_18_O2_inf-frac_18_O2_t0)+frac_18_O2_t0;
                
                content_annotation={sprintf('\\Re_0 = %0.2e mol/m^2/s',r_s_DO),sprintf('k_s = %0.3g m/s',ks_DO),sprintf('\\tau = %0.4g min',tau_calc_min),sprintf('delay correction = %0.3g s',delay_calc)};%,
                plot_fit_frac_18_O2=figure('Name','18O fraction: fitting of surface exchange rate','NumberTitle','off');
                plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
                hold on;
                plot(corr_time_min(low_index:high_index),frac_18_O2_calc_DenOtter(low_index:high_index),'r-','Linewidth',2)
                plot(corr_time_min,frac_18_O2_calc_DenOtter,'r--','Linewidth',2)
                xlabel('time (min)')
                ylabel('^1^8O fraction')
                axis([0 t_max_plot_scale 0 1])
                if background_correction==0
                title('Den Otter''s approach (surface exchange limitation)')
                elseif background_correction==1
                title('Background correction')   
                end
                legend('f_{18} (^*O_2)','best fit','Location','east')
                annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                
                
                % % % % % % % %             norm_frac_18_CO2=(frac_18_CO2-0)/(0.1737-0);
                % % % % % % % %
                % % % % % % % %             extract_tau_data={corr_time-delay_calc,norm_frac_18_CO2,15,size(corr_time)};
                % % % % % % % %             tau_calc2=fminsearch('exponential_decay_extract_tau',600)
                % % % % % % % %             tau_calc2_min=tau_calc2/60;
                % % % % % % % %
                % % % % % % % %              norm_calc=1-exp(-1/tau_calc2.*(corr_time-delay_calc));
                % % % % % % % %
                % % % % % % % %
                % % % % % % % %    plot_tau_CO2=figure('Name','CO2 evolution','NumberTitle','off');
                % % % % % % % % hold on;
                % % % % % % % % plot(corr_time_min-delay_calc/60,norm_frac_18_CO2,'b-','Linewidth',2)
                % % % % % % % % plot(corr_time_min-delay_calc/60,norm_calc,'r-','Linewidth',2)
                % % % % % % % %
                % % % % % % % %
                % % % % % % % %             xlabel('time (min)')
                % % % % % % % %             ylabel('dirty normalised f18 CO_2')
                % % % % % % % %             axis([0 t_max_plot_scale 0 1])
                % % % % % % % %             legend('f_{18}','Location','east')
                % % % % % % % %
                % % % % % % % %              annotation('textbox',[0.35 0.35 0.3 0.3],'String',sprintf('\\tau = %0.4g min',tau_calc2_min),'FitBoxToText','on')
                
                
                
                figure(plot_2steps_b4_fitting);% bring the figure back to foreground       
                custom_message='Optimise the fitting window for the 16O18O isotopologue:';
                custom_indexes=[1 size(corr_time_min)];
                [low_index high_index]=restrain_time_range(corr_time_min);
                close(plot_2steps_b4_fitting);
                %%fitting 1 parameter only (p)
                p_guess=1;
                %corr_time_minus_delay=corr_time-delay_calc;
                %two_steps_model_data={corr_time_minus_delay,norm_32,norm_34,norm_36,A,B,low_index,high_index};
                two_steps_model_data={corr_time,norm_32,norm_34,norm_36,A,B,low_index,high_index};
                f_mod=8;
                fitting_param=[0 0 delay_calc tau_calc f_mod]; 
                
                p_calc=fminsearch('Den_Otter_two_steps',p_guess,options_fit);                                
                p1=p_calc;
                p2=p_calc;
                tau_1= tau_calc;
                f_36_t0=norm_36_t0;
                f_18g_t0=frac_18_O2_t0;%f_36_t0+f_34_t0/2;
                f_18b_t0=f18_sample_t0;%0.00204;
                f_18_inf=frac_18_O2_inf;%f_36_inf+f_34_inf/2;%(A*f_18g_t0+B*f_18b_t0)/(A+B);
                beta_0=f_18_inf^2;
                beta_1=f_18_inf*((f_18b_t0-f_18g_t0)*(p1+p2)+2*(f_18g_t0-f_18_inf));
                beta_2=((f_18b_t0-f_18g_t0)*p1+f_18g_t0-f_18_inf)*((f_18b_t0-f_18g_t0)*p2+f_18g_t0-f_18_inf);
                tau_2=tau_1*(1-beta_1/(2*f_18_inf*(f_18g_t0-f_18_inf)));
                epsilon_1=2*f_18_inf*(f_18g_t0-f_18_inf);
                epsilon_2=beta_2*tau_1/(tau_1-2*tau_2);
                f_36=f_18_inf^2+epsilon_1*exp(-(corr_time-delay_calc)/tau_1)+epsilon_2*exp(-2*(corr_time-delay_calc)/tau_1)+(f_36_t0-f_18_inf^2-epsilon_1-epsilon_2)*exp(-(corr_time-delay_calc)/tau_2);
                f_18g=f_18_inf+(f_18g_t0-f_18_inf)*exp(-(corr_time-delay_calc)/tau_1);
                f_18b=f_18_inf+(f_18b_t0-f_18_inf)*exp(-(corr_time-delay_calc)/tau_1);
                f_34=2*(f_18g-f_36);
                f_32=1-f_34-f_36;
                
                for i=1:size(corr_time)
                    if corr_time(i)<delay_calc
                        f_36(i)=norm_36_t0;
                        f_34(i)=norm_34(1);
                        f_32(i)=norm_32(1);
                    end
                end
                
                kdis_DO=ks_DO*2/(p1+p2);%kdis according to Den Otter
                
                content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('Date : %s',experiment_date),sprintf('k_s = %0.3g m/s',ks_DO),sprintf('k_{dis} = %0.3g m/s',kdis_DO),sprintf('\\tau = %0.4g min , delay = %0.3g s',tau_calc_min,delay_calc),sprintf('p_1 = p_2 = %0.3g',p1)};
                %{sprintf('%s (fit)',raw_data_filename),sprintf('Date & time:%s - %s',experiment_date,experiment_time),sprintf('\\Re_s = %0.3g mol/m^2/s',r_s),sprintf('\\Re_{dis} = %0.3g mol/m^2/s',r_dis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis} = %0.3g m/s',kdis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg   |   \\rho = %0.3g g/cm^3',m_mg,rho_gcc),sprintf('SSA = %0.3g m^2/g',SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('M_{calc} = %0.4g g/mol   |   \\nu_{calc} = %0.3g',M_calc,nu_calc)};
                plot_2steps_fitting=figure('Name','fitting with the 2 steps model','NumberTitle','off');
                hold on;
                plot(corr_time_min,norm_32,'b-','Linewidth',2)
                plot(corr_time_min(low_index:high_index),f_32(low_index:high_index),'b--','Linewidth',1)
                plot(corr_time_min,norm_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
                plot(corr_time_min(low_index:high_index),f_34(low_index:high_index),'--','color',[1.0  0.0  1.0],'Linewidth',1)
                plot(corr_time_min,norm_36,'r-','Linewidth',2)
                plot(corr_time_min(low_index:high_index),f_36(low_index:high_index),'r--','Linewidth',1)
                
                plot(corr_time_min(1:low_index),f_32(1:low_index),'b:','Linewidth',1)
                plot(corr_time_min(1:low_index),f_34(1:low_index),':','color',[1.0  0.0  1.0],'Linewidth',1)
                plot(corr_time_min(1:low_index),f_36(1:low_index),'r:','Linewidth',1)
                plot(corr_time_min(high_index:end),f_32(high_index:end),'b:','Linewidth',1)
                plot(corr_time_min(high_index:end),f_34(high_index:end),':','color',[1.0  0.0  1.0],'Linewidth',1)
                plot(corr_time_min(high_index:end),f_36(high_index:end),'r:','Linewidth',1)
%                 
%                 plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
%                 plot(corr_time_min(low_index:high_index),frac_18_O2_calc_DenOtter(low_index:high_index),'c-','Linewidth',2)
%                 plot(corr_time_min,frac_18_O2_calc_DenOtter,'c--','Linewidth',2)
                
                title('Den Otter''s approach (surface exchange limitation)')
                xlabel('time (min)')
                ylabel('Isotopologue fraction f_i')
                axis([0 t_max_plot_scale 0 1])
                legend('f_{32}','f^{calc}_{32}','f_{34}','f^{calc}_{34}','f_{36}','f^{calc}_{36}','Location','north','NumColumns',3) %'f_{18} (^*O_2)','best fit'
                annotation('textbox',[0.35 0.35 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                
                
                
                % for export later:
                delay_DO=delay_calc;
                p_DO=p1;
                
                
                
                
                %% Independant fitting of p1, p2,tau and a delay 
                if advanced_DO_fitting==1
                    
                choice18 = questdlg('Further fitting to the 2 steps model ? (p1, p2, a delay and tau can be fitted as chosen in the parameters) !![testing]!!', ...
                    '2 steps model fitting [beta]', ...
                    'yes','no','no');
                switch choice18
                    case 'no'
                    case 'yes'
                        
                        
                        choice = questdlg('Optimise the fitting window?', ...
            'Input requested', ...
            'yes','same as before','yes');
        switch choice
            case 'same as before'
            case 'yes'            
                custom_message='Choose the fitting window:';
                custom_indexes=[1 size(corr_time_min)];
                [low_index high_index]=restrain_time_range(corr_time_min);
        end          
                        
         f_mod=0;
                if fit_p1==1
                    f_mod=bitset(f_mod,4);
                    p1_guess=p_calc;
                else                    
                    choice = questdlg('Setting p1:', ...
            'Input requested', ...
            'same as before','type','same as before');
        switch choice
            case 'same as before'
                  p1_guess=p_calc;
            case 'type'       
                p1_guess=str2double(inputdlg({'Set p1:'},'Setting p1'));
        end
                end
                
                if fit_p2==1
                    f_mod=bitset(f_mod,3);
                    p2_guess=p_calc;
                else
                   choice = questdlg('Setting p2:', ...
            'Input requested', ...
            'same as before','type','same as before');
        switch choice
            case 'same as before'
                  p2_guess=p_calc;
            case 'type' 
                      p2_guess=str2double(inputdlg({'Set p2:'},'Setting p2'));
                end     
                end
                
                if fit_delay_DO==1
                    f_mod=bitset(f_mod,2);
                    delay_guess=delay_calc;
                else    
                     choice = questdlg('Set the delay:', ...
            'Input requested', ...
            'same as before','set visually','type','same as before');
        switch choice
            case 'set visually'
                    message = sprintf('Select the delay');
                    uiwait(msgbox(message,'Notification'));
                    delay_coord=ginput(1);                
             delay_guess=delay_coord(1);
            case 'type'
                    delay_guess = str2double(inputdlg({'Set the delay (s):'},'Setting delay'));
            case 'same as before'
                delay_guess=delay_calc;
        end
                end
                
                if fit_tau_DO==1
                    f_mod=bitset(f_mod,1);
                    tau_guess=tau_calc;
                else                    
                    choice = questdlg('Set tau:', ...
        'Setting tau', ...
        'same as before','manual input','same as before');
    switch choice
        case 'same as before'
            tau_guess=tau_calc;
         case 'manual input'
         tau_guess = str2double(inputdlg({'particle radius (um)'},'Setting tau'));
    end      
                end             
                        fitting_param=[p1_guess p2_guess delay_calc tau_calc f_mod];  
                        
                % setting the guess for fminsearch
                n=1;
                clear('guess');
                if bitget(f_mod,4)==1
                    guess(n)=p1_guess;
                    n=n+1;
                end
                if bitget(f_mod,3)==1
                    guess(n)=p2_guess;
                    n=n+1;
                end
                if bitget(f_mod,2)==1
                    guess(n)=delay_guess;
                    n=n+1;
                end
                if bitget(f_mod,1)==1
                    guess(n)=tau_guess;
                end
                                               
                fitted_parameters=fminsearch('Den_Otter_two_steps',guess,options_fit); % fitting
                        
                n=1; % retrieving the fitted parameters
                if bitget(f_mod,4)==0
                    p1=p1_guess;
                else
                    p1=fitted_parameters(n);
                    n=n+1;
                end
                if bitget(f_mod,3)==0
                    p2=p2_guess;
                else
                    p2=fitted_parameters(n);
                    n=n+1;
                end
                if bitget(f_mod,2)==0
                    delay_calc=delay_guess;
                else
                    delay_calc=fitted_parameters(n);
                    n=n+1;
                end
                if bitget(f_mod,1)==0
                    tau_calc=tau_guess;
                else
                    tau_calc=fitted_parameters(n);
                end   
                        
                        tau_calc_min=tau_calc/60;                        
                        r_s_free=A*B/(tau_calc*S*(A+B));%mol/m^2/s
                        ks_free=r_s_free/c_o;
                        
                        f_36_t0=norm_36_t0;
                        f_18g_t0=frac_18_O2_t0;
                        f_18b_t0=f18_sample_t0;
                        f_18_inf=frac_18_O2_inf;
                        beta_0=f_18_inf^2;
                        beta_1=f_18_inf*((f_18b_t0-f_18g_t0)*(p1+p2)+2*(f_18g_t0-f_18_inf));
                        beta_2=((f_18b_t0-f_18g_t0)*p1+f_18g_t0-f_18_inf)*((f_18b_t0-f_18g_t0)*p2+f_18g_t0-f_18_inf);
                        tau_2=tau_1*(1-beta_1/(2*f_18_inf*(f_18g_t0-f_18_inf)));
                        epsilon_1=2*f_18_inf*(f_18g_t0-f_18_inf);
                        epsilon_2=beta_2*tau_1/(tau_1-2*tau_2);
                        f_36_bis=f_18_inf^2+epsilon_1*exp(-(corr_time-delay_calc)/tau_1)+epsilon_2*exp(-2*(corr_time-delay_calc)/tau_1)+(f_36_t0-f_18_inf^2-epsilon_1-epsilon_2)*exp(-(corr_time-delay_calc)/tau_2);
                        f_18g_bis=f_18_inf+(f_18g_t0-f_18_inf)*exp(-(corr_time-delay_calc)/tau_1);
                        f_18b_bis=f_18_inf+(f_18b_t0-f_18_inf)*exp(-(corr_time-delay_calc)/tau_1);
                        f_34_bis=2*(f_18g-f_36_bis);
                        f_32_bis=1-f_34-f_36_bis;
                        
                        for i=1:size(corr_time)
                            if corr_time(i)<delay_calc
                                f_36_bis(i)=norm_36_t0;
                                f_34_bis(i)=norm_34(1);
                                f_32_bis(i)=norm_32(1);
                                f_18g_bis(i)=frac_18_O2_t0(1);
                                f_18b_bis(i)=f18_sample_t0(1);
                            end
                        end
                        
                               kdis_DO_free=ks_free*2/(p1+p2);
                        
                
                content_annotation={sprintf('\\Re_0 = %0.2e mol/m^2/s',r_s_free),sprintf('k_s = %0.3g m/s',ks_free),sprintf('\\tau = %0.4g min',tau_calc_min),sprintf('delay correction = %0.3g s',delay_calc)};
                plot_fit_frac_18_O2=figure('Name','Further fitting of the 2 steps model [beta]','NumberTitle','off');
                plot(corr_time_min,frac_18_O2,'k-','Linewidth',2)
                hold on;
                plot(corr_time_min(low_index:high_index),f_18g_bis(low_index:high_index),'r-','Linewidth',2)
                plot(corr_time_min,f_18b_bis,'k-.','Linewidth',1)
                plot(corr_time_min,f_18g_bis,'r--','Linewidth',2)
           

                xlabel('time (min)')
                ylabel('^1^8O fraction')
                axis([0 t_max_plot_scale 0 1])
                if background_correction==0
                title('Den Otter''s approach (surface exchange limitation)')
                elseif background_correction==1
                title('Background correction')   
                end
                legend('f_{18} (^*O_2)','best fit','f^{bulk (calc)}_{18}','Location','east')
                annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                        
                        
                        
                        
                 
                        
                        content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('Date : %s',experiment_date),sprintf('k_s = %0.3g m/s',ks_free),sprintf('k_{dis} = %0.3g m/s',kdis_DO_free),sprintf('\\tau = %0.4g min , delay = %0.3g s',tau_calc_min,delay_calc),sprintf('p_1 = %0.3g | p_2 = %0.3g',p1,p2)};
                        
                        plot_2steps_fitting2=figure('Name','Further fitting of the 2 steps model [beta]','NumberTitle','off');
                        hold on;
                        plot(corr_time_min,norm_32,'b-','Linewidth',2)
                        plot(corr_time_min,f_32_bis,'b--','Linewidth',1)
                        plot(corr_time_min,norm_34,'-','color',[1.0  0.0  1.0],'Linewidth',2)
                        plot(corr_time_min,f_34_bis,'--','color',[1.0  0.0  1.0],'Linewidth',1)
                        plot(corr_time_min,norm_36,'r-','Linewidth',2)
                        plot(corr_time_min,f_36_bis,'r--','Linewidth',1)
                        title('Den Otter''s approach (surface exchange limitation) [full freedom]')
                        xlabel('time (min)')
                        ylabel('Isotopologue fraction f_i')
                        axis([0 t_max_plot_scale 0 1])
                        legend('f_{32}','f^{calc}_{32}','f_{34}','f^{calc}_{34}','f_{36}','f^{calc}_{36}','Location','north','NumColumns',3)
                        annotation('textbox',[0.35 0.35 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                        
                        
                        % for export later:
                        ks_DO_free=ks_free;
                        %kdis_DO_free;
                        delay_DO_free=delay_calc;
                        tau_DO_free=tau_calc;
                        p1_DO_free=p1;
                        p2_DO_free=p2;
                        
                        
                end
                end
        end
        
        
        
        
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% fitting the dissociation rate (Boukamp 1994) %%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %S=m*SSA;
        choice = questdlg('Fit k_dis (Boukamp 1994)?', ...
            'k_dis fitting', ...
            'yes','no','yes');
        switch choice
            case 'no'
                
            case 'yes'
                y=log((frac_18_O2-frac_18_O2_inf)/(frac_18_O2_t0-frac_18_O2_inf));%% NB: natural logarithm ln
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting k_s as a funciton of time
                if k_as_a_function_of_time==1
                    for i=diff_step_Boukamp+1:size(corr_time)-diff_step_Boukamp
                        deriv_y(i,1)=(y(i+diff_step_Boukamp)-y(i-diff_step_Boukamp))/(corr_time(i+diff_step_Boukamp)-corr_time(i-diff_step_Boukamp));
                    end
                    deriv_y=deriv_y(diff_step_Boukamp+1:end-1,1);
                    r_s_t=abs(-deriv_y/S/(1/A+1/B));%mol/m^2/s
                    ks_t=r_s_t/c_o;
                    ks_t_cms=ks_t*100;
                    corr_time_deriv_min=corr_time_min(diff_step_Boukamp+1:end-diff_step_Boukamp-1,1);
                    plot_k_s_t=figure('Name','Evolution of ks as a function of time','NumberTitle','off');
                    hold on;
                    plot(corr_time_deriv_min,ks_t,'b.','Linewidth',2)
                    xlabel('time (min)')
                    ylabel('k_s (m/s)')
                    axis([0 t_max_plot_scale 0 max(ks_t_cms(1:150))])
                end
                
                % local determination of k_s (usual approach)
                plot_lnplot=figure('Name','Determination of the overall exchange rate','NumberTitle','off');
                hold on;
                plot(corr_time_min,y,'b.','Linewidth',2);
                xlabel('time (min)')
                ylabel('$\mathrm{ln\frac{f_{18}(t)-f_{18}(\infty)}{f_{18}(0)-f_{18}(\infty)}}$','Interpreter','Latex','FontSize',16)
                axis([0 t_max_plot_scale -inf inf])
                custom_message='Choose the fitting window:';
                [low_index high_index]=basic_restrain_time_range(corr_time_min);
                close(plot_lnplot);
                [LCpolynom,quality]=polyfit(corr_time(low_index:high_index,1),y(low_index:high_index,1),1);
                regLC=LCpolynom(2)+LCpolynom(1)*corr_time(low_index:high_index,1);
                regLCfull=LCpolynom(2)+LCpolynom(1)*corr_time;
                r_s=abs(-LCpolynom(1)/S/(1/A+1/B));%mol/m^2/s
                ks=r_s/c_o;
                ks_cms=ks*100;
                %plot fit
                content_annotation={sprintf('\\Re_0 = %0.2e mol/m^2/s',r_s),sprintf('k_s = %0.3g m/s',ks)};
                plot_lnplot_fit=figure('Name','Determination of the overall exchange rate - fit','NumberTitle','off');
                hold on;
                plot(corr_time_min,y,'b.',corr_time_min(low_index:high_index,1),regLC,'r-','Linewidth',2);
                xlabel('time (min)')
                ylabel('$\mathrm{ln\frac{f_{18}(t)-f_{18}(\infty)}{f_{18}(0)-f_{18}(\infty)}}$','Interpreter','Latex','FontSize',16)
                axis([0 t_max_plot_scale -inf inf])
                title('Boukamp''s approach  (surface exchange limitation)')
                annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                
                
                % Force 50/50mix fitting or not
                r_dis_guess=0.5*r_s;%mol/m^3/s
                if fiftyfiftymix_option==1
                    choice = questdlg('Set f_{18}(0)-f_{18}(inf) to 0?', ...
                        '5050mix GPA', ...
                        'yes','no','no');
                    switch choice
                        case 'no'
                            delta_f18=frac_18_O2_t0-frac_18_O2_inf;
                        case 'yes'
                            delta_f18=0;
                    end
                else
                    delta_f18=frac_18_O2_t0-frac_18_O2_inf;
                end
                
                r_dis_and_delay_guess=[r_dis_guess 120];
                r_dis_and_delay_fit = fminsearch('Boukamp_k',r_dis_and_delay_guess,options_fit);
                r_dis_fit=r_dis_and_delay_fit(1);
                delay_Boukamp=r_dis_and_delay_fit(2);
                p=r_dis_fit*S/A;
                q=r_s*S*(1/A+1/B);
                norm_36_calc=frac_18_O2_inf^2+(norm_36_t0-frac_18_O2_inf^2)*exp(-p*(corr_time-delay_Boukamp))-2*frac_18_O2_inf*(frac_18_O2_t0-frac_18_O2_inf)*(exp(-p*(corr_time-delay_Boukamp))-exp(-q*(corr_time-delay_Boukamp)))-(1-r_s/r_dis_fit*(1+A/B))^2/(1-2*r_s/r_dis_fit*(1+A/B))*(frac_18_O2_t0-frac_18_O2_inf)^2*(exp(-p*(corr_time-delay_Boukamp))-exp(-q*(corr_time-delay_Boukamp)));
                norm_34_calc=2*(frac_18_O2-norm_36_calc);
                norm_32_calc=1-norm_34_calc-norm_36_calc;
                kdis=r_dis_fit/c_o;
                content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('Date & time:%s - %s',experiment_date,experiment_time),sprintf('\\Re_s = %0.2e mol/m^2/s',r_s),sprintf('\\Re_{dis} = %0.2e mol/m^2/s',r_dis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis} = %0.3g m/s',kdis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg   |   \\rho = %0.3g g/cm^3',m_mg_f,rho_gcc),sprintf('SSA = %0.3g m^2/g',SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('M_{calc} = %0.4g g/mol   |   \\nu_{calc} = %0.3g',M_calc_i,nu_calc),sprintf('Delay correction = %0.3g s',delay_Boukamp)};
                plot_fit_dis=figure('Name','Boukamp''s approach','NumberTitle','off');
                hold on;
                plot(corr_time_min,norm_36,'b.',corr_time_min(low_index:high_index),norm_36_calc(low_index:high_index),'r-',corr_time_min,norm_36_calc,'r--','Linewidth',2);
                xlabel('time (min)')
                ylabel('^1^8O_2 fraction')
                title('Boukamp''s approach  (surface exchange limitation)')
                axis([0 t_max_plot_scale 0 1])
                annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                
                
                %choose a better window for the fit if not satisfactory
                new_range_choice_fit=1;
                if new_range_choice_fit==1
                    choice = questdlg('Choose the range for better fitting?', ...
                        'Fit range', ...
                        'yes','no','yes');
                    switch choice
                        case 'no'
                        case 'yes'
                            custom_message='Choose the fitting window:';
                            [low_index high_index]=basic_restrain_time_range(corr_time_min);
                            %                         r_dis_fit = fminsearch('Boukamp_k',r_dis_fit,options_fit);
                            %                         p=r_dis_fit*S/A;
                            %                         q=r_s*S*(1/A+1/B);
                            %                         norm_36_calc=frac_18_O2_inf^2+(norm_36_t0-frac_18_O2_inf^2)*exp(-p*corr_time)-2*frac_18_O2_inf*(frac_18_O2_t0-frac_18_O2_inf)*(exp(-p*corr_time)-exp(-q*corr_time))-(1-r_s/r_dis_fit*(1+A/B))^2/(1-2*r_s/r_dis_fit*(1+A/B))*(frac_18_O2_t0-frac_18_O2_inf)^2*(exp(-p*corr_time)-exp(-q*corr_time));
                            %                         kdis=r_dis_fit/c_o;
                            %                         content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('Date & time:%s - %s',experiment_date,experiment_time),sprintf('\\Re_s = %0.3g mol/m^2/s',r_s),sprintf('\\Re_{dis} = %0.3g mol/m^2/s',r_dis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis} = %0.3g m/s',kdis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg   |   \\rho = %0.3g g/cm^3',m_mg_f,rho_gcc),sprintf('SSA = %0.3g m^2/g',SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('M_{calc} = %0.4g g/mol   |   \\nu_{calc} = %0.3g',M_calc_i,nu_calc)};
                            %                         plot_fit_dis2=figure('Name','Title','NumberTitle','off');
                            %                         hold on;
                            %                         plot(corr_time_min,norm_36,'b.',corr_time_min(low_index:high_index),norm_36_calc(low_index:high_index),'r-',corr_time_min,norm_36_calc,'r--','Linewidth',2);
                            %                         xlabel('time (min)')
                            %                         ylabel('^1^8O_2 fraction')
                            %                         axis([0 t_max_plot_scale 0 1])
                            %                         annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                            r_dis_and_delay_guess=[r_dis_guess 120];
                            r_dis_and_delay_fit = fminsearch('Boukamp_k',r_dis_and_delay_guess,options_fit);
                            r_dis_fit=r_dis_and_delay_fit(1);
                            delay_Boukamp=r_dis_and_delay_fit(2);
                            p=r_dis_fit*S/A;
                            q=r_s*S*(1/A+1/B);
                            norm_36_calc=frac_18_O2_inf^2+(norm_36_t0-frac_18_O2_inf^2)*exp(-p*(corr_time-delay_Boukamp))-2*frac_18_O2_inf*(frac_18_O2_t0-frac_18_O2_inf)*(exp(-p*(corr_time-delay_Boukamp))-exp(-q*(corr_time-delay_Boukamp)))-(1-r_s/r_dis_fit*(1+A/B))^2/(1-2*r_s/r_dis_fit*(1+A/B))*(frac_18_O2_t0-frac_18_O2_inf)^2*(exp(-p*(corr_time-delay_Boukamp))-exp(-q*(corr_time-delay_Boukamp)));
                            norm_34_calc=2*(frac_18_O2-norm_36_calc);
                            norm_32_calc=1-norm_34_calc-norm_36_calc;
                            kdis=r_dis_fit/c_o;
                            content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('Date & time:%s - %s',experiment_date,experiment_time),sprintf('\\Re_s = %0.2e mol/m^2/s',r_s),sprintf('\\Re_{dis} = %0.2e mol/m^2/s',r_dis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis} = %0.3g m/s',kdis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg   |   \\rho = %0.3g g/cm^3',m_mg_f,rho_gcc),sprintf('%s = %0.3g m^2/g',SSA_display,SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('M_{calc} = %0.4g g/mol   |   \\nu_{calc} = %0.3g',M_calc_i,nu_calc),sprintf('Delay correction = %0.3g s',delay_Boukamp)};
                            plot_fit_dis=figure('Name','Boukamp''s approach','NumberTitle','off');
                            hold on;
                            plot(corr_time_min,norm_36,'b.',corr_time_min(low_index:high_index),norm_36_calc(low_index:high_index),'r-',corr_time_min,norm_36_calc,'r--','Linewidth',2);
                            xlabel('time (min)')
                            ylabel('^1^8O_2 fraction')
                            title('Boukamp''s approach  (surface exchange limitation)')
                            axis([0 t_max_plot_scale 0 1])
                            annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                            
                            
                            
                            
                            
                    end
                end
                
                %Try exponential decay fit (usefull for quartz background /any catalyst dissociation kinetics)
                if quartz_background_fit==1
                    choice = questdlg('Try exponential decay fit (usefull for quartz background)?', ...
                        'Fit range', ...
                        'yes','no','no');
                    switch choice
                        case 'no'
                        case 'yes'
                            norm_norm_36=(norm_36-norm_36_t0)/(norm_36_inf-norm_36_t0);
            f_mod=3;
            fitting_param=[0 0 f_mod];
            
                            extract_tau_data={corr_time,norm_norm_36,low_index,high_index};
                            guess=[A/r_dis_fit/S 120];
                            tau_and_delay_fit = fminsearch('exponential_decay_extract_tau_and_delay',guess,options_fit);
                            tau=tau_and_delay_fit(1);
                            delay=tau_and_delay_fit(2);
                            norm_36_calc=norm_36_t0+(norm_36_inf-norm_36_t0)*(1-exp(-(corr_time-delay)./tau));
                            r_dis_surface_catalysis_fit=A/tau/S;
                            k_dis_surface_catalysis=r_dis_surface_catalysis_fit/c_o;
                            content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('\\Re_s = %0.2e mol/m^2/s',r_s),sprintf('\\Re_{dis}_(_c_a_t_) = %0.2e mol/m^2/s',r_dis_surface_catalysis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis}_(_c_a_t_) = %0.3g m/s',k_dis_surface_catalysis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg',m_mg_f),sprintf('\\rho = %0.3g g/cm^3',rho_gcc),sprintf('%s = %0.3g m^2/g',SSA_display,SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('\\nu_{calc} = %0.3g',nu_calc)};
                            plot_fit_dis3=figure('Name','Boukamp''s approach','NumberTitle','off');
                            hold on;
                            plot(corr_time_min,norm_36,'b.',corr_time_min,norm_36_calc,'r-','Linewidth',2);
                            xlabel('time (min)')
                            ylabel('^1^8O_2 fraction')
                            title('Boukamp''s approach - modified for quartz background fit')
                            axis([0 t_max_plot_scale 0 1])
                            annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                    end
                end
                
                % %                 %Choose kdis when the fit is still too far (hidden because ~deprecated)
                % %                 choice = questdlg('Choose \\Re_{dis} initial value for better fitting?', ...
                % %                     'Fit range', ...
                % %                     'yes','no','no');
                % %                 switch choice
                % %                     case 'no'
                % %                     case 'yes'
                % %                         r_dis_guess=input('\\Re_{dis} guess(mol/m^2/s):');
                % %                         r_dis_fit = fminsearch('Boukamp_k',r_dis_guess,options_fit);
                % %                         p=r_dis_fit*S/A;
                % %                         q=r_s*S*(1/A+1/B);
                % %                         norm_36_calc=frac_18_O2_inf^2+(norm_36_t0-frac_18_O2_inf^2)*exp(-p*corr_time)-2*frac_18_O2_inf*(frac_18_O2_t0-frac_18_O2_inf)*(exp(-p*corr_time)-exp(-q*corr_time))-(1-r_s/r_dis_fit*(1+A/B))^2/(1-2*r_s/r_dis_fit*(1+A/B))*(frac_18_O2_t0-frac_18_O2_inf)^2*(exp(-p*corr_time)-exp(-q*corr_time));
                % %                         kdis=r_dis_fit/c_o;
                % %                         content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('Date & time:%s - %s',experiment_date,experiment_time),sprintf('\\Re_s = %0.3g mol/m^2/s',r_s),sprintf('\\Re_{dis} = %0.3g mol/m^2/s',r_dis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis} = %0.3g m/s',kdis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg   |   \\rho = %0.3g g/cm^3',m_mg,rho_gcc),sprintf('SSA = %0.3g m^2/g',SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('M_{calc} = %0.4g g/mol   |   \\nu_{calc} = %0.3g',M_calc_i,nu_calc)};
                % %                         % content_annotation={sprintf('%s (fit)',raw_data_filename),sprintf('\\Re_s = %0.3g mol/m^2/s',r_s),sprintf('\\Re_{dis} = %0.3g mol/m^2/s',r_dis_fit),sprintf('k_s = %0.3g m/s',ks),sprintf('k_{dis} = %0.3g m/s',kdis),sprintf('P = %0.3g mbar   |   V = %0.3g mL',P_mbar,V_ml),sprintf('T_{i} = %0.1f %cC',T_C,char(176)),sprintf('m = %0.3g mg',m_mg),sprintf('\\rho = %0.3g g/cm^3',rho_gcc),sprintf('SSA = %0.3g m^2/g',SSA),sprintf('M_{input} = %0.4g g/mol   |   \\nu_{input} = %0.3g',M,nu),sprintf('\\nu_{calc} = %0.3g',nu_calc)};
                % %                         plot_fit_dis4=figure('Name','Title','NumberTitle','off');
                % %                         hold on;
                % %                         plot(corr_time_min,norm_36,'b.',corr_time_min,norm_36_calc,'r-','Linewidth',2);
                % %                         xlabel('time (min)')
                % %                         ylabel('^1^8O_2 fraction')
                % %                         axis([0 t_max_plot_scale 0 1])
                % %                         annotation('textbox',[0.35 0.6 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                % %                 end
                
        end %end of Boukamp approach
        
end

%save the figures and values
choice = questdlg('Save the figures, fitted values and data?', ...
    'Saving', ...
    'yes','no','yes');
switch choice
    case 'no'
    case 'yes'
        mkdir(fullfile(saving_path,sprintf('%s',timestamp)));%create a directory for each fit
        saving_directory=strcat(saving_path,sprintf('%s',timestamp),'\');%set the path for the saving directory
        
        %% %%%%%%%%%%%%%%%%%%%%% %%
        %% exporting the vectors %%
        %% %%%%%%%%%%%%%%%%%%%%% %%
        
        %% corrected & normalised raw data +18O fractions:
        if save_norm_data==1            
                potential_data_to_save={'t [min]' 'corr_time_min';
                    'norm_16' 'norm_16';
                    'norm_17' 'norm_17';
                    'norm_18' 'norm_18';
                    'norm_19' 'norm_19';
                    'norm_20' 'norm_20';
                    'norm_21' 'norm_21';
                    'norm_22' 'norm_22';
                    'norm_23' 'norm_23';
                    'norm_28' 'norm_28';
                    'norm_30' 'norm_30';
                    'norm_32' 'norm_32';
                    'norm_34' 'norm_34';
                    'norm_36' 'norm_36';
                    'norm_40' 'norm_40';
                    'norm_44' 'norm_44';
                    'norm_46' 'norm_46';
                    'norm_48' 'norm_48';
                    'frac_18_O2' 'frac_18_O2';
                    'frac_18_O' 'frac_18_O'};         
         export_fit_path=fullfile(saving_directory,'norm_data_and_O18_fraction.csv');
         analysis_log(potential_data_to_save,export_fit_path);
        end
        
        %% Den Otter fit:
        if exist('f_32')==1
potential_data_to_save={'t [min]' 'corr_time_min';
                'norm_32' 'norm_32';
                'norm_34' 'norm_34';
                'norm_36' 'norm_36';
                'f18_O2' 'frac_18_O2';
                'f18_O2_calc' 'frac_18_O2_calc_DenOtter';
                '32_calc' 'f_32';
                '34_calc' 'f_34';
                '36_calc' 'f_36';
                'p1p2_32_calc' 'f_32_bis';
                'p1p2_34_calc' 'f_34_bis';
                'p1p2_36_calc' 'f_36_bis' };            
         export_fit_path=fullfile(saving_directory,'Den_Otter_fit.csv');
         analysis_log(potential_data_to_save,export_fit_path);
        end
        
        %% Mizusaki fit:
        if exist('frac_18_O2_calc_Mizu_D')==1 | exist('frac_18_O2_calc_Mizu_kD')==1
            
            
            potential_data_to_save={'t [min]' 'corr_time_min';
                'f18_O2' 'frac_18_O2';
                'f18_O2_calc_Mizu_D' 'frac_18_O2_calc_Mizu_D';
                'frac_18_O2_calc_Mizu_kD' 'frac_18_O2_calc_Mizu_kD'};
         export_fit_path=fullfile(saving_directory,'Mizusaki_fit.csv');
         analysis_log(potential_data_to_save,export_fit_path);
        end
        
        
        %% Boukamp fit:
        if exist('y')==1 | exist('regLC')==1 | exist('norm_36_calc')==1
            potential_data_to_save={'t [min]' 'corr_time_min';
                'norm_36' 'norm_36';
                'frac_18_O2' 'frac_18_O2';
                'norm_36_calc' 'norm_36_calc';
                'log_F' 'y';
                'log_F_fit' 'regLCfull'};
         export_fit_path=fullfile(saving_directory,'Boukamp_fit.csv');
         analysis_log(potential_data_to_save,export_fit_path);
        end
        
        %% saving all in one:
        if save_all_in_one==1
            potential_data_to_save={'t [min]' 'corr_time_min';
                    'norm_16' 'norm_16';
                    'norm_17' 'norm_17';
                    'norm_18' 'norm_18';
                    'norm_19' 'norm_19';
                    'norm_20' 'norm_20';
                    'norm_21' 'norm_21';
                    'norm_22' 'norm_22';
                    'norm_23' 'norm_23';
                    'norm_28' 'norm_28';
                    'norm_30' 'norm_30';
                    'norm_32' 'norm_32';
                    'norm_34' 'norm_34';
                    'norm_36' 'norm_36';
                    'norm_40' 'norm_40';
                    'norm_44' 'norm_44';
                    'norm_46' 'norm_46';
                    'norm_48' 'norm_48';
                    'frac_18_O2' 'frac_18_O2';
                    'frac_18_O' 'frac_18_O';
                    'frac_18_O2_calc_Mizu_D' 'frac_18_O2_calc_Mizu_D';
                    'frac_18_O2_calc_Mizu_kD' 'frac_18_O2_calc_Mizu_kD';
                    'frac_18_O2_calc_DenOtter' 'frac_18_O2_calc_DenOtter';
                    'DenOtter_p_32_calc' 'f_32';
                    'DenOtter_p_34_calc' 'f_34';
                    'DenOtter_p_36_calc' 'f_36';
                    'DenOtter_p1p2_32_calc' 'f_32_bis';
                    'DenOtter_p1p2_34_calc' 'f_34_bis';
                    'DenOtter_p1p2_36_calc'  'f_36_bis';
                    'Boukamp_norm_32_calc' 'norm_32_calc';
                    'Boukamp_norm_34_calc' 'norm_34_calc';
                    'Boukamp_norm_36_calc' 'norm_36_calc';
                    'Boukamp_log_F' 'y';
                    'Boukamp_log_F_fit' 'regLCfull'};  
         export_fit_path=fullfile(saving_directory,'data_and_fits.csv');
         analysis_log(potential_data_to_save,export_fit_path);
        end
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        %% exporting fitted parameters %%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        option_export_fit_parameters=1;
        if option_export_fit_parameters==1
            %saves the experimental values and fitted parameters in the log
            %file in the local log folder; these can be saved to the main
            %database later
            
            
potential_parameters_to_save={'timestamp' 'timestamp';
                                          'sample' 'raw_data_filename';
                                          'm_f [g]' 'm_f';
                                          'm_i [g]' 'm_i';
                                          'm_f (calc) [g]' 'm_f_calc';
                                          'm_i (calc) [g]' 'm_i_calc';
                                          'M [g/mol]' 'M';
                                          'nu' 'nu';
                                          'rho [g/m³]' 'rho';
                                          'SSA [m²/g]' 'SSA_log';
                                          'SSA (calc) [m²/g]' 'SSA_calc';
                                          'V [m³]' 'V';
                                          'P [Pa]' 'P';
                                          'T [K]' 'T';
                                          'nOgas [mol]' 'A';
                                          'nOsolid [mol]' 'B' ;
                                          'cO [mol/m³]' 'c_o';
                                          'particle_radius [m]' 'radius_log';
                                          'particle_radius (calc) [m]' 'radius_calc';
                                          'sphericity'  'sphericity';
                                          'Mizu_D(D only) [m²/s]' 'D_Mizu_D';
                                          'delay_Mizu [s]' 'delay_Mizu_D';
                                          'Mizu_k [m/s]'  'k_Mizu_kD';
                                          'Mizu_D [m²/s]'  'D_Mizu_kD';
                                          'Mizu_Lc [um]' 'L_C_Mizu_um';
                                          'Mizu_delay [s]'  'delay_Mizu_kD';
                                          'Mizu_radius [um]' 'radius_Mizu_kD';
                                          'ks_DO [m/s]' 'ks_DO';
                                          'kdis_DO [m/s]' 'kdis_DO';
                                          'delay_DO [s]' 'delay_DO';
                                          'p_DO' 'p_DO';
                                          'ks_DO_free [m/s]' 'ks_DO_free';
                                          'kdis_DO_free [m/s]' 'kdis_DO_free';
                                          'delay_DO_free [s]' 'delay_DO_free';
                                          'p1_DO_free' 'p1_DO_free';
                                          'p2_DO_free' 'p2_DO_free';
                                          'ks_B [m/s]'  'ks';
                                          'rs_B [mol/m²/s]' 'r_s';
                                          'kdis_B [m/s]' 'kdis';
                                          'rdis_B [mol/m²/s]'   'r_dis_fit';
                                          'nu_calc' 'nu_calc';
                                          'f18_t0'  'frac_18_O2_t0';
                                          'f18_inf'  'frac_18_O2_inf';
                                          'f36_t0'  'norm_36_t0';
                                          'f36_inf'  'norm_36_inf';
                                          'f18_t0(sample)' 'f18_sample_t0'}; 
savloarch_log(potential_parameters_to_save,fullfile(saving_path,'fitting_log.csv'),'archive');
            
 
        
        %% %%%%%%%%%%%%%%%%%%%%% %%
        %% exporting the figures %%
        %% %%%%%%%%%%%%%%%%%%%%% %%
        option_export_figures=1;
        if option_export_figures==1
            fig_res=sprintf('-r%d',fig_res);
            message=msgbox('Saving figures, please wait...','Notification'); 
            
            list_of_figures_to_save={'raw_data' 'plot_raw323436';
            'f18' 'plot_frac_18_O2';
            'Boukamp_1' 'plot_lnplot_fit';
            'Boukamp_2' 'plot_fit_dis';            
            'DO_tau' 'plot_charac_time';
           'DO_f18fit' 'plot_fit_frac_18_O2';                     
           'Den_Otter' 'plot_2steps_fitting';           
           'Den_Otter_2' 'plot_2steps_fitting2';           
           'thvsxp' 'plot_thvsexp';           
           'deltavsdiff' 'plot_deltavsdiff';           
           'Boukamp_2_2' 'plot_fit_dis2';           
           'Boukamp_2_3' 'plot_fit_dis3';           
           'Boukamp_2_4' 'plot_fit_dis4';   
           'water_and_co' 'plot_frac_17_19_20_21_22' ;          
           '14-40' 'plot_frac_14_40';           
           'carbonates_and_co' 'plot_frac_28_30_44_46_48';
           'Boukamp_k_s_t' 'plot_k_s_t';           
           'Mizusaki_D_an_k' 'plot_fit_Mizusaki_D_an_k';           
           'Mizusaki_D_limited' 'plot_fit_Mizusaki_D_limited'};            
            
            %};%name of files
            %          eval...   if exist('plot_raw323436')==1
            %                 % saveas(plot_raw323436,strcat(saving_directory,sprintf('%s-raw.jpg',raw_data_filename)),'jpg');
            %                 print('-dpng',fig_res,plot_raw323436,strcat(saving_directory,sprintf('%s-raw.png',raw_data_filename)));
            %                 saveas(plot_raw323436,strcat(saving_directory,sprintf('%s-raw.fig',raw_data_filename)),'fig');
            %             end
            
            if exist('plot_raw323436')==1
                % saveas(plot_raw323436,strcat(saving_directory,sprintf('%s-raw.jpg',raw_data_filename)),'jpg');
                print('-dpng',fig_res,plot_raw323436,strcat(saving_directory,sprintf('%s-raw.png',raw_data_filename)));
                saveas(plot_raw323436,strcat(saving_directory,sprintf('%s-raw.fig',raw_data_filename)),'fig');
            end
            
%             if exist('plot_norm323436shift')==1
%                 print('-dpng',fig_res,plot_norm323436shift,strcat(saving_directory,sprintf('%s-norm.png',raw_data_filename)));
%                 saveas(plot_norm323436shift,strcat(saving_directory,sprintf('%s-norm.fig',raw_data_filename)),'fig');
%             end
            
            if exist('plot_frac_18_O2')==1
                print('-dpng',fig_res,plot_frac_18_O2,strcat(saving_directory,sprintf('%s-f18.png',raw_data_filename)));
                saveas(plot_frac_18_O2,strcat(saving_directory,sprintf('%s-f18.fig',raw_data_filename)),'fig');
            end
            
            %         if exist('plot_lnplot')==1
            %             print('-dpng',fig_res,plot_lnplot,strcat(saving_directory,sprintf('%s-Boukamp_1.png',raw_data_filename)));
            %             saveas(plot_lnplot,strcat(saving_directory,sprintf('%s-Boukamp_1.fig',raw_data_filename)),'fig');
            %         end
            
            if exist('plot_lnplot_fit')==1
                print('-dpng',fig_res,plot_lnplot_fit,strcat(saving_directory,sprintf('%s-Boukamp_1.png',raw_data_filename)));
                saveas(plot_lnplot_fit,strcat(saving_directory,sprintf('%s-Boukamp_l.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_fit_dis')==1
                print('-dpng',fig_res,plot_fit_dis,strcat(saving_directory,sprintf('%s-Boukamp_2.png',raw_data_filename)));
                saveas(plot_fit_dis,strcat(saving_directory,sprintf('%s-Boukamp_2.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_charac_time')==1
                print('-dpng',fig_res,plot_charac_time,strcat(saving_directory,sprintf('%s-DO_tau.png',raw_data_filename)));
                saveas(plot_charac_time,strcat(saving_directory,sprintf('%s-DO_tau.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_fit_frac_18_O2')==1
                print('-dpng',fig_res,plot_fit_frac_18_O2,strcat(saving_directory,sprintf('%s-DO_f18fit.png',raw_data_filename)));
                saveas(plot_fit_frac_18_O2,strcat(saving_directory,sprintf('%s-DO_f18fit.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_2steps_fitting')==1
                print('-dpng',fig_res,plot_2steps_fitting,strcat(saving_directory,sprintf('%s-Den_Otter.png',raw_data_filename)));
                saveas(plot_2steps_fitting,strcat(saving_directory,sprintf('%s-Den_Otter.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_2steps_fitting2')==1
                print('-dpng',fig_res,plot_2steps_fitting2,strcat(saving_directory,sprintf('%s-Den_Otter_2.png',raw_data_filename)));
                saveas(plot_2steps_fitting2,strcat(saving_directory,sprintf('%s-Den_Otter_2.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_thvsexp')==1
                print('-dpng',fig_res,plot_thvsexp,strcat(saving_directory,sprintf('%s-thvsxp.png',raw_data_filename)));
                saveas(plot_thvsexp,strcat(saving_directory,sprintf('%s-thvsxp.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_deltavsdiff')==1
                print('-dpng',fig_res,plot_deltavsdiff,strcat(saving_directory,sprintf('%s-deltavsdiff.png',raw_data_filename)));
                saveas(plot_deltavsdiff,strcat(saving_directory,sprintf('%s-deltavsdiff.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_fit_dis2')==1
                print('-dpng',fig_res,plot_fit_dis2,strcat(saving_directory,sprintf('%s-Boukamp_2_2.png',raw_data_filename)));
                saveas(plot_fit_dis2,strcat(saving_directory,sprintf('%s-Boukamp_2_2.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_fit_dis3')==1
                print('-dpng',fig_res,plot_fit_dis3,strcat(saving_directory,sprintf('%s-Boukamp_2_3.png',raw_data_filename)));
                saveas(plot_fit_dis3,strcat(saving_directory,sprintf('%s-Boukamp_2_3.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_fit_dis4')==1
                print('-dpng',fig_res,plot_fit_dis4,strcat(saving_directory,sprintf('%s-Boukamp_2_4.png',raw_data_filename)));
                saveas(plot_fit_dis4,strcat(saving_directory,sprintf('%s-Boukamp_2_4.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_frac_17_19_20_21_22')==1
                print('-dpng',fig_res,plot_frac_17_19_20_21_22,strcat(saving_directory,sprintf('%s-water_and_co.png',raw_data_filename)));
                saveas(plot_frac_17_19_20_21_22,strcat(saving_directory,sprintf('%s-water_and_co.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_frac_14_40')==1
                print('-dpng',fig_res,plot_frac_14_40,strcat(saving_directory,sprintf('%s-14-40.png',raw_data_filename)));
                saveas(plot_frac_14_40,strcat(saving_directory,sprintf('%s-14-40.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_frac_28_30_44_46_48')==1
                print('-dpng',fig_res,plot_frac_28_30_44_46_48,strcat(saving_directory,sprintf('%s-carbonates_and_co.png',raw_data_filename)));
                saveas(plot_frac_28_30_44_46_48,strcat(saving_directory,sprintf('%s-carbonates_and_co.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_k_s_t')==1
                print('-dpng',fig_res,plot_k_s_t,strcat(saving_directory,sprintf('%s-Boukamp_k_s_t.png',raw_data_filename)));
                saveas(plot_k_s_t,strcat(saving_directory,sprintf('%s-Boukamp_k_s_t.fig',raw_data_filename)),'fig');
            end
            
            if exist('plot_fit_Mizusaki_D_an_k')==1
                print('-dpng',fig_res,plot_fit_Mizusaki_D_an_k,strcat(saving_directory,sprintf('%s-Mizusaki_D_an_k.png',raw_data_filename)));
                saveas(plot_fit_Mizusaki_D_an_k,strcat(saving_directory,sprintf('%s-Mizusaki_D_an_k.fig',raw_data_filename)),'fig');
            end
            if exist('plot_fit_Mizusaki_D_limited')==1
                print('-dpng',fig_res,plot_fit_Mizusaki_D_limited,strcat(saving_directory,sprintf('%s-Mizusaki_D_limited.png',raw_data_filename)));
                saveas(plot_fit_Mizusaki_D_limited,strcat(saving_directory,sprintf('%s-Mizusaki_D_limited.fig',raw_data_filename)),'fig');
            end
            
            delete(message);
        end

%save the log to the main db (by defaut in the same folder than the current script)
if option_log_to_database==1
    choice = questdlg('log to the main database?', ...
        'Saving', ...
        'yes','no','yes');
    switch choice
        case 'no'
        case 'yes'
            if exist(GPA_log_db_path,'file') == 2
                GPA_log_db=readtable(GPA_log_db_path);
                GPA_log_db=[GPA_log_db;new_log_line];%append the line
            else
                GPA_log_db = [new_log_line];%write the 1st line in the log
            end
            writetable(GPA_log_db,GPA_log_db_path,'Delimiter',',');%write the fitting log file
    end
end

end
% saving the experimental parameters that may have change
%save_conditions();
savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save
end

if option_update_weigth_and_nu==1
    choice = questdlg('update the experimental parameters with the new molar weight nu and particle size?', ...
        'Update', ...
        'yes','no','no');
    switch choice
        case 'no'
        case 'yes'
            nu=nu_calc;
            M=M_calc_i;
            radius=radius_calc;
            %save_conditions();
            savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save

    end
end

if xor(calculate_SSA_from_radius==1,calculate_radius_from_SSA==1)==1
    choice = questdlg('Update the experimental parameters with the current SSA and particle radius?', ...
        'Update', ...
        'yes','no','no');
    switch choice
        case 'no'
        case 'yes'
            SSA_log=SSA;
            radius_log=radius;
            %save_conditions();
            savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save

    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%
%% INTERNAL ROUTINES %%
%%%%%%%%%%%%%%%%%%%%%%%

function [low_index  high_index] = restrain_time_range(restricable_time)
% Returns the lower and higher indexes of a time array by selecting
% visually its lower and higher time values.
% A custom message (must be assigned to a global variable "custom_message"
% can be displayed instead of the default one.
% Custom indexes (must be assigned to a global variable "custom_indexes"
% can be chosen instead of the default ones for the standard choice.
global custom_message custom_indexes
if isempty(custom_message)
    message= 'Restraint of the range:';
else
    message=custom_message;
end
if isempty(custom_indexes)
    low=10;
    high=size(restricable_time,1);
else
    low=custom_indexes(1);
    high=custom_indexes(2);
end
question = questdlg(message, ...
    'Range restraint', ...
    'std','none','by hand','std');
switch question
    case 'std'
        low_index=low;% can be good in some cases to change it...
        high_index=high;
    case 'none'
        low_index=1;
        high_index=size(restricable_time,1);
    case 'by hand'
        message = sprintf('Select the fitting range (2 points along x axis).');
        uiwait(msgbox(message));
        fitting_range=ginput(2);
        mini=min(fitting_range(:,1));
        maxi=max(fitting_range(:,1));
        %determines the upper and lower indexes
        s=size(restricable_time,1);
        low_index=1;
        for a =1:s
            if restricable_time(a)>mini
                break
            else
                low_index=a;
            end
        end
        high_index=s;
        for a = s:-1:1
            if restricable_time(a)<maxi
                break
            else
                high_index=a;
            end
        end
end
%reset the customised parameters:
custom_message='';
custom_indexes=[];
end





function [low_index  high_index] = basic_restrain_time_range(restricable_time)
% Returns the lower and higher indexes of a time array by selecting
% visually its lower and higher time values.
global custom_message
if isempty(custom_message)
    message = sprintf('Select the fitting range (2 points along x axis).');
else
    message=custom_message;
end
uiwait(msgbox(message));
fitting_range=ginput(2);
mini=min(fitting_range(:,1));
maxi=max(fitting_range(:,1));
%determines the upper and lower indexes
s=size(restricable_time,1);
low_index=1;
for a =1:s
    if restricable_time(a)>mini
        break
    else
        low_index=a;
    end
end
high_index=s;
for a = s:-1:1
    if restricable_time(a)<maxi
        break
    else
        high_index=a;
    end
end
end


