% ExampleExecutive.m
% This example shows a model setup for simulation of
% a series of photochemical chamber experiments.
% Read comments in each section for a guided tour.
%
% This has been updated to run the Washington Aerosol Module (WAM)
% 1/23/17 ELD
% 7/20/19 JAT

%% Notes


%% 
clear

%% Chemistry Inputs (AJC)


test_dat = [1   500    0.03   3    15     100    0    0    27.0    0.5     9.3 	0.4     3.6 	0.1     25.4];
 
test_vars = {'Run','LIMONENE_ppbv','HONO_ppmv','NO_ppbv','NOx_ppbv','AS_vol_um3_cm3','err_AS_vol','O3_max_ppbv','dOM_max','err_dOM','dOM_final','err_dOM_final','Y','err_Y','T_C'};



%% Set up run conditions - use solar zenith angle (SZA) as switch  
expts = [1]; %[2 9 10 11 12 13 14]; %list experiment runs from above table to model

for i = 1:length(expts) 

%% OPTIONS
    %{
"Verbose" can be set from 0-3; this just affects the level of detail printed to the command
  window regarding model progress.
"EndPointsOnly" is set to 0 because we want output to include all concentrations along each model step.
"LinkSteps" is set to 0 because each step is fully independent.
"Repeat" is set to 1 because we only want to go through each step once.
"IntTime" is the integration time for each step. Let's pretend each experiment lasts 3 hours.
"SavePath" is commented out, so output will be saved in a dated folder in the \Runs\ folder.
"Do_GP" calls WAM to calculate gas-particle partitioning
"GPupdateTime" is number of model seconds to calculate gas-particle
partitioning
"Do_Wall_Partitioning" uses Krechmer et al parameterization of c* driven
wall loss. use with care. If enabled set k_wall_lossV = 0.
    %}
    
    ModelOptions.Verbose       = 1;
    ModelOptions.EndPointsOnly = 0;
    ModelOptions.LinkSteps     = 0;
    ModelOptions.Repeat        = 1;
    ModelOptions.IntTime       = 10622.00; %seconds
    ModelOptions.SavePath      = '';
    ModelOptions.Do_GP         = 1 ; % run GP partitioning code, 0 calls regular F0AM_ModelCore
    ModelOptions.GPupdateTime = 1; %in seconds. add to Model Options!!
    ModelOptions.Do_Wall_Partitioning = 0;   
    
    expt = expts(i);
    numruns = length(expt); %number of individual chamber experiments
    
    SZA = 50.*ones(1,numruns);%note this SZA matches reported j_H2O2 within a few percent, but causes j_NO2 to be
    % 43% higher than reported from measured photon flux. JNO2 will be
    % overridden with measured
    
    %set "input concentrations" - i.e. those before lights on in chamber.
    %iC5H8 = test_dat(expt,2);%ppb
    iLIMONENE = test_dat(expt,2); 
    
    iHONO = test_dat(expt,3);
    iHONO = iHONO.*1000; %convert to ppb
%     
    iNO = test_dat(expt,4);%ppb
    iNO2 = (test_dat(expt,5)-iNO); %ppb
    
    KT = test_dat(expt,end)+273.13;
    KP = 1013.*ones(1,numruns); %ASSUME 1ATM pressure
    
    %% WAM particle and wall initialization
    
    %can make these inputs
    Nseed = 1e3; %number seed per cm3
    Rp = 25e-7; %seed particle radius (cm) (50nm diameter seed particles used in expt)
    SAi = Nseed.*4*pi*(Rp)^2; %initial seed surface area cm2/cm3
    Vseed = SAi.*Rp./3; %initial seed volume cm3/cm3
    
    MWoa_init = 175; %guess at average OA molecular mass (gets updated in model)
    Coa_ugm3 = .002; %initial OA (ug/m3). going down to .001 throws errors, this is the lowest we can go
    iCoa = ugm3_to_mr(Coa_ugm3,MWoa_init).*1e9; %ppb. assume average MW of OA constituents 175g/mol
    
    %Note, if you change gamma (e.g. mass accommodation coefficient) you should also change ModelOptions.GPupdateTime to 10
    %for gama = 1, recommed use ModelOptions.GPupdateTime = 10, otherwise computation errors result from stiff equations.
    gama = 1; 
    
    wall_params.wcstar_thresh = 1E6; %cstar threshold to invoke wall partitioning
    
    %equivalent absorbing organic concentration of the wall material C_wall
    %in ug/m3. 
    
    %Use C_wall = [] if you want to use the Krechmer et al parameterization based on c*'s. See F0AM_WallPartitioning_Generator.m
    %Otherwise the C_wall value you enter here will be used for all
    %species.
    
    wall_params.C_wall = []; %30000; %equivalent wall absorbing mass concentration ug/m3
    
    wall_params.kwall_transport = 1e-5; %per second timescale to mix vapors to/from wall surface 
    
    %% Meteorology
    Met = {...
        %   names       values
        'P'         KP; %1013.*ones(1,numruns)       ;...    %Pressure, mbar
        'T'         KT; %kroll_2006_dat(1:numruns,end)+273.13        ;...    %Temperature, K
        'RH'        1;         ;...    %Relative Humidity, percent
        'SZA'       SZA;...  % 50 Gives roughly right jvalue for H2O2 in Caltech 2006 chamber, need to update jno2
        %     'LFlux'     'ExampleLightFlux.txt'     ; %Text file for radiation spectrum
        'jcorr'     3.3e-4 ; % for kroll chamber %ones(1,numruns);...           %light attenuation factor
        'kwall_loss' 0.*ones(1,numruns);...% 6E-5  ;...                         %wall loss rate of particles, per second. From Crump&Seinfeld 1983 AS&T (Table 2) for ~100nm particles=1E-5, John suggests 3E-5 (0.1/hr)
        'kwall_lossV' 0;%0.*ones(1,numruns);... % sometimes 1e-5 used                         %wall loss rate of condensable vapors, per second. 2E-4 wall loss reasonable (E-4 fastest rate)
        'Cstar_threshold' 100.*ones(1,numruns);...   %Cstar threshold to filter compounds that can partition. ug/m3
        'kdil'      0.*ones(1,numruns)          ;...%Dilution rate constant, per second
        'EF_C5H8'   0  ;...                              %emission factor to supply isoprene, continuous flow
        'EF_H2O2'   0  ;...                              %emission factor to supply H2O2, continuous flow
        'EF_NO'     0  ;...                                  %emission factor to supply NOx, continuous flow
        'J3'        0.000175./600;...                           %Kroll state 0.00029/60,have used 0.00025/60 NOTE: 3.7e-6 s^-1 for Krechmer
        'J4'        0.29./60;... %was 0.29/60                           %jvalue in Kroll seems too high (based on O3 production) Kroll 2006 %Override calcuated J-value for NO2 + hv -->
        'J5'        2E-4;...                                 %Override calculated J-value for NO3+hv--> NO   original=.02485
        'J6'        4.6E-4;...                               %Override calculated J-value for NO3+hv--> NO2+O3   original=0.1747
        'LFlux'     '340nm_update.txt'     ; %Text file for radiation spectrum
        };
    
    %% CHEMICAL CONCENTRATIONS
    %{
Imagine that we are simulating a series of three isoprene oxidation experiments, each with
  different NOx conditions. We fill the bag with isoprene, H2O2 (our HOx source), and NO2,
  then turn the lights on and let it run.
Let's assume that initial isoprene and H2O2 are the same for all three runs,
  so we can input them as scalars.
All HoldMe values are set to 0 so that concentrations will evolve.
If you'd rather use CH3ONO as your OH source, add it in below and comment out H2O2.
  Note that this will make more NOx.
    %}
    
    InitConc = {...
        %   names       conc(ppb)           HoldMe
        'LIMONENE'       iLIMONENE                0;
        'HONO'       iHONO                0;
        'NO'         iNO                  0;
        'NO2'       iNO2                  0;...
        'O3'        500                     0;
        'CO'        0                      0; 
        'OH'        0                       0;
        'ttlOA'      iCoa                 0;
        'OAinit'     iCoa                 1;...
        };
    
    seedValues = {...
        'Rp'         Rp                    0;...    %initial seed radius (cm)
        'SA'         SAi                   0;...    %initial seed surface area (cm2/cm3)
        'Vseed'      Vseed                 0;...    %initial seed volume (cm3/cm3)
        'Nseed'      Nseed                 0;...    %initial number of seed particles
        'MWoa_init'   MWoa_init             1;...
        };
    
    
    %% CHEMISTRY
    %{
The ChemFiles input is a cell array of strings
 specifying functions and scripts for the chemical mechanism.
THE FIRST CELL is always a function for generic K-values.
THE SECOND CELL is always a function for J-values (photolysis frequencies).
All other inputs are scripts for mechanisms and sub-mechanisms.
Here we give an example using MCMv3.3.1.  Note that this mechanism was extracted from the MCM website for
the specific set of initial species included above (isoprene and inorganics).
    %}
    
    ChemFiles = {...
        'MCMv331_K(Met)';
        'MCMv331_J(Met,1)'; %Jmethod flag of 0 specifies using "MCM" J-value method.
       % 'MCMv331_Inorg_Isoprene';
        'limonene';
        'limonene_particle';
       % 'Isoprene_DHDHP_Scheme_V3old';
        };
    
    %create ParticleChemFiles list - these will execute after Gas Particle
    %Partitioning module created
    
    ParticleChemFiles= {...
        'limonene_particle'
        };%
    
    %% Vapor Pressures // Molecular Mass // SMILES
    % To run gas-particle partitioning module, you will need to have a
    % saturation vapor concentration ("Cstar") for each specie you want to
    % participate in the gas-particle partitioning. You can set a threshold
    % Cstar above which the species is excluded from gas-particle partitioning.
    % This helps the model run faster by limiting the number of species requiring the more intensive
    % calculations
    
    %If this is the first time running the model, you will need to find and
    %possibly create two files (1) a .mat file containing info on each species
    %"...SpeciesInfo.mat" and (2) an excel file containing the vapor pressures for
    %each species. For the latter you may need to run
    %GetSMILESforVaporPressure.m first which will create a third file (SMILES.txt). These files are usually in
    % ...\Tools\SMILES\.
    
    matfile = 'MCMv331SpeciesInfo.mat'; %full path to the SpeciesInfo file
    Vpfile = 'SMILES_vapPress_Compernolle.xlsx'; %full path to the vapor pressure file
    
    %% DILUTION CONCENTRATIONS
    % We are not diluting the chamber air, so this input is irrelevant (but still necessary).
    
    BkgdConc = {...
        %   names           values
        'DEFAULT'       zeros(1,length(numruns));   %0 for all zeros, 1 to use InitConc
        };
    
 
    
    %% MODEL RUN
    % Now we call the model.
    % Output will be saved in the "SavePath" above and will also be written to the structure S.
    % Let's also throw away the inputs (don't worry, they are saved in the output structure).
    
    runnum = i;
    runtot = length(expts);
    
    if ModelOptions.Do_GP
        GPstruct.seedValues = seedValues;
        GPstruct.gama = gama;
        GPstruct.wall_params = wall_params;
        GPstruct.Vpfile = Vpfile;
        GPstruct.matfile = matfile;

    
    S = F0AM_ModelCore_GP(runnum,runtot,GPstruct,Met,InitConc,ChemFiles,ParticleChemFiles,BkgdConc,ModelOptions);
    clear Met InitConc ChemFiles ParticleChemFiles BkgdConc ModelOptions
    
    else
        
        ModelOptions = rmfield(ModelOptions,{'Do_GP','GPupdateTime','Do_Wall_Partitioning'});
        
        %may need to remove other fields in Met or alter F0AM_ModelCore
        %Met "Fields" section
        rmi = strmatch('kwall_loss',Met(:,1)); rmj = strmatch('Cstar_threshold',Met(:,1)); rm = union(rmi,rmj);
        Met = Met(~ismember(1:length(Met),rm),:);
        display(['Entering ',num2str(runnum),' of ',num2str(runtot)])
        S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);
        clear Met InitConc ChemFiles ParticleChemFiles BkgdConc ModelOptions
    end
   

end

decay = {'HONO';'NO'; 'NO2';'O3'};
OH={'OH'};
PlotConcGroup(decay,S,4,'ptype','line')
PlotConcGroup(OH,S,1,'ptype','line')
PlotConcGroup({'LIMONENE'},S,1,'ptype','line')
PlotConcGroup({'ttlOA'},S,1,'ptype','line')
PlotRates('LIMONENE', S, 3)
