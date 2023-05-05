function PlotConcMass(M2Plot, S, n2plot)
% function PlotConcMass(M2Plot, S) 
% Generates a time series of concentrations for all species at specific masses.
% INPUTS:
% M2Plot: list of masses to plot (includes particle phase I think)
% S: structure of model outputs. Must contain the following fields:
%    Time, Cnames, Conc, MW_out
% n2plot: number of individual masses to plot.
%         The first n2plot masses in the structure will be plotted.
%         The rest will be grouped into an "other" category.
%         If this is empty or not specified, all masses will be plotted.



% maybe add later!!
% varargin: One can  specify several options as name-value pairs:
% 
%           PlotConcGroup(...,'ptype',value)
%               Indicates type of plot.
%               Values can include 'fill', 'bar' or 'line'.
%               Default: 'fill'
%
%           PlotConcGroup(...,'unit',value)
%               Changes the concentration unit.
%               Value may be 'ppb', 'ppt', or 'percc' (the last is number density).
%               Default: 'ppb'
% 
%           PlotConcGroup(...,'scale',value)
%               Specifies a multiplier for concentrations.
%               Value is a scalar. For example, 1E-6 might be used when plotting OH number density.
%               SPECIAL CASE: setting this to 0 causes values to be normalized by sum.
%               Default: 1
%
%           PlotConcGroup(...,'name',value)
%               Specifies the name of your group.
%               Value is a string, used for the y axis label
%               Default: '[X]'
%
%           PlotConcGroup(...,'sortem',value)
%               Flag for whether or not to sort the group.
%               Value is a scalar, 0 (no) or 1 (yes).
%               Default: 1
%
% OUTPUTS are optional.
% X2plot: matrix of species plotted.
% Xnames: names of species plotted.
% if desired, you can use the "breakin" command to make variables from these.
%
% 20120724 GMW  Creation Date.
% 20151102 GMW  Modified for F0AMv3, including addition of name-value option inputs.
%               Also added outputs.
% 20160321 GMW  Removed option to input Sp2plot as an index, for consistency w/other function.
%               Added ability to input sub-groups a la PlotReactivity.
%               Added "sortem" option.
% 20230211 ISA  adapted PlotConcGroup to plot groups at masses

%%%%%INPUT CHECKING%%%%%
if ~iscell(M2Plot)
    error('PlotConcGroup: input "M2Plot" must be a cell array.')
end

if nargin<3 || isempty(n2plot)
    n2plot = length(M2Plot);
end


%%%%%OPTIONS%%%%%
varInfo = {...
    %name       %default    %valid
    'unit'      'ppb'       {'ppb','ppt','percc','ppbv','pptv'};...
    'scale'     1           [];...
    'ptype'     'fill'      {'fill','bar','line'};...
    'name'      '[X]'       [];...
    'sortem'    1           [0 1];...
    };
%ParsePairs(varargin,varInfo);

%only gas phase


%create sub-group (SpSubgroup) for each mass that is then saved in a list called Sp2plot
Sp2plot = {};
for m = M2Plot
    SpSubgroup = {};
    for i = 1:length(S.MW_out)
        if m{1}==floor(S.MW_out(i))
            cname=S.Cnames(i);
            SpSubgroup= [SpSubgroup, cname{1}];
            %disp(SpSubgroup);
        end  
    end
    if isempty(SpSubgroup)
        idx = find([M2Plot{:}] == m{1});
        M2Plot(idx) = [];
        Sp2plot = Sp2plot;
    else
        Sp2plot= [Sp2plot, {SpSubgroup}];
    end
end

disp(Sp2plot{1});


%%%%SUM SUB-GROUPS%%%%%
for i=1:length(Sp2plot)
    Sp = Sp2plot{i};
    if iscell(Sp)
        group = ExtractSpecies(Sp,S);
        S.Conc.(Sp{1}) = sum(breakin(group),2); %temporary summed sub-group
        Sp2plot{i} = Sp{1};
    end
end
S.Cnames = fieldnames(S.Conc);

%%%%%GET CONCENTRATIONS%%%%%
Conc = ExtractSpecies(Sp2plot,S, 1); % get concs
[Xmat,Xnames] = breakin(Conc);
Time = S.Time;
XnamesAll = Xnames;

X2plot = Xmat(:,1:n2plot); %assumes already sorted

% create mass legend in order of sorting
massLegend = {};
for x = 1:length(XnamesAll)
    for i = 1:length(S.Cnames)
        cname =S.Cnames{i};
        if XnamesAll{x} == string(cname)
            idx = i;
            massLegend{x}=floor(S.MW_out(idx));
        end
    end
end

%account for "other" in mass legend
massLegendFinal = [string(massLegend(1:n2plot)), 'other'];


%group "other"
if n2plot<length(Xnames)
    Xother = Xmat(:,n2plot+1:end);
    Xothernames = Xnames(n2plot+1:end);
    
    X2plot(:,end+1) = sum(Xother,2);
    Xnames = [Xnames(1:n2plot); 'Other'];
    
    [XoMax,iMax] = max(max(Xother)); %identify largest contribution in "other"
    disp(['PlotConcGroup: ' ])
    disp(['  Number of masses in "other" = ' num2str(length(Xothernames))])
    disp(['  Largest contribution in "other" is ' string(massLegend(n2plot+1)) ' at ' num2str(XoMax) 'ppb'])
end


% create mass legend
massLegend = {};
for x = 1:length(XnamesAll)
    for i = 1:length(S.Cnames)
        cname =S.Cnames{i};
        if XnamesAll{x} == string(cname)
            idx = i;
            massLegend{x}=floor(S.MW_out(idx));
        end
    end
end

massLegendFinal = [string(massLegend(1:n2plot)), 'other'];

% flip it
% X2plot = fliplr(X2plot);
% Xnames = flipud(Xnames);

%%%%%PLOT%%%%%
figure

plot(Time,X2plot,'LineWidth',2)

%Plot decorations
load fillcolors.mat
colormap(fillcolors);

legend(massLegendFinal,'Location','NorthEast')
ylabel('ppb')
xlabel('Model Time')
xlim([min(Time) max(Time)])
purtyPlot

%%%%%PLOT 2%%%%%
figure

plot(Time,log(X2plot),'LineWidth',2)

%Plot decorations
load fillcolors.mat
colormap(fillcolors);

legend(massLegendFinal,'Location','NorthEast')
ylabel('ln(ppb)')
xlabel('Model Time')
xlim([min(Time) max(Time)])
purtyPlot

%%%%%OUTPUT%%%%%
if nargout==0, clear X2plot Xnames; end
    

