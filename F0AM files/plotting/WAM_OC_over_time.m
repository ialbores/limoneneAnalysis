function k = WAM_OC_over_time(Sin,volume)
%Sin is the FOAM output structure to use
%volume is the volume of your bag (in m^3)
%contributors to predicted SOA mass conc.

S = Sin;

% Convert structure to cell
cell = struct2cell(S.Conc);
% Construct an array
alconc = [cell{:}];

%get particle concentrations
ps = (alconc(:, S.iOA))';
ps(isnan(ps))=0;

%convert to mol
%assuming MW_out is correct for finding molar mass
molar_mass = S.particle.mass;
for i = 1:length(S.iOA) 
    for j = 1:length(ps)
        mol(i, j) = ps(i,j)/1e6*volume/molar_mass(i);
    end
end

%get particle names
particle_names = S.particle.names;
%remove p from end of particle names
for i = 1:length(particle_names)
    particle_names{i} = particle_names{i}(1:end-1);
end
%load all MCM names
MCMnames = load('/Users/isabel/Documents/MATLAB/F0AM-WAM-2020-distribute/Cyclic/SpeciesInfoANE.mat');
MCMnames = struct2cell(MCMnames);
particle_names = convertCharsToStrings(particle_names);
MCMnames = convertCharsToStrings(MCMnames);
%find indices of particle names in MCM names
[~,~,iMCMnames] = intersect(particle_names, MCMnames{1});
%find particle species in smiles using MCM indices
smiles = load('/Users/isabel/Documents/MATLAB/F0AM-WAM-2020-distribute/Cyclic/SpeciesInfoANE.mat', 'SMILES');
smiles = struct2cell(smiles);
smiles= smiles{1};
particle_smiles = smiles(iMCMnames);

%use chem trans to get formulas from smiles
formula = chemtrans(particle_smiles, 'formula');
formula =formula(2:end,2);

%use molecular formulas to find mol O and mol C in aerosol

for i = 1:length(formula)
    c_number(i) = sscanf(formula{i},'C%d');
    h_number(i) = sscanf(formula{i},append('C',string(c_number(i)),'H%d'));
    o_test = extractAfter(formula{i}, 'O');
    if o_test
        o_number(i) = str2double(o_test);
    else 
        o_number(i) = str2double('1');
    end
    n_test = contains(formula{i}, 'N');
    no_test = contains(formula{i}, 'NO');
    if n_test 
        if no_test
            n_number(i) = 1;
        else 
            n_number(i) = sscanf(formula{i},append('C',string(c_number(i)),'H',string(h_number(i)),'N%d'));%extractBetween(formula{i}, 'N', 'O');
        end
    else
        n_number(i) = 0;
    end
    %correction to take out oxygens from nitrate groups, can comment out if
    %you want to include those!
    o_number(i) = o_number(i) -3*n_number(i);
end

%sum mol O and mol C
for i=1:length(c_number)
    for j = 1:length(ps)
        mol_c_by_cpd(i, j) = mol(i, j) * c_number(i);
        mol_h_by_cpd(i, j) = mol(i, j) * h_number(i);
        mol_o_by_cpd(i, j) = mol(i, j) * o_number(i);
        mol_n_by_cpd(i, j) = mol(i, j) * n_number(i);
    end
end


total_mol_c = sum(mol_c_by_cpd);
total_mol_h = sum(mol_h_by_cpd);
total_mol_o = sum(mol_o_by_cpd);
total_mol_n = sum(mol_n_by_cpd);

%take ratio
o_c_ratio = total_mol_o./total_mol_c;
h_c_ratio = total_mol_h./total_mol_c;
n_c_ratio = total_mol_n./total_mol_c;


avg_o_c_ratio = mean(o_c_ratio);
avg_h_c_ratio = mean(h_c_ratio);
avg_n_c_ratio = mean(n_c_ratio);


 display(avg_o_c_ratio);
 display(avg_h_c_ratio);
 display(avg_n_c_ratio);


%can choose to plot
plot(o_c_ratio);
hold on 
plot(h_c_ratio);
plot(n_c_ratio);
legend('O/C','H/C','N/C')
axis([0 360 0 3]); % change y-axis value depending on time of run


end

