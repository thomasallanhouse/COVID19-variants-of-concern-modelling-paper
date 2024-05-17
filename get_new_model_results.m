%% get_new_model_results.m:
% Script to run the simple vaccines model using the jupyter code
%--------------------------------------------------------------------------

% Changelog
% put change_days and R_changes_UK_without_immunity to single values so we
% don't have a relaxation roadmap
clear

%% Set global flag variables
make_mex_flag = true;

%% Make mex: run this the first time to make the mex file
if make_mex_flag == true
    %%
    clear changed_parameters
    changed_parameters.VOC_imp_size = 0;
    parameters = make_parameters(changed_parameters);
    codegen run_simple_vaccines -args {parameters}
end

%% Run the model with a VOC introduction date after the resident wave has finished
% to get the immunity profile for the Jupyter program, and then something
% to compare to after that
VOC_introduction_date = 200;

clear changed_parameters
% changed_parameters.s_VOC = 0.6; % susceptibility to VOC variant for UK recovereds
changed_parameters.VOC_imp_date = datenum(2021,5,17)+VOC_introduction_date;
changed_parameters.s_VOC = 0.6; % susceptibility of unvaccinated, previously-infected, against the new strain
changed_parameters.e_pVOC = 0.4; % (or e_aVOC) susceptibility of vaccinated, not previously infected, against the new strain
changed_parameters.e_aVOC = 0.4;
% changed_parameters.VOC_imp_size = 10;
% changed_parameters.e_pVOC_scaling = 0; % efficacy of Pfizer vaccine for VOC variant, proportional scaling of resident variants

parameters = make_parameters(changed_parameters);
% parameters.e_aVOC = parameters.e_pVOC; % Make the vaccine efficacies the same, for simplicity
% this gives the fully population output: for_jupyter_pop_out(UK,VOC,Vaccine)
% so e.g. for_jupyter_outputs(3,4,4) are people who are infected with UK resident variants, recovered from VOC and vaccinated with new vaccine
% it also saves the inputted parameters and for_jupyter_outputs for ease of plotting
[t,for_jupyter_pop_out,for_jupyter_parameters,for_jupyter_outputs] = run_simple_vaccines_mex(parameters);

figure; plot(t,for_jupyter_outputs.I_UK)
hold on; plot(t,for_jupyter_outputs.I_VOC)

%% Write required outputs to file to be used by the Python program
% Find the proportions in each group
VOC_intro_pop = for_jupyter_pop_out(:,:,:,VOC_introduction_date);

% for vaccinated groups we add up AZ and Pfizer (having put the efficacies
% to be the same)
prop_sus_no_vac = sum(VOC_intro_pop(1,:,1),'all'); % proportion not previously infected and not vaccinated
prop_rec_no_vac = sum(VOC_intro_pop(2:4,:,1),'all'); % proportion previously infected and not vaccinated
prop_sus_vacc = sum(VOC_intro_pop(1,:,2:3),'all'); % proportion not previously infected and vaccinated
prop_rec_vacc= sum(VOC_intro_pop(2:4,:,2:3),'all'); % proportion previously infected and vaccinated
writematrix([prop_sus_no_vac,prop_rec_no_vac,prop_sus_vacc,prop_rec_vacc],'prop_vec_in.csv')

% find the susceptibility of different groups
sus_ur = parameters.s_VOC; % susceptibility of unvaccinated, previously-infected, against the new strain
sus_vu = 1-parameters.e_pVOC; % (or e_aVOC) susceptibility of vaccinated, not previously infected, against the new strain
sus_vr = min(sus_ur,sus_vu);
writematrix([1,sus_ur, sus_vu, sus_vr],'const_vec_in.csv');

%% Now we stop and run the Jupyter model in multitype_matlab_outputs.ipynb
% to obtain Outputs_for_matlab/FPT_params_R....csv
% Either put the path to your python distribution here and run directly
% from matlab, or else go to the python script and run it there
python_path = 'C:\Users\dysonl\anaconda3\';
system([python_path,'python run_multitype_matlab_outputs.py'])

%% and then come back here to use them
params = readmatrix('Outputs_for_matlab/FPT_params_R=4.0.csv');

% Rename the parameters for use in generating samples of the first passage time: growth rate
growth_rate = params(1); 
% Diffusion variance
variance_all = params(2);
% Scaling for eigenvector approximation
evec_scaling = params(3);
% First Passage Time to Z^* = upper_limit
upper_limit = params(4);
% beta for the VOC
VOC_beta = params(5);

time = 0:1:1000;
cdf_chisq = zeros(length(time), 1);

% Generate the cdf
for T = 2:length(time)
    cdf_chisq(T) = 1- integral(@(x) chisq(x, time(T), growth_rate, variance_all), eps, upper_limit);
end

first_passage_times = inverse_sampling(1000, cdf_chisq, time);

% Then we also need the eigenvector so we know what states to start with
% states = (unvac not-previously infected, unvac previously infected, vac not-prev, vac prev-inf)
% evector is exposed states followed by infectious states
evector = readmatrix('Outputs_for_matlab\dominant_eigenvector_R=4.0.csv');
evector = evector/parameters.UK_popn_size; % not sure this is right
evector = num2cell(evector);

% parameters.VOC_imp_distribution should be in the form
% (Resident variants, VOC, vaccination) = (S/R, E/I, U/A/P/N)
VOC_imp_distribution = zeros(2,2,4); % for using an initial VOC prevalence distributed between classes
S = 1; R = 2;
E = 1; I = 2;
not_vac = 1; vac = 2;
% first four entries in the evector are for VOC exposed
[VOC_imp_distribution(S,E,not_vac),VOC_imp_distribution(R,E,not_vac),VOC_imp_distribution(S,E,vac),VOC_imp_distribution(R,E,vac)] = evector{1:4};
% and the second four are VOC infectious
[VOC_imp_distribution(S,I,not_vac),VOC_imp_distribution(R,I,not_vac),VOC_imp_distribution(S,I,vac),VOC_imp_distribution(R,I,vac)] = evector{5:8};

%%
clear changed_parameters
% changed_parameters.s_VOC = 1-0.4; % susceptibility to VOC variant for UK recovereds
save_outputs = zeros(366,length(first_passage_times));
changed_parameters.specify_distribution = true; % by default we don't use the distribution
changed_parameters.VOC_imp_distribution = VOC_imp_distribution;
changed_parameters.beta_VOC_changes = VOC_beta*ones(1,5);
for i=1:length(first_passage_times)
    changed_parameters.VOC_imp_date = datenum(2021,5,17)+VOC_introduction_date+first_passage_times(i);
    changed_parameters.date1 = changed_parameters.VOC_imp_date;
    parameters = make_parameters(changed_parameters);
    parameters.e_aVOC = parameters.e_pVOC; % Make the vaccine efficacies the same, for simplicity
    [~,~,~,VOCintro_outputs] = run_simple_vaccines_mex(parameters);
    save_outputs(:,i) = VOCintro_outputs.I_VOC;
    save_dates(:,i) = VOCintro_outputs.dates;
end

%%
figure; plot(for_jupyter_outputs.dates,for_jupyter_outputs.I_UK*for_jupyter_parameters.UK_popn_size)
hold on
plot(for_jupyter_outputs.dates,for_jupyter_outputs.I_VOC*for_jupyter_parameters.UK_popn_size)
plot(save_dates(:,1:100),save_outputs(:,1:100)*parameters.UK_popn_size)

% PDF of non-central chi**2 with 0 degrees of freedom
function chisq_pdf = chisq(x, t, growth_rate, variance_all)
    
   x_scale = 2.*growth_rate.*x./(((variance_all./2)).*(exp(growth_rate*t) - 1));
   
   lamb = 2*growth_rate*exp(growth_rate*t)./((variance_all./2)*(exp(growth_rate*t) - 1));
   chisq_pdf = growth_rate./((variance_all./2).*(exp(growth_rate*t) - 1)) .* sqrt(exp(growth_rate*t)./x) .* exp(- 1./2 * (lamb + x_scale)) .* besseli(1, sqrt(x_scale*lamb)) ./ ((1-exp(-lamb./2)));
end


% Sample from approximate distribution given by chisq_pdf
function sample = inverse_sampling(n, cdf, t)
    runiform = unifrnd(0, 1, n);
    sample = zeros(n, 1);
    for i = 1:n
        sample(i) = t(find(cdf>=runiform(i),1,'first'));
    end
    
end