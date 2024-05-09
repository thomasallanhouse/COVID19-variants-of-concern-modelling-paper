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

%--------------------------------------------------------------------------
%% NO VOC RUN %%
%--------------------------------------------------------------------------

VOC_introduction_date = 200;

% run the model with no VOC, to get the immunity profile for the Jupyter
% program
clear changed_parameters
% changed_parameters.VOC_imp_size = 0;
changed_parameters.s_VOC = 0.6; % susceptibility to VOC variant for UK recovereds
changed_parameters.VOC_imp_date = 738293+VOC_introduction_date;
% changed_parameters.VOC_imp_size = 10;
% changed_parameters.e_pVOC_scaling = 0; % efficacy of Pfizer vaccine for VOC variant, proportional scaling of resident variants



parameters = make_parameters(changed_parameters);
[t,no_VOC_pop_out,no_VOC_parameters,no_VOC_outputs] = run_simple_vaccines_mex(parameters);

figure; plot(t,no_VOC_outputs.I_UK)
hold on; plot(t,no_VOC_outputs.I_VOC)
% no_VOC_pop_out is no_VOC_pop_out(i,j,k)
% UK VOC V
% S  S  0
% E  E  VA
% I  I  VP
% R  R  VN
% so e.g. no_VOC_pop_out(3,4,4) are people who are infected with UK resident variants, recovered from VOC and vaccinated with new vaccine
% we want: 
% # Proportion of population in each type
% p_vac = 0.7
% p_no_vac = 1- p_vac 
% 
% p_rec = 0.2 # % recovered from previous strains
% p_sus = 1-p_rec # % never infected with any coronavirus
% 
% # Reduced susceptibility for each type based on infection/vaccine immunity
% const_vec = np.array((1., sus_ur, sus_vu, sus_vr))
% prop_vec = np.array((p_sus*p_no_vac, p_rec*p_no_vac, p_sus*p_vac, p_rec*p_vac))
%%%%%%%% Q: this then calculates scale = const_vec*prop_vec but isn't used
%%%%%%%% Q: do we care about p_vac and p_rec on their own or only const_vec
%%%%%%%% and prop_vec?

%% Write required outputs to file to be used by the Python program
% Find the proportions in each group
VOC_intro_pop = no_VOC_pop_out(:,:,:,VOC_introduction_date);
N = sum(VOC_intro_pop,'all');
prop_sus_no_vac = sum(VOC_intro_pop(1,:,1),'all')/N; % proportion not previously infected and not vaccinated
prop_rec_no_vac = sum(VOC_intro_pop(2:4,:,1),'all')/N; % proportion previously infected and not vaccinated
prop_sus_vacc = sum(VOC_intro_pop(1,:,2:3),'all')/N; % proportion not previously infected and vaccinated
%%%%%%%% Q: which vaccination to use? Or work out a mixture?
prop_rec_vacc= sum(VOC_intro_pop(2:4,:,2:3),'all')/N; % proportion previously infected and vaccinated
writematrix([prop_sus_no_vac,prop_rec_no_vac,prop_sus_vacc,prop_rec_vacc],'prop_vec_in.csv')

% find the susceptibility of different groups
sus_ur = parameters.s_VOC; % susceptibility of unvaccinated, previously-infected, against the new strain
sus_vu = parameters.e_pVOC; % (or e_aVOC) susceptibility of vaccinated, not previously infected, against the new strain
sus_vr = min(sus_ur,sus_vu);
writematrix([1,sus_ur, sus_vu, sus_vr],'const_vec_in.csv');

% parameters.s_UK = 0; % susceptibility to UK variant for VOC recovereds
% parameters.s_VOC = 0; % susceptibility to VOC variant for UK recovereds
% 
% parameters.e_aUK = (1-0.65); % proportion of susceptibility remaining after AZ vaccine (1 - efficacy of AZ vaccine for resident variants) 
% parameters.e_pUK = (1-0.75); % proportion of susceptibility remaining after Pfizer vaccine (1 - efficacy of Pfizer vaccine for resident variants) 
% parameters.e_nUK = (1-0.65); % proportion of susceptibility remaining after new vaccine (1 - efficacy of new vaccine for resident variants) 
% 
% parameters.e_aVOC_scaling = 1; % efficacy of AZ vaccine for VOC variant, proportional scaling of resident variants
% parameters.e_pVOC_scaling = 1; % efficacy of Pfizer vaccine for VOC variant, proportional scaling of resident variants
% parameters.e_nVOC_scaling = 1; % efficacy of new vaccine for VOC variant, proportional scaling of resident variants
% 
% 
% % Get effective R for each variant
% [~,ReffUK_no_VOC] = get_Reff(no_VOC_pop_out,no_VOC_parameters);
% 
% % Get infectious temporal profiles from outputs data
% I_UK_no_VOC = no_VOC_outputs.I_UK;
% 
% % Get vaccination coverage data
% vacc_coverage_data = squeeze(sum(no_VOC_pop_out(:,:,2:4,:),[1,2,3]))*100;
