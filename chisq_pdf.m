close all

% Read outputs from Python code
params = csvread('Outputs_for_matlab/FPT_params_R=4.0.csv');
gillespie_samples = csvread('./Gillespie_times_epsilon=0.001.csv');
gillespie_samples_R = gillespie_samples(6, 2:end);

% Take non-zero (i.e. NaN) first-passage times
gillespie_samples_R = gillespie_samples_R(gillespie_samples_R > 0);
growth_rate = params(1);
% Diffusion variance
variance_all = params(2);
% Scaling for eigenvector approximation
evec_scaling = params(3);
% First Passage Time to Z^* = upper_limit
upper_limit = params(4);

time = 0:1:1000;
chisq_integral = zeros(length(time), 1);
xvec = eps:1:20000;
cdf_chisq = zeros(length(time), 1);


% Show comparison between samples and theoretical result
for T = 1:length(time)
    chi_sq_integral = integral(@(x) chisq(x, time(T), growth_rate, variance_all), eps, upper_limit);
    prob =  (1-(chi_sq_integral));
    cdf_chisq(T) = prob;
end

%PDF from CDF
pdf_fpt_chisq = gradient(cdf_chisq);
plot(time, pdf_fpt_chisq)

% plot
hold on 
xlim([0, 200])
ylim([0, 0.04])
samples = inverse_sampling(10000, cdf_chisq, time);
histogram(samples, 'Normalization','pdf')
histogram(gillespie_samples_R, 'Normalization','pdf')


figure(2)
ecdf(samples)
hold on
ecdf(gillespie_samples_R-4)
plot(time, cdf_chisq)
xlim([0, 200])

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
        sample(i) = t(min(find(cdf>=runiform(i))));
    end
    
end