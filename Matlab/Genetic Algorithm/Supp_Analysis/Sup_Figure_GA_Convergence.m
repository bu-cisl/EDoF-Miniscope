%First load the correct dir
addpath('BIG GA Test from Remote Desktop - Supp Figure');

%Save directory
save_dir = 'Sup Figure GA Convergence/';

%What was the sweep params?
Pop = 10:10:100; %change in pop
Gen = 2:2:20; %Change in gen
iter = 10; %Num trials per condition

%First we can use the saved data to get figures
load('Cost_Mean_Table.mat');
load('Cost_STD_Table.mat');
load('Runtime_Mean_Table.mat');
load('Runtime_STD_Table.mat');
load('CPUtime_Mean_Table.mat');
load('CPUtime_STD_Table.mat');

%Note that will probably commented out - I ended this one early so let me
%remove the last two rows until I merge another result to these files
cost_mean_table = cost_mean_table(1:8,:);
cost_std_table = cost_std_table(1:8,:);
cputime_mean_table = cputime_mean_table(1:8,:);
cputime_std_table = cputime_std_table(1:8,:);
runtime_mean_table = runtime_mean_table(1:8,:);
runtime_std_table = runtime_std_table(1:8,:);

%Now do some plotting!
figure
surferror(Gen,Pop(1:8),cost_mean_table,cost_std_table) %No axis labels - manually put later
colormap jet
axis('square')
title('Cost vs Gen & Pop','Fontsize',16)
colorbar

figure
imagesc(Gen,Pop(1:8),cost_mean_table) %No axis labels - manually put later
colormap jet
axis('square')
title('Cost vs Gen & Pop','Fontsize',16)
ylabel('Population','fontsize',14)
xlabel('Generation','fontsize',14)
colorbar

figure
imagesc(Gen,Pop(1:8),cost_std_table./abs(cost_mean_table))
colormap jet
axis('square')
title('Volatility vs Gen & Pop','Fontsize',16)
ylabel('Population','fontsize',14)
xlabel('Generation','fontsize',14)
colorbar

figure
imagesc(Gen,Pop(1:8),log10(cputime_mean_table))
axis('square')
title('Log10(CPU Runtime) vs Gen & Pop','Fontsize',16)
ylabel('Population','fontsize',14)
xlabel('Generation','fontsize',14)
colorbar

%% Lets analyze the convergence parameters - i.e. learned masks!
iter_step = 2; %This is a lot of data, we can sample every iter_step itertions instead of make it easier to digest

%Iterate through all .mat objects representing the points
learned_params = zeros(numel(dir(fullfile('BIG GA Test from Remote Desktop - Supp Figure/Learned Mask', '**', '*.mat')))/iter_step,3);

%While a triple loop is messy it makes it super easy to generate the
%correct name
for i = 1:length(Pop)-2
    for j = 1:length(Gen)
        for k = 1:iter_step:iter
            load(['BIG GA Test from Remote Desktop - Supp Figure/Learned Mask/' ...
                'GA_Opt_Params_Pop_' num2str(Pop(i)) '_Gen_' num2str(Gen(j)) '_iter_' num2str(k) '.mat'])
        learned_params(i+j+ceil(k/iter_step),:) = optvars;
        end
    end
end