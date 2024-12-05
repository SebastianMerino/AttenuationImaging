clear, clc
BaseDir = 'P:\smerino\simulation_acs\rf_data\24_12_04_liver';

%% Sample 9 -> no effect
load(fullfile(BaseDir,'dens_08_2024b.mat'))

[M,N] = size(liver11);
liver_map = imresize(liver11,[M/2,N], 'nearest');
liver_map = [liver_map; 6*ones(M/2,N)];

figure, tiledlayout(1,2)
nexttile, imagesc(liver11)
axis image
nexttile, imagesc(liver_map)
axis image

save(fullfile(BaseDir,'dens_09_2024b.mat'),"liver_map")


load(fullfile(BaseDir,'att_08_2024b.mat'))

[M,N] = size(liver11);
liver_map = imresize(liver11,[M/2,N], 'nearest');
liver_map = [liver_map; liver11(end,end)*ones(M/2,N)];

figure, tiledlayout(1,2)
nexttile, imagesc(liver11)
axis image
nexttile, imagesc(liver_map)
axis image

save(fullfile(BaseDir,'att_09_2024b.mat'),"liver_map")

%% Sample 10 -> let's see
load(fullfile(BaseDir,'dens_08_2024b.mat'))
[M,N] = size(liver11);
liver_map = [liver11(1:370,:); liver11(615:end,:)];
liver_map = [liver_map; liver11(end,end)*ones(M - size(liver_map,1),N)];
figure, tiledlayout(1,2)
nexttile, imagesc(liver11)
axis image
nexttile, imagesc(liver_map)
axis image
save(fullfile(BaseDir,'dens_10_2024b.mat'),"liver_map")

load(fullfile(BaseDir,'att_08_2024b.mat'))
[M,N] = size(liver11);
liver_map = [liver11(1:370,:); liver11(615:end,:)];
liver_map = [liver_map; liver11(end,end)*ones(M - size(liver_map,1),N)];
figure, tiledlayout(1,2)
nexttile, imagesc(liver11)
axis image
nexttile, imagesc(liver_map)
axis image
save(fullfile(BaseDir,'att_10_2024b.mat'),"liver_map")

%% Sample 12 -> let's see
in = load(fullfile(BaseDir,'dens_10_2024b.mat'));
[M,N] = size(in.liver_map);
liver_map = imresize(in.liver_map,[M*1.5,N], 'nearest');
liver_map = liver_map(1:M,:);
figure, tiledlayout(1,2)
nexttile, imagesc(in.liver_map)
axis image
nexttile, imagesc(liver_map)
axis image
save(fullfile(BaseDir,'dens_11_2024b.mat'),"liver_map")

in = load(fullfile(BaseDir,'att_10_2024b.mat'));
[M,N] = size(in.liver_map);
liver_map = imresize(in.liver_map,[M*1.5,N], 'nearest');
liver_map = liver_map(1:M,:);
figure, tiledlayout(1,2)
nexttile, imagesc(in.liver_map)
axis image
nexttile, imagesc(liver_map)
axis image
save(fullfile(BaseDir,'att_11_2024b.mat'),"liver_map")

%% Sample 12 -> let's see
load(fullfile(BaseDir,'dens_08_2024b.mat'))
liver_map = liver11;
save(fullfile(BaseDir,'dens_12_2024b.mat'),"liver_map")

load(fullfile(BaseDir,'att_08_2024b.mat'))
liver_map = liver11;
save(fullfile(BaseDir,'att_12_2024b.mat'),"liver_map")
