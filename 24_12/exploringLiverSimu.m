setup

% dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\liver_RED';

% dataDir = 'P:\smerino\simulation_acs\rf_data\liver';
% targetFiles = dir(fullfile(dataDir,"rf_qus_livernew202410_AC_test.mat"));

% dataDir = 'P:\smerino\simulation_acs\rf_data\24_12_02_liver';
% targetFiles = dir(fullfile(dataDir,"rf_qus_livernew.mat"));

dataDir = 'P:\smerino\simulation_acs\rf_data\24_12_04_liver';
targetFiles = dir(fullfile(dataDir,"rf*.mat"));

% refDir = fullfile(dataDir,'ref');
refDir = 'P:\smerino\simulation_acs\rf_data\24_12_02_liver\ref';
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = fullfile(dataDir,'test','mySim');
mkdir(resultsDir)
%%
c0 = 1540;
freqL = 3.5e6; freqH = 7.5e6;
% freqL = 3e6; freqH = 8e6;
wl = 2*c0/(freqL + freqH);

blockParams.xInf = -2;
blockParams.xSup = 2;
blockParams.zInf = 2.9;
blockParams.zSup = 5.5;
blockParams.blocksize = [20 20]*wl;
blockParams.overlap = 0.8;

NptodB = log10(exp(1))*20;
iAcq = 2;
%%
for iAcq = 2:length(targetFiles)
    %% Loading sample
    fprintf("Simulation no. %i, %s\n",iAcq,targetFiles(iAcq).name);
    out = load(fullfile(dataDir,targetFiles(iAcq).name));
    xBm = out.x*1e2; % [cm]
    zBm = out.z'*1e2; % [cm]
    c0 = 1540;
    sam1 = out.rf(:,:);
    fs = out.fs;
    kgrid = out.kgrid;
    attenuation_map = out.attenuation_map;
    liverAcs = attenuation_map(end,end);

    % Plot region of interest B-mode image
    bMode = db(hilbert(sam1));
    bMode = bMode - max(bMode(:));

    % get Spectra
    [Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,blockParams);

    % Plotting spectra
    spectrumSamzf = db(squeeze(mean(Sp/2+Sd/2, 2)));
    spectrumSamzf = spectrumSamzf - max(spectrumSamzf,[],2);
    figure,
    imagesc(f,z_ACS,spectrumSamzf, [-50 0]),
    ax = gca; ax.YDir = 'reverse';
    hold on
    xline(freqL/1e6, 'k--'), xline(freqH/1e6, 'k--')
    hold off
    colorbar
    xlim([0 12])
    xlabel('f [MHz]')
    ylabel('z [cm]')
    title('Sample power spectrum by depth')

    %% Generating Diffraction compensation
    % Generating references
    clear att_ref_map rfRef

    for ff = 1:length(refFiles)
        out = load(fullfile(refDir,refFiles(ff).name));
        rfRef(:,:,ff) = out.rf(:,:,end);
    end
    alpha0Ref = double(mean(out.attenuation_map(:))); gammaRef = 1;
    att_ref_map(1,1,:) = alpha0Ref*f.^gammaRef/NptodB;

    [SpRef,SdRef,~,~,~] = getSpectrum(rfRef,xBm,zBm,fs,blockParams);

    % Plotting spectra
    spectrumRefzf = db(squeeze(mean(SpRef/2+SdRef/2, 2)));
    spectrumRefzf = spectrumRefzf - max(spectrumRefzf,[],2);
    figure,
    imagesc(f,z_ACS,spectrumRefzf, [-50 0]),
    ax = gca; ax.YDir = 'reverse';
    hold on
    xline(freqL/1e6, 'k--'), xline(freqH/1e6, 'k--')
    hold off
    colorbar
    xlim([0 12])
    xlabel('f [MHz]')
    ylabel('z [cm]')
    title('Reference power spectrum by depth')

    %% Setting up
    L = (z_ACS(2) - z_ACS(1))/(1 - blockParams.overlap)/2;   % (cm)
    sld = log(Sp) - log(Sd);
    sldRef = log(SpRef) - log(SdRef);
    compensation = sldRef - 4*L*att_ref_map;

    range = f>freqL/1e6 & f<freqH/1e6;
    b = sld(:,:,range) - compensation(:,:,range);
    ufr = f(range);

    [m,n,p] = size(b);
    A1 = kron( 4*L*ufr , speye(m*n) );
    A2 = kron( ones(size(ufr)) , speye(m*n) );
    mask = ones(m,n,p);
    tol = 1e-3;

    %% SLD by depth
    nLines = 4;
    depthPoints = floor(m/nLines);
    figure, hold on
    for ii = 1:nLines
        line = mean(b(ii*depthPoints-depthPoints+1:ii*depthPoints,:,:),[1 2]);
        zline = mean(z_ACS(ii*depthPoints-depthPoints+1:ii*depthPoints));
        plot(ufr,squeeze(line*NptodB/4/L), 'LineWidth',2)
        leg{ii} = "z = "+round(zline,2)+"cm";
    end
    plot(ufr,liverAcs*ufr, 'k--')
    hold off
    grid on
    leg{ii+1} = 'Ideal';
    legend(leg, 'Location','northwest')
    xlabel('Frequency [MHz]')
    ylabel('Attenuation [dB/cm]')


    %% SLD
    A = [A1,A2];
    [u,~] = cgs(A'*A,A'*b(:),1e-6,200);
    BSLD = (reshape(u(1:m*n)*NptodB,m,n));
    CSLD = (reshape(u(1+m*n:end)*NptodB,m,n));

    figure,
    plot(z_ACS,mean(BSLD,2))
    grid on
    xlabel('Depth [cm]')
    ylabel('ACS [dB/cm/MHz]')
    yline(liverAcs,'k--')
    title('SLD, ACS by depth')
    ylim([-1 3])


    %% RSLD-TV
    muRsld = 10^1;

    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
    toc
    BRSLD = (reshape(Bn*NptodB,m,n));
    CRSLD = (reshape(Cn*NptodB,m,n));

    % Creating masks and ideal map
    [Xq,Zq] = meshgrid(xBm,zBm);
    attIdeal = interp2(kgrid.y*100, kgrid.x*100 - kgrid.x(1)*100, attenuation_map, Xq, Zq);
    roi = ones(size(bMode));

    % Plotting constants
    dynRange = [-60,0];
    attRange = [0.2,1.6];
    bsRange = [-5,5];

    figure('Units','centimeters', 'Position',[5 5 24 8]);
    tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

    t1 = nexttile;
    imagesc(xBm,zBm,bMode,dynRange)
    axis image
    c = colorbar(t1, 'westoutside');
    c.Label.String = 'dB';
    title('B-mode')
    ylabel('Axial [cm]')
    xlabel('Lateral [cm]')

    t2 = nexttile;
    imagesc(xBm,zBm,attIdeal,attRange)
    xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
    % axis equal
    % xlim([x_ACS(1) x_ACS(end)]),
    % ylim([z_ACS(1) z_ACS(end)]),
    axis image
    title('Ideal')

    t3 = nexttile;
    [~,hB,hColor] = imOverlayInterp(bMode,BRSLD,dynRange,attRange,1,...
        x_ACS,z_ACS,roi,xBm,zBm);
    title('RSLD, ACS')
    % colorbar off
    xlabel('Lateral [cm]')

    t4 = nexttile;
    [~,hB,hColor] = imOverlayInterp(bMode,CRSLD,dynRange,bsRange,1,...
        x_ACS,z_ACS,roi,xBm,zBm);
    title('RSLD, \Delta BSC')
    % colorbar off
    xlabel('Lateral [cm]')

    colormap(t1,gray)
    colormap(t2,turbo)
    colormap(t3,turbo)
    colormap(t4,parula)

    %%
    fit = [z_ACS,ones(m,1)]\mean(BRSLD,2);
    figure,
    plot(z_ACS,mean(BRSLD,2))
    grid on
    xlabel('Depth [cm]')
    ylabel('ACS [dB/cm/MHz]')
    yline(liverAcs,'k--')
    title(sprintf("RSLD, ACS = %.2fz + %.2f",fit(1),fit(2)))
    ylim([0 2])
    hold on,
    plot(z_ACS,fit(1)*z_ACS + fit(2), 'r--')

    %%
    save_all_figures_to_directory(resultsDir,char("sam"+iAcq+"fig"))
    close all
    
end

