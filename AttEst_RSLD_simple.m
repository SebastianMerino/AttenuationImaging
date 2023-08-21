%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pontifica Universidad Católica del Perú
% Attenuation Estimator RSLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code starts here
function [] = AttEst_RSLD_simple()

clear all; close all; clc;

addpath('./functions_att');

set(0,'DefaultTextFontName','Arial');
set(0,'DefaultTextFontSize',12);
font = 42;
%db_size = [10 20 30 40 50];   % CAMBIAR
%font = 42;
db_size = [20];   % CAMBIAR

frame_list = 1;
angle_list = 0;

% LOOPING FOR EACH FRAME AND FOR EACH SS (?) 
for ff = 1:length(frame_list)
    frame = frame_list(ff);
    angle = angle_list(ff);
    ACS_estimation = 1;
    
    for ss = 1:length(db_size)
        blocksize = db_size(ss);
        % Directory wherefigures ans result will be saved
        output_dir = [pwd,'\figures\frame',num2str(ff),'\FIG',num2str(ss)];

        rf_class = 102;
        homogeneous = 0;
        ACS_estimation  = 1;
        winsize         = 0.5;
        window_type     = 5;
        saran_layer     = 0;
        spectra_estimation = 1;
        type_inclusion  = 1;

        % parameters
        load parameters.mat
        mu = mu_values;
        %att_ref_dB = attenuation_reference_dB;
        bottom = ACS_range_maps(1);
        top = ACS_range_maps(2);
        x_inf = lateral_limits(1);
        x_sup = lateral_limits(2);
        z_inf = axial_limits(1);
        z_sup = axial_limits(2);
        c1x = center_inclusion(1);
        c1y = center_inclusion(2);
        c2x = center_background(1);
        c2y = center_background(2);
        db_size = data_block_sizes;
        dyn_range = dynamic_range_bmode;
        freq_L = frequency_limits(1);
        freq_H = frequency_limits(2);
        ground_truth = ground_truth_values;
        overlap_pc = overlap_percentage;

        % Reference phantom method by default
        load sam1.mat rf fs c0 x z
        sam1 = rf;  % sample 1 is RF data

        ind_x =  overall_limits(1) <= x & x <= overall_limits(2); % 8th
        %ind_z = 0.042<= z & z <= 0.123; % 7th
        x=x(ind_x);
        ind_z =  overall_limits(3) <= z & z <= overall_limits(4); % 8th
        %ind_z = 0.042<= z & z <= 0.123; % 7th
        z=z(ind_z);
        sam1 = sam1(ind_z,ind_x);

        %% Plotting B-mode

        % Lateral (x) and axial (z) deltas
        dx=(x(end)-x(1))/(length(x)-1);
        dz=(z(end)-z(1))/(length(z)-1);

        % x and z at the starting time
        x_ori = x;
        z_ori = z;

        % x and z in [cm]
        x = x*1e2;
        z = z*1e2;

        % Plot B-mode image
        Im = abs(hilbert(sam1));
        Im_db = 20*log10(Im/max(Im(:)));
        Im_db_ori = Im_db;
        figure(1); set(1,'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'FontSize',font);
        imagesc(x,z,Im_db); axis image; colormap gray; caxis([0-dyn_range 0]);
        hb1 = colorbar;
        ylabel(hb1,'dB','FontSize', font)
        %title('B-mode image');
        xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
        set(gca,'FontSize',font);
                

        %% Region of interest selection

        % Cut until the limits for ACS estimation
        ind_x = x_inf <= x & x <= x_sup;
        x = x(ind_x);
        ind_z = z_inf <= z & z <= z_sup;
        z = z(ind_z);
        sam1 = sam1(ind_z,ind_x);

        freq_L = freq_L*1e6;   % (Hz)
        freq_H = freq_H*1e6;   % (Hz)

        wl = c0/mean([freq_L freq_H]);   % Wavelength (m)
        disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
            num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);

        % Number of repetitions of datablocks within a single
        % datablock
        rpt = 1/(1-overlap_pc);   % r = rep = 1/(1-overlap)
        rpt = round(rpt);
        overlap = 1 - (1/rpt);

        % Part without lateral overlap
        wx  = round((blocksize*wl)/(dx*rpt));
        % Number of lines laterally = repetitions * block without repetition
        nx  = rpt*wx;
        new_blocksize = round(nx*dx/(wl));

        % RF data columns
        L2   = size(sam1,2);
        % Blocksize colums
        ncol = floor((L2-(rpt-1)*wx)/wx);
        sam1 = sam1(:,1:wx*(ncol+rpt-1));
        % Actual rf data columns
        L2 = size(sam1,2);
        x  = x(1:L2);

        xi = 1;
        xf = L2;
        x0 = (xi:wx:xf+1-nx);
        x_ACS = x(1,x0+round(nx/2));
        nx_ACS = nx*dx;

        n  = length(x0);

        wz = floor(nx*dx/(dz*rpt));
        nz = rpt*wz;

        % winsize: Percentage of window (max 0.5)
        % nw: Samples of each window axially
        nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;
        L = (nz - nw)*dz*100;   % (cm)

        NFFT = 2^(nextpow2(nw)+2);
        band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
        rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
        f  = band(rang)*1e-6; % [MHz]
        L3 = length(rang);
        p = L3;

        L1   = size(sam1,1);   % RF data: rows
        nrow = floor((L1-(rpt-1)*wz)/wz);        % Blocksize rows
        sam1 = sam1(1:wz*(nrow+rpt-1),:);
        L1   = size(sam1,1);   % RF data: rows
        z    = z(1:L1);

        zi = 1;
        zf = L1;
        z0 = (zi:wz:zf+1-nz);
        m  = length(z0);

        z_ACS = z(z0+round(nz/2));
        nz_ACS = nz*dz;

        z0p = z0 + (nw-1)/2;            % Proximal coordinate
        z0d = z0 + (nz-1) - (nw-1)/2;   % Distal coordinate

        %% Plot region of interest B-mode image
        %db_lim = input('Dyn. range (ejm: 50): ');   % Frame
        Im=abs(hilbert(sam1));   % envelope calculation
        Im_db=20*log10(Im/max(Im(:)));   % log scale
        figure(2); set(2,'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'FontSize',font);
        imagesc(x,z,Im_db); axis image; colormap gray; caxis([0-dyn_range 0]);
        hb2=colorbar; ylabel(hb2,'dB','FontSize', font)
        %title('B-mode image');
        xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
        set(gca,'FontSize',font);

        %% CALCULATE SPECTRA

        disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
            num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
        disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
        disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
        disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
        disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);


        %bottom = input('Bottom ACS (ejm: 0): ');
        %top = input('Top ACS (ejm: 1.5): ');

        %window_type = input('Windowing type (Tukey-0.25 (5), Hamming (6), Boxcar (7)): ');
        switch window_type
            case 5
                windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
            case 6
                windowing = hamming(nw);   % Hamming
            case 7
                windowing = rectwin(nw);   % Boxcar

        end

        % Windowing neccesary before Fourier transform
        windowing = windowing*ones(1,nx);

        figure(3); set(3,'units','normalized','outerposition',[0 0 1 1]);
        figure(6); set(6,'units','normalized','outerposition',[0 0 1 1]);
        for jj=1:n
            for ii=1:m

                xw = x0(jj) ;   % x window
                zp = z0p(ii);
                zd = z0d(ii);

                sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
                sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

                blockP(ii,jj) = effective_lines( sam1(zp-(nw-1)/2:2:zp+(nw-1)/2,xw:xw+nx-1) );
                blockD(ii,jj) = effective_lines( sam1(zd-(nw-1)/2:2:zd+(nw-1)/2,xw:xw+nx-1) );
                blockT(ii,jj) = effective_lines( sam1(zp-(nw-1)/2:2:zd+(nw-1)/2,xw:xw+nx-1) );

                [tempSp,temp_psnr_Sp] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
                [tempSd,temp_psnr_Sd] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);

                Sp(ii,jj,:) = (tempSp(rang));
                Sd(ii,jj,:) = (tempSd(rang));

                psnr_Sp(ii,jj) = temp_psnr_Sp;
                psnr_Sd(ii,jj) = temp_psnr_Sd;

                if jj==floor(3*n/6)
                    if ii==ceil(m/6)
                        figure(3)
                        set(gca,'FontSize',font);
                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k');
                        title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                        axis([0 25 -70 0]);
                        set(gca,'FontSize',font);

                        figure(6)
                        spmax = max(tempSp((1:NFFT/2+1)));
                        set(gca,'FontSize',font);
                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'k');
                        %title('Spectrum');
                        xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                        axis([0 25 -70 0]);
                        set(gca,'FontSize',font);

                        figure(7); set(7,'units','normalized','outerposition',[0 0 1 1]);
                        s1 = sub_block_p(:,round(end/2));
                        smax = max(abs(s1));
                        set(gca,'FontSize',font);
                        plot(s1/smax,'k');
                        %title('Gated ultrasound signal');
                        xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
                        %axis([0 25 -70 0]);
                        set(gca,'FontSize',font);
                        axis tight;
                        ylim([-1 1])
                    elseif ii==round(m/2)
                        figure(3)
                        hold on;
                        set(gca,'FontSize',font);
                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'r');
                        title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                        axis([0 25 -70 0]);
                        set(gca,'FontSize',font);

                        figure(6)
                        hold on;
                        set(gca,'FontSize',font);
                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'r');
                        %title('Spectrum');
                        xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                        axis([0 25 -70 0]);
                        set(gca,'FontSize',font);

                        figure(8); set(8,'units','normalized','outerposition',[0 0 1 1]);
                        s2 = sub_block_p(:,round(end/2));
                        %hold on;
                        set(gca,'FontSize',font);
                        plot(s2/smax,'r');
                        %title('Gated ultrasound signal');
                        xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
                        %axis([0 25 -70 0]);
                        set(gca,'FontSize',font);
                        axis tight;
                        ylim([-1 1])

                    elseif ii==floor(5*m/6)
                        figure(3)
                        hold on;
                        set(gca,'FontSize',font);
                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'b');   %
                        title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                        axis([0 25 -70 0]);
                        set(gca,'FontSize',font);
                        %pause
                        legend('Top','Half','Bottom');

                        figure(6)
                        hold on;
                        set(gca,'FontSize',font);
                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'b');   %
                        %title('Spectrum');
                        xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity normalized (dB)');
                        axis([0 25 -70 0]);
                        set(gca,'FontSize',font);
                        %pause
                        legend('1 cm above focal dist.','At focal distance','1 cm below focal dist.');

                        figure(9); set(9,'units','normalized','outerposition',[0 0 1 1]);
                        %hold on;
                        s3 = sub_block_p(:,round(end/2));
                        set(gca,'FontSize',font);
                        plot(s3/smax,'b');
                        %title('Gated ultrasound signal');
                        xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
                        %axis([0 25 -70 0]);
                        set(gca,'FontSize',font);
                        %legend('Top','Half','Bottom');
                        axis tight;
                        ylim([-1 1])

                    end

                end

            end

        end
        figure(3)
        hold on;
        plot(band((1:NFFT/2+1))*1e-6,-20*ones(NFFT/2+1,1),'k--');
        plot(band((1:NFFT/2+1))*1e-6,-15*ones(NFFT/2+1,1),'k--');
        plot(band((1:NFFT/2+1))*1e-6,-10*ones(NFFT/2+1,1),'k--');

        %% Diffraction compensation (IDEAL, ANALYTICAL EQUATION)
        load sam1.mat transducer
        for ii=1:m
            zp = z0p(ii);
            zd = z0d(ii);
            zap = z(zp)*1e-2;
            zad = z(zd)*1e-2;
            D_p(ii,1,:) = diffraction_JR(transducer.fn,transducer.fl,band(rang),c0,zap);
            D_d(ii,1,:) = diffraction_JR(transducer.fn,transducer.fl,band(rang),c0,zad);
        end

        for ii=1:n
            D_p(:,ii,:) = D_p(:,1,:);
            D_d(:,ii,:) = D_d(:,1,:);
        end

        diffraction_compensation = log(D_p) - log(D_d);


        save([output_dir,'\spectra_data.mat']);

        %% Au = b
        b = (log(Sp) - log(Sd)) - (diffraction_compensation);
        % watch curves
        temp_ratio = 0; count_ratio = 0;
        for ii = 1
            for jj = round(m/2)
                %for ii = 1:round(m/5):m
                %for jj = 1:round(n/5):n
                figure(5); set(5,'units','normalized','outerposition',[0 0 1 1]); box on;
                set(gca,'FontSize',font);
                temp_ratio = temp_ratio + squeeze(b(ii,jj,:));
                count_ratio = count_ratio + 1;
                plot(f,squeeze(b(ii,jj,:)));
                %title('Spectrum REFERENCE');
                xlabel('\bfFrequency (MHz)'); ylabel('\bfSpectral log ratio (dB/cm)');
                axis tight
                set(gca,'FontSize',font);
                hold on
            end
        end
        hold off;
        figure(55); set(55,'units','normalized','outerposition',[0 0 1 1]); box on;
        plot(f,temp_ratio/count_ratio,'r');

        
        A1 = kron( 4*L*f , speye(m*n) );
        A2 = kron( ones(size(f)) , speye(m*n) );
        A = [A1 A2];
        
        % IMAGESC- WITHOUT REGULARIZATION
        [u,~] = cgs2(A'*A,A'*b(:),1e-6,20);

        % Standard SLD
        % BS: Beta. Attenuation coefficient slopes of blocks.
        % CS: Constants of blocks.
        BS = u(1:end/2); %CS = u(end/2+1:end);
        BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
        BS = reshape(BS,m,n);

        BS_interp = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BS, x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100, 1, [bottom top], '(a)');


        %%         Masks for TV. Noramlly just ones....
        mask1 = ones(m,n);
        mask4 = speye(m*n,m*n);
        mask_diff = speye(2*m*n,2*m*n);

        for kk=1:p
            mask(:,:,kk)=mask1(:,:);
        end

        %%
        disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
            num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
        disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
        disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
        disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
        disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);


        % Inclusion region in SLD
        R1_SLD = square_center(BS_interp,x,z,c1x,c1y,square_L);   % inclusion

        % Background region in SLD
        if rf_class == 102 && homogeneous == 0,
            R2_SLD = square_center(BS_interp,x,z,c2x,c2y,square_L);
        else
            R2_SLD = square_center(BS_interp,x,z,c2x,c2y,square_L);
        end

        % Mean and standard deviation percentage error: MPE, SDPE

        SLD.result(:,:,1) = [nanmean(reshape(R1_SLD(:,:,1),[],1)) nanstd(reshape(R1_SLD(:,:,1),[],1))];
        SLD.result(:,:,2) = [nanmean(reshape(R2_SLD(:,:,1),[],1)) nanstd(reshape(R2_SLD(:,:,1),[],1))];

        SLD.error{1} = (reshape(R1_SLD(:,:,1),[],1) - ground_truth(1))/ground_truth(1);
        SLD.error{2} = [abs(nanmean(SLD.error{1})) nanstd(SLD.error{1}) 100*abs(nanmean(SLD.error{1})) 100*nanstd(SLD.error{1}) ];

        SLD.error{3} = (reshape(R2_SLD(:,:,1),[],1) - ground_truth(2))/ground_truth(2);
        SLD.error{4} = [abs(nanmean(SLD.error{3})) nanstd(SLD.error{3}) 100*abs(nanmean(SLD.error{3})) 100*nanstd(SLD.error{3}) ];

        SLD.bias(1,1) = SLD.error{2}(3);
        SLD.bias(1,2) = SLD.error{4}(3);

        disp(['Standard SLD. '...
            'R1: ',num2str(SLD.result(1,1,1),'%3.2f'),' +/- ',num2str(SLD.result(1,2,1),'%3.2f'),' dB.cm^{-1}.MHz^{-1}. '...
            'R2: ',num2str(SLD.result(1,1,2),'%3.2f'),' +/- ',num2str(SLD.result(1,2,2),'%3.2f'),' dB.cm^{-1}.MHz^{-1}.']);

        SLD.cnr_delta = 0.005;
        SLD.cnr = abs( SLD.result(:,1,1) - SLD.result(:,1,2) )./( SLD.cnr_delta + sqrt( (SLD.result(:,2,1)).^2 + (SLD.result(:,2,2)).^2 ) );
        disp(['CNR_SLD: ',num2str(log10((SLD.cnr)),'%3.2f')])
        save ([output_dir,'\stats_SLD.mat'],'SLD','R1_SLD','R2_SLD');

        %% Multi-frequency. Julien denoising
        % bj = b. For Julien
        jr.mf = 0;
        jr.tol = 1e-3;
        jr.minimask = mask(:,:,1);
        jr.minimask(isnan(jr.minimask))=0;

        if jr.mf == 1

            jr.factor = 1:0.5:3;
            %jr.factor = 1;
            jr.factor = 1;
            for qq = 1:length(jr.factor)

                for kk = 1:p
                    jr.bNoise = b(:,:,kk);   % a single frequency
                    jr.mean(kk,1) = abs(nanmean(jr.bNoise(:)));
                    jr.sd(kk,1) = nanstd(jr.bNoise(:));
                    jr.mu(kk,qq) = jr.factor(qq)*jr.sd(kk,1)/jr.mean(kk,1);
                    jr.bDenoise = IRLS_ANIS_TV(jr.bNoise(:),speye(length(jr.bNoise(:))),jr.mu(kk,qq),m,n,jr.tol,mask,jr.minimask(:));
                    jr.b(:,:,kk) = reshape(jr.bDenoise,m,n);

                end

                for ii = 1:round(m/5):m
                    for jj = 1:round(n/5):n
                        figure(6); set(6,'units','normalized','outerposition',[0 0 1 1]); box on;
                        plot(f,squeeze(jr.b(ii,jj,:)));
                        hold on
                    end
                end
                hold off;


                [jr.u{qq},~] = cgs2(A'*A,A'*jr.b(:),1e-6,20);

                % Standard SLD
                % BS: Beta. Attenuation coefficient slopes of blocks.
                % CS: Constants of blocks.
                BJtemp = jr.u{qq}(1:end/2); %CS = jr.u{qq}(end/2+1:end);
                BJtemp = 8.686*BJtemp;   % [dB.cm^{-1}.MHz^{-1}]
                BJ(:,:,qq) = reshape(BJtemp,m,n);
                BJ(:,:,qq) = BJ(:,:,qq).*mask1;

                BJ_interp(:,:,qq) = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BJ(:,:,qq), x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100+qq, 1, [bottom top], '(a)');




                % Inclusion region in SLD
                R1_JSLD = square_center(BJ_interp,x,z,c1x,c1y,square_L);   % inclusion

                % Background region in SLD
                if rf_class == 102 && homogeneous == 0,
                    R2_JSLD = square_center(BJ_interp,x,z,c2x,c2y,square_L);
                else
                    R2_JSLD = square_center(BJ_interp,x,z,c2x,c2y,square_L);
                end

                % Mean and standard deviation percentage error: MPE, SDPE

                JSLD.result(:,:,1) = [nanmean(reshape(R1_JSLD(:,:,1),[],1)) nanstd(reshape(R1_JSLD(:,:,1),[],1))];
                JSLD.result(:,:,2) = [nanmean(reshape(R2_JSLD(:,:,1),[],1)) nanstd(reshape(R2_JSLD(:,:,1),[],1))];

                JSLD.error{1} = (reshape(R1_JSLD(:,:,1),[],1) - ground_truth(1))/ground_truth(1);
                JSLD.error{2} = [abs(nanmean(JSLD.error{1})) nanstd(JSLD.error{1}) 100*abs(nanmean(JSLD.error{1})) 100*nanstd(JSLD.error{1}) ];

                JSLD.error{3} = (reshape(R2_JSLD(:,:,1),[],1) - ground_truth(2))/ground_truth(2);
                JSLD.error{4} = [abs(nanmean(JSLD.error{3})) nanstd(JSLD.error{3}) 100*abs(nanmean(JSLD.error{3})) 100*nanstd(JSLD.error{3}) ];

                JSLD.bias(1,1) = JSLD.error{2}(3);
                JSLD.bias(1,2) = JSLD.error{4}(3);

                disp(['JR SLD. '...
                    'R1: ',num2str(JSLD.result(1,1,1),'%3.2f'),' +/- ',num2str(JSLD.result(1,2,1),'%3.2f'),' dB.cm^{-1}.MHz^{-1}. '...
                    'R2: ',num2str(JSLD.result(1,1,2),'%3.2f'),' +/- ',num2str(JSLD.result(1,2,2),'%3.2f'),' dB.cm^{-1}.MHz^{-1}.']);

                JSLD.cnr_delta = 0.005;
                JSLD.cnr = abs( JSLD.result(:,1,1) - JSLD.result(:,1,2) )./( JSLD.cnr_delta + sqrt( (JSLD.result(:,2,1)).^2 + (JSLD.result(:,2,2)).^2 ) );

                save ([output_dir,'\stats_JSLD.mat'],'JSLD');


            end

            %b = jr.b;
        end

        %keyboard
        %% Regularization: Au = b
        tol = 1e-3;
        %mu = logspace(2,3,3)';
        %mu = logspace(1,2,3)';
        %mu = [0 1];

        clear mask
        mask1 = ones(m,n);
        mask4 = speye(m*n,m*n);
        mask_diff = speye(2*m*n,2*m*n);

        for kk=1:p
            mask(:,:,kk)=mask1(:,:);
        end

        for mm = 1:length(mu)

            mu1 = mu(mm);
            mu2 = mu1;
            [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

            BR(:,:,mm) = (reshape(Bn*8.686,m,n));
            BR_interp(:,:,mm) = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BR(:,:,mm), x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100+10*mm, 1, [bottom top], ' (b) ');

            %BR_interp(:,:,mm) = attROT(BR_interp(:,:,mm), x, z, angle, rot_lat, rot_axi, bottom, top);


        end

        save([output_dir,'\stats_ACS.mat'])

        %elseif ACS_estimation == 0,

        %   load([output_dir,'\stats_ACS.mat'])

        %end

        %% Selection of regions with differeent attenuation: layers or inclusions
        % Type 1: One circular inclusion
        % R1: Region 1: Inclusion (Standards SLD, Regularized SLD)
        % R2: Region 2: Background (Standards SLD, Regularized SLD)
        % Type_inclusion = input('Type of regions (No Inc. (1), Line(2), Inc_Back(3), Inc_Back_V2 (4)): ');

        switch type_inclusion

            case 1

                for ii = 1:size(BR_interp,3)
                    R1_RSLD(:,:,ii) = square_center(BR_interp(:,:,ii),x,z,c1x,c1y,square_L);
                    R1_RSLDg(:,:,ii) = square_centerg(BR_interp(:,:,ii),x,z,c1x,c1y,square_L);
                    if rf_class == 102 && ground_truth(1) == 1.0,
                        R2_RSLD(:,:,ii) = square_center(BR_interp(:,:,ii),x,z,c2x,c2y,square_L);
                    else
                        R2_RSLD(:,:,ii) = square_center(BR_interp(:,:,ii),x,z,c2x,c2y,square_L);
                    end
                    figure(70+ii); imagesc(R1_RSLDg(:,:,ii)); colormap('jet') ;colorbar;
                    if ii == ceil(size(BR_interp,3)/2)
                        [m2,n2,~ ]= size(R1_RSLD);
                        map_print( R1_RSLD(:,:,ii) , m2 , n2 , 20+ii,  x , z, 2, [bottom top], '(b)')
                        map_print( R2_RSLD(:,:,ii) , m2 , n2 , 30+ii,  x , z, 2, [bottom top], '(b)')
                    end
                end
        end

        for ii = 1:size(BR_interp,3)
            RSLD.result(ii,:,1) = [mu(ii) nanmean(reshape(R1_RSLD(:,:,ii),[],1)) nanstd(reshape(R1_RSLD(:,:,ii),[],1))];
            RSLD.result(ii,:,2) = [mu(ii) nanmean(reshape(R2_RSLD(:,:,ii),[],1)) nanstd(reshape(R2_RSLD(:,:,ii),[],1))];

            disp(['Regularized SLD. ',num2str(RSLD.result(ii,1,1),'%8.3e'),'. '...
                'R1: ',num2str(RSLD.result(ii,2,1),'%3.2f'),' +/- ',num2str(RSLD.result(ii,3,1),'%3.2f'),' dB.cm^{-1}.MHz^{-1}. '...
                'R2: ',num2str(RSLD.result(ii,2,2),'%3.2f'),' +/- ',num2str(RSLD.result(ii,3,2),'%3.2f'),' dB.cm^{-1}.MHz^{-1}.']);
            % inclusion
            RSLD.error{ii,1} = (reshape(R1_RSLD(:,:,ii),[],1) - ground_truth(1))/ground_truth(1);
            RSLD.error{ii,2} = [abs(nanmean(RSLD.error{ii,1})) nanstd(RSLD.error{ii,1}) 100*abs(nanmean(RSLD.error{ii,1})) 100*nanstd(RSLD.error{ii,1}) ];
            RSLD.vrr(ii,1) = 100*( 1 - nanvar(reshape(R1_RSLD(:,:,ii),[],1))/ nanvar(reshape(R1_SLD(:,:,1),[],1)) );
            % background
            RSLD.error{ii,3} = (reshape(R2_RSLD(:,:,ii),[],1) - ground_truth(2))/ground_truth(2);
            RSLD.error{ii,4} = [abs(nanmean(RSLD.error{ii,3})) nanstd(RSLD.error{ii,3}) 100*abs(nanmean(RSLD.error{ii,3})) 100*nanstd(RSLD.error{ii,3}) ];
            RSLD.vrr(ii,2) = 100*( 1 - nanvar(reshape(R2_RSLD(:,:,ii),[],1))/ nanvar(reshape(R2_SLD(:,:,1),[],1)) );

            RSLD.bias(ii,1) = RSLD.error{ii,2}(3);
            RSLD.bias(ii,2) = RSLD.error{ii,4}(3);

        end

        RSLD.cnr_delta = 0.005;
        RSLD.cnr = abs( RSLD.result(:,2,1) - RSLD.result(:,2,2) )./( RSLD.cnr_delta + sqrt( (RSLD.result(:,3,1)).^2 + (RSLD.result(:,3,2)).^2 ) );

        for ffff = 1:length(mu),
            disp(['CNR_RSLD: ',num2str(log10((RSLD.cnr(ffff))),'%3.2f')])
        end

        % table only with mu = 10^2.5 results of MPE, SDPE and CNR
        RSLD.table.inc(ss,:)  = [SLD.result(:,:,1) , SLD.error{2}(1:2) , RSLD.result(2,2:3,1), RSLD.error{2,2}(1:2) ];
        RSLD.table.back(ss,:) = [SLD.result(:,:,2) , SLD.error{4}(1:2) , RSLD.result(2,2:3,2), RSLD.error{2,4}(1:2) ];
        RSLD.table.cnr(ss,:)  = [SLD.cnr RSLD.cnr(2,1)];
        RSLD.table.vrr(ss,1)  = RSLD.vrr(2,1);

        save ([output_dir,'\stats_RSLD.mat'],'RSLD','R1_RSLD','R2_RSLD','BR_interp')

        % Plot regularized results mean and sd for R1 and R2
        figure(10); set(10,'units','normalized','outerposition',[0 0 1 1]); box on;
        %subplot(121)
        set(gca,'FontSize',font,'xscale','log');
        %set(gca,'FontSize',font,'xscale','log');
        errorbar(mu, RSLD.result(:,2,1), RSLD.result(:,3,1),'r-s','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);

        hold on;
        plot(mu, ground_truth(1)*ones(size(mu)),'r--s');
        xlim([RSLD.result(1,1,1) RSLD.result(end,1,1)]);
        axis tight
        ylim([bottom top]);
        %title('(a)');
        xlabel('\bfRegularization parameter \it{\mu}'); ylabel('\bfMean and SD (dB.cm^{-1}.MHz^{-1})');
        set(gca,'FontSize',font);
        %legend('Regularized SLD', 'Through-transmission')
        set(gca,'XScale','log');
        %errorbarlogx();

        figure(11); set(11,'units','normalized','outerposition',[0 0 1 1]); box on;
        %subplot(122)
        set(gca,'FontSize',font,'xscale','log');
        errorbar(mu, RSLD.result(:,2,2), RSLD.result(:,3,2),'b-s','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);
        %errorbarlogx();
        hold on;
        plot(mu, ground_truth(2)*ones(size(mu)),'b--s');
        xlim([RSLD.result(1,1,1) RSLD.result(end,1,1)]);
        axis tight;
        ylim([bottom top]);
        %title('(b)');
        xlabel('\bfRegularization parameter \it{\mu}'); ylabel('\bfMean and SD (dB.cm^{-1}.MHz^{-1})');
        set(gca,'FontSize',font);
        set(gca,'XScale','log');
        %legend('Regularized SLD', 'Through-transmission')

        %%
        figure(12); set(12,'units','normalized','outerposition',[0 0 1 1]); box on;
        %subplot(121)
        set(gca,'FontSize',font,'xscale','log');
        %set(gca,'FontSize',font,'xscale','log');
        errorbar(mu, RSLD.result(:,2,2), RSLD.result(:,3,2),'r-.d','LineWidth',4,'MarkerFaceColor','r','MarkerSize',25);

        %xlim([10^-5 10^11]);

        hold on;
        errorbar(mu, RSLD.result(:,2,1), RSLD.result(:,3,1),'b-s','LineWidth',4,'MarkerFaceColor','b','MarkerSize',25);

        %plot(mu1, ground_truth(1)*ones(size(mu1)),'r--s','LineWidth',4);
        plot(mu, ground_truth(1)*ones(size(mu)),'b--','LineWidth',4);
        plot(mu, ground_truth(2)*ones(size(mu)),'r--','LineWidth',4);
        plot(ones(1,17)*10^1.25, linspace(bottom,top,17), 'k--','LineWidth',4);
        plot(ones(1,17)*10^3.75, linspace(bottom,top,17), 'k--','LineWidth',4);

        %xlim([results_BR(1,1,1) results_BR(end,1,1)]);
        axis tight
        xlim([RSLD.result(1,1,1) RSLD.result(end,1,1)]);
        ylim([bottom top]);
        %title('(a)');
        xlabel('\bfRegularization parameter \mu'); ylabel('\bfMean and SD (dB.cm^{-1}.MHz^{-1})');
        set(gca,'FontSize',font);
        %legend('Regularized SLD', 'Through-transmission')
        set(gca,'XScale','log');
        %errorbarlogx();
        h = legend('Background','Inclusion','Location','NorthEast');
        set(h,'fontsize',font);

        % Plot Bias and CNR
        figure(13); set(13,'units','normalized','outerposition',[0 0 1 1]); box on;
        %subplot(211);
        set(gca,'FontSize',font,'xscale','log');
        %xlim([10^0 10^5]); % semilogx;
        hold on;
        plot(mu, RSLD.cnr,'b-s','LineWidth',4,'MarkerFaceColor','b','MarkerSize',25);
        plot(ones(1,17)*10^1.25, linspace(bottom,5.5,17), 'k--','LineWidth',4);
        plot(ones(1,17)*10^3.75, linspace(bottom,5.5,17), 'k--','LineWidth',4);
        axis tight
        %ylim([0.2 1.4]);
        %ylim([0 2.0]);
        xlabel('\bfRegularization parameter \mu'); ylabel('\bfCNR_{RSLD}');
        set(gca,'FontSize',font);
        % legend('Inclusion','Background','Location','NorthEastOutside');
        ylim([0 4.5])

        figure(14); set(14,'units','normalized','outerposition',[0 0 1 1]); box on;
        %subplot(212);
        set(gca,'FontSize',font,'xscale','log');
        %xlim([10^0 10^5]); % semilogx;
        hold on;
        plot(mu, RSLD.bias(:,2),'r-d','LineWidth',4,'MarkerFaceColor','r','MarkerSize',25);
        plot(mu, RSLD.bias(:,1),'b-s','LineWidth',4,'MarkerFaceColor','b','MarkerSize',25);
        plot(ones(1,17)*10^1.25, linspace(0,50,17), 'k--','LineWidth',4);
        plot(ones(1,17)*10^3.75, linspace(0,50,17), 'k--','LineWidth',4);
        %ylim([0.2 1.4]);
        %ylim([0 2.0]);
        xlabel('\bfRegularization parameter \mu');
        ylabel('\bfMPE (%)');
        set(gca,'FontSize',font);
        h = legend('Background','Inclusion','Location','NorthWest');
        set(h,'fontsize',font);
        plot(mu, 20*ones(size(mu)),'k--','LineWidth',4);
        plot(mu, 10*ones(size(mu)),'k--','LineWidth',4);

        %         close 5 6 7 8 9 10 11 12 13 14 55
        %figure(33); set(gca,'FontSize',20);

        figure(3);
        legend('Top','Half','Bottom'); title('SPECTRA')

        %pause
        %keyboard
        %keyboard
        %         save_all_figures_to_directory(output_dir);
        %
        %         clearvars -except frame_list angle_list db_size font frame ff ss angle
        %         close all;
        
    end    
    
end

% load handel.mat;
% nBits = 16;
% sound(y/1.5,Fs,nBits);

return



% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %


%ACS = 0;
%R1 = 0;
%R2 = 0;

for ff = 1:length(frame_list)
   file_name = [pwd,'\figures\frame',num2str(ff),'\FIG1\stats_RSLD.mat'];
   load(file_name); 
   ACS(:,:,:,ff) =  BR_interp;
   R1(:,:,:,ff) = R1_RSLD;
   R2(:,:,:,ff) = R2_RSLD; 
    
end

load('D:\IUS2017\RSLD_SC\code\figures\frame1\FIG1\stats_ACS.mat','x','z');
ACS = nanmean(ACS,4);
R1 = nanmean(R1,4);
R2 = nanmean(R2,4);

figure; imagesc(x,z,ACS(:,:,1)); caxis([0 1.2]); axis image; colorbar

for kk=1:size(R1,3);

    in(kk,:) = [ nanmean(reshape(R1(:,:,kk),[],1)) nanstd(reshape(R1(:,:,kk),[],1)) ];
    bg(kk,:) = [ nanmean(reshape(R2(:,:,kk),[],1)) nanstd(reshape(R2(:,:,kk),[],1)) ];
end

in
bg

bnew = 0
temp = 0;
for ff = 1:length(frame_list)
   file_name = [pwd,'\figures\frame',num2str(ff),'\FIG1\spectra_data.mat'];
   load(file_name); 
   b = (log(Sp) - log(Sd)) - (diffraction_compensation);
   bnew = bnew + b;
   figure(95); set(95,'units','normalized','outerposition',[0 0 1 1]); box on;
   hold on
   plot(f,squeeze(b(round(m/2),round(n/2),:)));  
end

bnew = bnew/5;
figure(95); set(95,'units','normalized','outerposition',[0 0 1 1]); box on;
plot(f,squeeze(bnew(round(m/2),round(n/2),:)),'r');  

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

[u,~] = cgs2(A'*A,A'*bnew(:),1e-6,20);

% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); %CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);

BS_interp = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BS, x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100, 1, [bottom top], '(a)');

b = bnew; mask = ones(m,n,p); mask1 = mask(:,:,1);

%% Multi-frequency. Julien denoising
            % bj = b. For Julien
            jr.mf = 0;
            jr.tol = 1e-3;
            jr.minimask = mask(:,:,1);
            jr.minimask(isnan(jr.minimask))=0;
            
            if jr.mf == 1

                jr.factor = 1:0.5:3;
                %jr.factor = 1;
                jr.factor = 1;
                for qq = 1:length(jr.factor)
                    
                    for kk = 1:p
                        jr.bNoise = b(:,:,kk);   % a single frequency
                        jr.mean(kk,1) = abs(nanmean(jr.bNoise(:)));
                        jr.sd(kk,1) = nanstd(jr.bNoise(:));
                        jr.mu(kk,qq) = jr.factor(qq)*jr.sd(kk,1)/jr.mean(kk,1);
                        jr.bDenoise = IRLS_ANIS_TV(jr.bNoise(:),speye(length(jr.bNoise(:))),jr.mu(kk,qq),m,n,jr.tol,mask,jr.minimask(:));
                        jr.b(:,:,kk) = reshape(jr.bDenoise,m,n);
                        
                    end
                    
                for ii = 1:round(m/5):m
                    for jj = 1:round(n/5):n
                    figure(6); set(6,'units','normalized','outerposition',[0 0 1 1]); box on;
                    plot(f,squeeze(jr.b(ii,jj,:)));
                    hold on
                    end
                end
                hold off;              
                    
                                        
                    [jr.u{qq},~] = cgs2(A'*A,A'*jr.b(:),1e-6,20);
                    
                    % Standard SLD
                    % BS: Beta. Attenuation coefficient slopes of blocks.
                    % CS: Constants of blocks.
                    BJtemp = jr.u{qq}(1:end/2); %CS = jr.u{qq}(end/2+1:end);
                    BJtemp = 8.686*BJtemp;   % [dB.cm^{-1}.MHz^{-1}]
                    BJ(:,:,qq) = reshape(BJtemp,m,n);
                    BJ(:,:,qq) = BJ(:,:,qq).*mask1;
                    
                    BJ_interp(:,:,qq) = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BJ(:,:,qq), x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100+qq, 1, [bottom top], '(a)');
           
                    
                    
                    
                    % Inclusion region in SLD
            R1_JSLD = square_center(BJ_interp,x,z,c1x,c1y,square_L);   % inclusion
           
            % Background region in SLD
            if rf_class == 102 && homogeneous == 0,
                R2_JSLD = two_square_center(BJ_interp,x,z,c2x,c2y,square_L);  
            else 
                R2_JSLD = square_center(BJ_interp,x,z,c2x,c2y,square_L);    
            end
            
            % Mean and standard deviation percentage error: MPE, SDPE
            
            JSLD.result(:,:,1) = [nanmean(reshape(R1_JSLD(:,:,1),[],1)) nanstd(reshape(R1_JSLD(:,:,1),[],1))];
            JSLD.result(:,:,2) = [nanmean(reshape(R2_JSLD(:,:,1),[],1)) nanstd(reshape(R2_JSLD(:,:,1),[],1))];
            
            JSLD.error{1} = (reshape(R1_JSLD(:,:,1),[],1) - ground_truth(1))/ground_truth(1);
            JSLD.error{2} = [abs(nanmean(JSLD.error{1})) nanstd(JSLD.error{1}) 100*abs(nanmean(JSLD.error{1})) 100*nanstd(JSLD.error{1}) ];
            
            JSLD.error{3} = (reshape(R2_JSLD(:,:,1),[],1) - ground_truth(2))/ground_truth(2);
            JSLD.error{4} = [abs(nanmean(JSLD.error{3})) nanstd(JSLD.error{3}) 100*abs(nanmean(JSLD.error{3})) 100*nanstd(JSLD.error{3}) ];
            
            JSLD.bias(1,1) = JSLD.error{2}(3);
            JSLD.bias(1,2) = JSLD.error{4}(3);  
            
            disp(['JR SLD. '...
                'R1: ',num2str(JSLD.result(1,1,1),'%3.2f'),' +/- ',num2str(JSLD.result(1,2,1),'%3.2f'),' dB.cm^{-1}.MHz^{-1}. '...
                'R2: ',num2str(JSLD.result(1,1,2),'%3.2f'),' +/- ',num2str(JSLD.result(1,2,2),'%3.2f'),' dB.cm^{-1}.MHz^{-1}.']);
            
            JSLD.cnr_delta = 0.005;
            JSLD.cnr = abs( JSLD.result(:,1,1) - JSLD.result(:,1,2) )./( JSLD.cnr_delta + sqrt( (JSLD.result(:,2,1)).^2 + (JSLD.result(:,2,2)).^2 ) );
            
            save ([output_dir,'\stats_JSLD.mat'],'JSLD');
           
                end  
            
                %b = jr.b;
            end
            
            
            
            mu = logspace(2,3,3)';
tol = 1e-3;
%mu = logspace(1,2,3)';
%mu = [0 1];
angle = 0;
clear mask
mask1 = ones(m,n);
mask4 = speye(m*n,m*n);
mask_diff = speye(2*m*n,2*m*n);

for kk=1:p
    mask(:,:,kk)=mask1(:,:);
end

for mm = 1:length(mu)
    
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,bnew(:),mu1,mu2,m,n,tol,mask(:));
    
    BR2(:,:,mm) = (reshape(Bn*8.686,m,n));
    BR2_interp(:,:,mm) = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BR2(:,:,mm), x_ACS, z_ACS, nx_ACS, nz_ACS, m, n, 100+10*mm, 1, [bottom top], ' (b) ');
    
    BR2_interp(:,:,mm) = attROT(BR2_interp(:,:,mm), x, z, angle, rot_lat, rot_axi, bottom, top);
    
    
end



aa = 3;
%sumar


load PRE.mat CNRfinal_SLD CNRfinal_RSLD MPEfinal_SLD MPEfinal_RSLD_incl MPEfinal_RSLD_back VRRfinal_RSLD_incl VRRfinal_RSLD_back line axis_x

%font = 36;

font = 42;
figure(99999); set(99999,'units','normalized','outerposition',[0 0 1 1]); box on; hold on
%subplot(121)
set(gca,'FontSize',font);
plot(axis_x{1,1},line{1,1},'k-s','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);
plot(axis_x{2,1},line{2,1},'r-o','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);
plot(axis_x{3,1},line{3,1},'b-d','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);
plot(axis_x{4,1},line{4,1},'m->','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);
plot(axis_x{5,1},line{5,1},'g-<','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);
%loglog(Lcurve_x_axis,Lcurve_y_axis,'k-s','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);
%title('(a)');
axis tight;
xlabel('\bfLateral distance (cm)');
ylabel('\bfdB.cm^{-1}.MHz^{-1}');
%xlabel('\bf||Au-b||_2');
%ylabel('\bf||Lu||_2');
set(gca,'FontSize',font);
legend('50','40','30','20','10','Location','EastOutside');
%legend('30','20','10','Location','EastOutside');
ylim([0 1.2]);

figure(99999)
L = length(axis_x{5,1});
yy1 = 1.04*ones(L,1);
[~,tt] = find((axis_x{5,1}>-0.5)&(axis_x{5,1}<1.5));
yy1(tt) = 1.04;
plot(axis_x{5,1},yy1,'b--','LineWidth',2);


save JOURNAL_TABLE.mat CNRfinal_SLD CNRfinal_RSLD MPEfinal_SLD MPEfinal_RSLD_incl MPEfinal_RSLD_back VRRfinal_RSLD_incl VRRfinal_RSLD_back line axis_x

end

%% Diffraction (fn,fl,band,c0,za)

function D = diffraction_JR(fn,fl,band,c0,za)

% Calcula el factor de compensacion por difraccion en la profundidad
% de analisis za a la frecuencia fa (X.Chen)
%
% Entradas:     fn: #focal del transductor
%               fl: Longitud focal del tx (m)
%               band: Frecuencia de analisis (Hz)
%               c: Velocidad del sonido en el medio (m/s)
%               za: Profundidad de analisis (m)
%
% Salidas:      D: Factor de compensacion por difraccion

a  = (fl/(fn*2));                        % Radio del transductor.
Gp = ( (pi*band/c0)*(a^2) ) / fl;           % Preasure gain factor.

for i = 1 : length(Gp)
    
    if 1/(1+pi/Gp(i))<=za/fl && za/fl<= 1/(1-pi/Gp(i))
        D(i) = ((pi*a^2)/(za^2))* 0.46*exp(-(0.46/pi)*Gp(i).^2*((fl/za)-1).^2);
    else
        %D(i) = ((pi*a^2)/(za.^2))* 1.07*(Gp(i)*((fl/za)-1)).^-2;
        D(i) = ((pi*a^2)/(za.^2))* 1.00*(Gp(i)*((fl/za)-1)).^-2;
    end
    % To PLOT
    %D(i)= D(i)/((pi*a^2)/(za^2));
    
end

end


%% Map print
% v: data vector
% m: number of rows
% n: number of columns
% i: figure
% x: lateral distance
% z: axial distance
% k: type of colormap
% c: caxis range
% t: title

%function [] = map_print2( Bmode, v , m , n , i,  x , z, k, c, t)
function [swe_Vi] = map_print2( Im_db_ori, x_ori, z_ori, Im_db, x, z, BS, x_ate, z_ate, nx_ate, nz_ate, m, n, i, k, c, t)

%sl = 60;
sl = 42;

Bmode = Im_db;

SWS_PD = BS;
%Properties.Width_B = x_ori*1e3;
%Properties.Depth_B = z_ori*1e3;
Properties.Width_B = x*10;
Properties.Depth_B = z*10;

Properties.Width_S = x*10;
Properties.Depth_S = z*10;

rd = c;
% BMode with SWS map superposed
% close all
BW=ones(size(Bmode));
[X1,Y1] = meshgrid(1:size(SWS_PD,2),1:size(SWS_PD,1));
[X2,Y2] = meshgrid(linspace(1,size(SWS_PD,2),size(Bmode,2)),linspace(1,size(SWS_PD,1),size(Bmode,1)));

%[XX1,YY1] = meshgrid(x_ate,z_ate);
%[XX2,YY2] = meshgrid(x,z);

swe_Vi = interp2(X1,Y1,SWS_PD,X2,Y2);


%swe_Vi = interp2(XX1,YY1,SWS_PD,XX2,YY2,'linear');
%figure; imagesc(swe_Vi); colorbar;
%figure; imagesc(SWS_PD); colorbar; caxis([0,1])

% a=Properties.Depth_S;
%a=Properties.Depth_S-Properties.Depth_S(1);
a = Properties.Depth_S;
%figure(500)
figure(i);
set(i,'units','normalized','outerposition',[0 0 1 1]);

transparency=0.55;
alpha_mat = (1-transparency)*ones(size(Bmode));

Bmode = Im_db_ori;
Properties.Width_B = x_ori*1e3;
Properties.Depth_B = z_ori*1e3;

caxis_bmode = [-60 0];
breast = 10;
if breast == 1,
    im_breast = imread('D:\JOURNAL2016\code\data\Breast\sam4\18-20-39 - journal.png');
    %figure; imshow(im_breast);
    im_breast = rgb2gray(im_breast);
    Bmode = imresize(double(im_breast),size(Bmode));
    % figure; imagesc(Bmode); colormap gray;
    %keyboard
    % x_ori = linspace(x_ori(1),x_ori(end),size(im_breast,2));
    % z_ori = linspace(z_ori(1),z_ori(end),size(im_breast,1));
    % Properties.Width_B = x_ori*1e3;
    % Properties.Depth_B = z_ori*1e3;
    %Bmode = double(im_breast);
    caxis_bmode = [0 255];
end


%subimage(1e-1*Properties.Width_B, 1e-1*Properties.Depth_B, 64*mat2gray(Bmode,[-60 0]), gray(64));

subimage(1e-1*Properties.Width_B, 1e-1*Properties.Depth_B, 64*mat2gray(Bmode,caxis_bmode), gray);
%subimage(1e-1*Properties.Width_B, 1e-1*Properties.Depth_B, 64*mat2gray(Bmode), gray);

% imagesc(1e3*Properties.Width_B, 1e3*Properties.Depth_B,Bmode)
hold on
%h=subimage(1e-1*Properties.Width_S, 1e-1*a, 64*mat2gray(BW.*swe_Vi, [rd]), jet(64));
h=subimage(1e-1*Properties.Width_S, 1e-1*a, 64*mat2gray(BW.*swe_Vi, [rd]), jet);

% h=subimage(1e3*Properties.Width_S, 1e3*Properties.Depth_S, swe_Vi);
%set(gca,'FontSize',sl);


%set( h, 'AlphaData', alpha_mat) ; % .5 transparency
set(gca,'FontSize',sl);

%xlabel('Width [mm]','fontsize',16);ylabel('Depth [mm]','fontsize',sl)
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
set(gca,'FontSize',sl);
%title(['Vibration Frequency ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)

%h2 = colorbar;
%set(get(h2,'YLabel'),'String','SWS [m/s]','FontSize',sl);
h2 = colorbar;
ylabel(h2,'dB.cm^{-1}.MHz^{-1}','FontSize', sl);

n_tick = rd(end);  % number of tick marks desired on color axis
%n_tick = 3;
% get the current colorbar axis limits
cbar_lim = get(h2, 'YLim');
% label colorbar correctly
set(h2, 'YTick',cbar_lim(1):range(cbar_lim)/n_tick:cbar_lim(2));
set(h2,'YTickLabel', rd(1):rd(end),'FontSize',sl);
set(gca,'FontSize',sl);


hold off

end


%% Map print
% v: data vector
% m: number of rows
% n: number of columns
% i: figure
% x: lateral distance
% z: axial distance
% k: type of colormap
% c: caxis range
% t: title

function [] = map_print( v , m , n , i,  x , z, k, c, t)

sl = 36;   % Size letter for figures
sl = 60;
figure(i); set(i,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'FontSize',sl);
imagesc(x,z,reshape(v,m,n)); axis image; caxis(c);
hb1=colorbar; ylabel(hb1,'dB.cm^{-1}.MHz^{-1}','FontSize', sl)
%title(t);
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
set(gca,'FontSize',sl);

switch k
    case 1
        colormap hot;
    case 2
        colormap jet;
    case 3
        colormap pink;
        
    otherwise
        disp('Error 4357!');
        
end

end

%% Region type: One circular inclusion

% Find center
% input: coordinates of three points
% output: center

function [c_row,c_col,radius]=find_center(x,y)
% Line 1 - Perpendicular to the line (x1,y1),(x2,y2)
% a1*x+b1*y=c1
mp1=[(x(1)+x(2))/2 (y(1)+y(2))/2];   % mid-point
a1=x(2)-x(1);
b1=y(2)-y(1);
c1=a1*mp1(1)+b1*mp1(2);

% Line 2 - Perpendicular to the line (x2,y2),(x3,y3)
% a2*x+b2*y=c2
mp2=[(x(2)+x(3))/2 (y(2)+y(3))/2];   % mid-point
a2=x(3)-x(2);
b2=y(3)-y(2);
c2=a2*mp2(1)+b2*mp2(2);

% Ax=b
A= [a1 b1; a2 b2]; b= [c1; c2];
[center]=cgs(A,b,1e-6,20);
c_col=center(1);
c_row=center(2);

radiusX(1,1) = sqrt((c_col - x(1))^2 + (c_row-y(1))^2);
radiusX(2,1) = sqrt((c_col - x(1))^2 + (c_row-y(1))^2);
radiusX(3,1) = sqrt((c_col - x(1))^2 + (c_row-y(1))^2);
radius = mean(radiusX);

end

%% square inclusion
function img_inside = square_center(img,x,z,c1x,c1y,L)
%L = 0.6;   % [cm]
u = size(img);
mx = c1x;
my = c1y;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2),
    for ii=1:u(1),
        if abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) - mx) > L/2 || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > L/2
            %if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
            img(ii,jj) = NaN;
        end
    end
end
img_inside = img;
end

function img_inside = square_centerg(img,x,z,c1x,c1y,L)
%L = 0.6;   % [cm]
u = size(img);
mx = c1x;
my = c1y;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2),
    for ii=1:u(1),
        if abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) - mx) > L/2 || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > L/2
            %if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
            %img(ii,jj) = NaN;
        end
    end
end
img_inside = img;
end


%% square inclusion
function img_inside = two_square_center(img,x,z,c1x,c1y,L)
%L = 0.6;   % [cm]
u = size(img);
mx = c1x;
my = c1y;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
img1 = img; img2 = img;
for jj=1:u(2),
    for ii=1:u(1),
        if  abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) - mx) > L/4 || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > L/2 ,
            img1(ii,jj) = 0;
        end
        
        if  abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) + mx) > L/4  || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > L/2 ,
            img2(ii,jj) = 0;
        end
    end
end
img = img1 + img2;
img(img==0) = NaN;
img_inside = img;
end


%% square inclusion
function img_inside = rectangular_width(img,x,z,c1x,c1y,L)
%L = 0.6;   % [cm]
u = size(img);
mx = c1x;
my = c1y;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2),
    for ii=1:u(1),
        if abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) - mx) > L || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > L/2
            %if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
            img(ii,jj) = NaN;
        end
    end
end
img_inside = img;
end


%% square inclusion
function img_inside = rectangular_high(img,x,z,c1x,c1y,L)
%L = 0.6;   % [cm]
u = size(img);
mx = c1x;
my = c1y;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2),
    for ii=1:u(1),
        if abs(x(1) + (jj-1)/(u(2)-1)*(x(end)-x(1)) - mx) > L/2 || abs(z(1) + (ii-1)/(u(1)-1)*(z(end)-z(1)) - my) > L
            %if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
            img(ii,jj) = NaN;
        end
    end
end
img_inside = img;
end




%% circle region
% function img_inside = circle_center(img,x,z,c1x,c1y,radius)
% u = size(img);
% mx = c1x;
% my = c1y;
% %fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
% %fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
% for jj=1:u(2),
%     for ii=1:u(1),
%         if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
%             img(ii,jj) = NaN;
%         end
%     end
% end
% img_inside = img;
% end


% circle inclusion
function img_inc = circle_inclusion(img,x,z,c1x,c1y)
[c_row,c_col,radius] = find_center(c1x,c1y);
u = size(img);
mx = c_col;
my = c_row;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2),
    for ii=1:u(1),
        if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 > radius^2,
            img(ii,jj) = NaN;
        end
    end
end
img_inc = img;
end


% circle background
% function img_back = circle_background(img,x,z,c1x,c1y,c2x,c2y)
function img_back = circle_background(img,x,z,c1x,c1y)

[c_row,c_col,radius] = find_center(c1x,c1y);
u = size(img);
mx = c_col;
my = c_row;
%fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
%fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
for jj=1:u(2),
    for ii=1:u(1),
        if (jj/u(2)*(x(end)-x(1))+ x(1)-mx)^2 + (ii/u(1)*(z(end)-z(1))+ z(1)-my)^2 <= radius^2,
            img(ii,jj) = NaN;
        end
    end
end
%
% [c2_row,c2_col,radius2] = find_center(c2x,c2y);
% v = size(img);
% mx2 = c2_col;
% my2 = c2_row;
% %fx = 1/u(2)*(x(end)-x(1))*1e3+ x(1)*1e3;
% %fy = 1/u(1)*(z(end)-z(1))*1e3+ z(1)*1e3;
% for jj=1:v(2),
%     for ii=1:v(1),
%         if (jj/v(2)*(x(end)-x(1))+ x(1)-mx2)^2 + (ii/v(1)*(z(end)-z(1))+ z(1)-my2)^2 > radius2^2,
%             img(ii,jj) = NaN;
%         end
%     end
% end

img_back = img;

end

%% Differences for attenuations B and constants K
% m: number of rows
% n: number of columns (it must include the columns for B and C)
% Wx: weight-horizontal
% Wy: weight-vertical
function [B, Bt, BtB] = DiffOper(m,n,Wx,Wy)
D1 = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1); % col
D1(:,1) = [];
D1(1,1) = 0;
D2 = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1); % row
D2(:,1) = [];
D2(1,1) = 0;
%B = [ kron(speye(m),D) ; kron(D,speye(m)) ];
B = [ Wx*kron(D1,speye(m)) ; Wy*kron(speye(n),D2)];

%B = [ kron(speye(m),D) ]; %Vertical
% B = [ kron(D,speye(m)) ]; % Horizontal
%B = [B B];
%B = [B sparse(size(B,1),size(B,2))];
Bt = B';
BtB = Bt*B;
end

% rosin 1 normal, rosin 2 para L curve
function [threshold,ind2]=rosin(X,Y)

%clear all; close all; clc; home;

% rosin 1
%[Y1, ind] = max(Y);
%Y2 = Y(end);
%X1 = X(ind);
%X2 = X(end);

% rosin 2
ind = 1;
Y1 = Y(ind);
Y2 = Y(end);
X1 = X(ind);
X2 = X(end);

% Line Equation
% y = a*x +b
% y -a*x - b = 0
a = (Y2-Y1)/(X2-X1);
b = Y1 - a*X1;
% Distance from (x0,y0)
% abs(y0 - a*x0 -b)/sqrt(1+a^2)

Y=Y(ind:end);
X=X(ind:end);

distance_rosin = abs(Y - a*X - b)/sqrt(1+a^2);
[~, ind2] = max(distance_rosin);
threshold = X(ind2);

end

function save_all_figures_to_directory(dir_name)

figlist=findobj('type','figure');

for i=1:numel(figlist)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.fig']));
    figure(figlist(i))
    set(gcf,'PaperPositionMode','auto')
    %pause(2)
    saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i).Number) '.png']));
    %pause(2)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.eps']));
    
end

end

function [B,C] = AlterOpti_ADMM(A1,A2,b,mu1,mu2,m,n,tol,mask)

p = length(mask)/(m*n);
minimask = reshape(mask,[m n p]);
minimask = minimask(:,:,1);
minimask = minimask(:);
b = mask.*b;
b(isnan(b)) = 0;
A = [A1 A2];
[u,~] = cgs2(A'*A,A'*b,1e-6,20);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);
%figure(109); imagesc(8.686*reshape(B,m,n)); colormap pink; caxis([0 1.2])
B = B(:);
C = C(:);
D = 0;
v = 0;

F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);

ite  = 0;
error = 1;

while abs(error) > tol && ite < 20;
    ite = ite + 1;
    
    rho = 1;
    % First part of ADMM algorithm: B
   % B = IRLS_ANIS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);    
    B = IRLS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);    
    %F(2*ite,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);
    %error = F(2*ite,1) - F(2*ite-1,1);
    
    % Second part of ADMM algorithm: C
   % C = IRLS_ANIS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
    C = IRLS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
    %F(2*ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);
    %error = F(2*ite+1,1) - F(2*ite,1);
    
    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;
    
    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;    
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);

end

end


% TV Andres Leonel Coila
function [TV] = TVcalc(B,M,N,mask)
% Returns the Total Variation (isotropic) 
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%  
% Outputs:
%       TV              number

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = Dh.^2 + Dv.^2;
P = sqrt(P);
TV = norm(P(:).*mask,1);
end

% TV Andres Leonel Coila - ANISOTROPIC
function [TV] = TVcalc2(B,M,N,mask)
% Returns the Total Variation (anisotropic) 
% Inputs: 
%       B               vector containing an image
%       M,N             image size
%       mask            binary mask of analyzed region
%  
% Outputs:
%       TV              number
mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = abs(Dh) + abs(Dv); % I don't get this
%P = sqrt(P);
TV = norm(P(:).*mask,1);

end

% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,~,minimask)
[u,~] = cgs2(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol
    
    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];
    
    %Dx*X(:) - Dh(:);
    %Dy*X(:) - Dv(:);
    
    P = Dh.^2 + Dv.^2;
    eps = 0.01;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    W = kron(speye(2),omega);
    
    AtA = A'*A;
    %mu=5000;
    %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
    [u] = cgs2( AtA + mu*D'*W*D, A'*b, 1e-6 , 20, u );
    
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));
    
end
end

% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_ANIS_TV(b,A,mu,M,N,tol,~,minimask)
% Iteratively reweighted least squares (anisotropic)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc2(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol && ite_irls < 20
    
    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];
    
    %Dx*X(:) - Dh(:);
    %Dy*X(:) - Dv(:);
    
    P = Dh.^2 + Dv.^2;
    eps = 0.1;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    %W = kron(speye(2),omega);
    
    Px = abs(Dh + eps);
    Px = 1./Px;
    Px = Px(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Px,0,omega);
    Wx = kron(speye(1),omega);
    
    Py = abs(Dv + eps);
    Py = 1./Py;
    Py = Py(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Py,0,omega);
    Wy = kron(speye(1),omega);
    
    AtA = A'*A;
    %mu=5000;
    %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
    [u] = cgs2( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , A'*b, 1e-6 , 20, u );
    
    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc2(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));
    
end

%figure(909); plot(1:length(G),G);

end


% A*x = b
function [x,ite_cgs] = cgs2(A,b,tol,maxit,varargin)
% Solves the system of linear equations A*x = b for x using the 
% CONJUGATE GRADIENTS SQUARED METHOD
% Inputs: 
%       A,b         inputs
%       tol         tolerance for error norm
%       maxit       maximum iterations
%       varargin    initial guess for x
%  
% Outputs:
%       x           column vector containing the answer
%       ite_cgs     number of iterations

if length(varargin) == 1
    x = varargin{1};
else
    x= zeros(size(A,2),1);
end
r = A*x - b;
p = -r;
ite_cgs = 0;

while norm(r,2) > tol && ite_cgs < maxit
    alpha = (r'*r)/(p'*(A*p));
    x = x + alpha*p;
    rn = r + alpha*(A*p);
    beta = (rn'*rn)/(r'*r);
    p = -rn + beta*p;
    r = rn;
    ite_cgs = ite_cgs + 1;
end

end

%%
function Neff = effective_lines(data_block)

[~, N] =size(data_block);
RHO = corrcoef(data_block);

rho = diagSumCalc(RHO,1); % sum of the main diagonal
rho = rho(1:N-1)./(1:N-1);

val = (rho.^2)*(1:N-1)';
Neff = N./( 1 + 2/N*( val ) );

end


%%
function [diagSum] = diagSumCalc(squareMatrix, LLUR0_ULLR1)
%
% Input: squareMatrix: A square matrix.
%        LLUR0_ULLR1:  LowerLeft to UpperRight addition = 0
%                      UpperLeft to LowerRight addition = 1
%
% Output: diagSum: A vector of the sum of the diagnols of the matrix.
%
% Example:
%
% >> squareMatrix = [1 2 3;
%                    4 5 6;
%                    7 8 9];
%
% >> diagSum = diagSumCalc(squareMatrix, 0);
%
% diagSum =
%
%       1 6 15 14 9
%
% >> diagSum = diagSumCalc(squareMatrix, 1);
%
% diagSum =
%
%       7 12 15 8 3
%
% Written by M. Phillips
% Oct. 16th, 2013
% MIT Open Source Copywrite
% Contact mphillips@hmc.edu fmi.
%

if (nargin < 2)
    disp('Error on input. Needs two inputs.');
    return;
end

if (LLUR0_ULLR1 ~= 0 && LLUR0_ULLR1~= 1)
    disp('Error on input. Only accepts 0 or 1 as input for second condition.');
    return;
end

[M, N] = size(squareMatrix);

if (M ~= N)
    disp('Error on input. Only accepts a square matrix as input.');
    return;
end

diagSum = zeros(1, M+N-1);

if LLUR0_ULLR1 == 1
    squareMatrix = rot90(squareMatrix, -1);
end

for i = 1:length(diagSum)
    if i <= M
        countUp = 1;
        countDown = i;
        while countDown ~= 0
            diagSum(i) = squareMatrix(countUp, countDown) + diagSum(i);
            countUp = countUp+1;
            countDown = countDown-1;
        end
    end
    if i > M
        countUp = i-M+1;
        countDown = M;
        while countUp ~= M+1
            diagSum(i) = squareMatrix(countUp, countDown) + diagSum(i);
            countUp = countUp+1;
            countDown = countDown-1;
        end
    end
end

end


function imOut = attROT(im,x,z,iAngle,x_center,y_center,bottom, top)


m = z(end) - z(1);
n = x(end) - x(1);

M = round(m*1e3);
N = round(n*1e3);

m0 = round( (y_center - z(1))*1e3 );
n0 = round( (x_center - x(1))*1e3 );

BR_fr = im;
%BR_fr = BS;

BR_fr = imresize(BR_fr,[M N]);

figure(16); imagesc(x,z,BR_fr); colorbar; axis image; caxis([bottom top])

BR_output = rotateAround(BR_fr,m0,n0,iAngle);
BR_output(BR_output==0) = NaN;

BR_output = imresize(BR_output,[ length(z) length(x) ]);

figure(17); imagesc(x,z,BR_output); colorbar; axis image; caxis([bottom top])

imOut = BR_output;

end

%
% function matrix_rotation =
%
% x = 1:10;
% y = 1:10;
% % choose a point which will be the center of rotation
% x_center = x(3);
% y_center = y(3);
% % define a 60 degree counter-clockwise rotation matrix
% theta = -12; % pi/3 radians = 60 degrees
% R = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
% %Define affine transformation for translation
% a = [1 0 x_center;0 1 y_center; 0 0 1];
% c = [1 0 -x_center;0 1 -y_center; 0 0 1];
% M = a*R*c;
% %M =  c * R * a
% for i=1:10
% rot(:,i) = M*[x(i) y(i) 1]';
% end
% % pick out the vectors of rotated x- and y-data
% x_rot = rot(1,:);
% y_rot = rot(2,:);
% % make a plot
% plot(x,y,'k-',x_rot, y_rot, 'r-', x_center, y_center, 'bo');
% axis equal
%
