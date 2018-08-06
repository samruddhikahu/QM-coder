clc;
clear all;
close all;
%Tstrt = tic;
tic;
IPFolder = fullfile('D:\MATLAB\JPEG with ArithmeticQM','Input Images');
OPFolder = fullfile('D:\MATLAB\JPEG with ArithmeticQM','Reconstructed Images');
if ~exist(OPFolder, 'dir')
    mkdir(OPFolder)
end
OPXFileName = fullfile(OPFolder,'Data1.xls');
toc
for k10 = 1%:25               %%% The input image number which is to be encoded.
    tic;
    clear A CY C1 CCb CCr B I Ig Cres Crgb CYrec CCbrec CCbrs CCrrec CCrrs;
    clear Iy Seq_dcy Bout_dcy Seq_acy Bout_acy Seq_dcyr Seq_acyr dc_y Ir_y Irev_y;
    clear Ib Seq_dcb Bout_dcb Seq_acb Bout_acb Seq_dcbr Seq_acbr dc_b Ir_b Irev_b;
    clear Ir Seq_dcr Bout_dcr Seq_acr Bout_acr Seq_dcrr Seq_acrr dc_r Ir_r Irev_r;
    
    clc;
    close all;
    IPBaseFileName = sprintf('Image%d.bmp',k10);
    IPFullFileName = fullfile(IPFolder,IPBaseFileName);
    A = imread(IPFullFileName);
    imshow(A);
    title('Original Image');
    
    C1 = rgb2ycbcr(A);
    CY = C1(:,:,1);
    CCb = C1(:,:,2);
    CCr = C1(:,:,3);
    
    [H, W] = size(C1(:,:,1));
    sz = [H W];
    
    load('Nmat.mat');
    load('Cnmat.mat');
    %load('Prob_Est_Table1.mat');
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENCODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T1 = tic;
    %%%%%%%%%%%%%%%%%% Coding the Y Component %%%%%%%%%%%%%%%%%%%%%
    tic;
    s = 1;
    %I = zeros(m*n,64);
    for i = 1:8:H-7
        for j = 1:8:W-7
            b = CY(i:i+7,j:j+7);
            b = double(b);
            % Level shifting by -128
            bshifted = b-128;
            % Calculating 2D-DCT:
            bdct = dct2(bshifted);
            %Edct(s) = sum(sum(bdct.^2))/64;
            % Quantizing using the normalization matrix defined above:
            bq = round(bdct./Nmat);
            %Eq(s) = sum(sum(bq.^2))/64;
            Iy(s,:) = zigzag(bq);
            s = s + 1;
        end
    end
    toc
    
    % DC Coefficient Sequence Generation:
    tic;
    Seq_dcy = dc_seq(Iy(:,1)');
    toc
    
    % Encoding using Arithmetic QM coding:
    tic;
    Bout_dcy = QMcoder(Seq_dcy);
    toc
    
    %%%% Encoding the AC Coefficients:
    %Generating binary sequence from the AC coefficients:
    tic;
    Seq_acy = ac_seq(Iy);
    toc
    
    %%% Encoding using Arithmetic QM Coding:
    tic;
    Bout_acy = QMcoder(Seq_acy);
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Coding Cb Component %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Subsampling ratio used 4:2:0.
    clear CCbs;
    CCbs = zeros(H/2,W/2);
    k = 1;
    for i = 1:2:H
        l = 1;
        for j = 1:2:W
            G = CCb(i:i+1,j:j+1);
            Avg = (sum(sum(G)))/4;
            CCbs(k,l) = Avg;
            l = l + 1;
        end
        k = k + 1;
    end
    
    tic;
    s = 1;
    %I = zeros(m*n,64);
    for i = 1:8:(H/2)-7
        for j = 1:8:(W/2)-7
            b = CCbs(i:i+7,j:j+7);
            b = double(b);
            % Level shifting by -128
            bshifted = b-128;
            % Calculating 2D-DCT:
            bdct = dct2(bshifted);
            %Edct(s) = sum(sum(bdct.^2))/64;
            % Quantizing using the normalization matrix defined above:
            bq = round(bdct./Cnmat);
            %Eq(s) = sum(sum(bq.^2))/64;
            Ib(s,:) = zigzag(bq);
            s = s + 1;
        end
    end
    toc
    
    % DC Coefficient Sequence Generation:
    tic;
    Seq_dcb = dc_seq(Ib(:,1)');
    toc
    
    % Encoding using Arithmetic QM coding:
    tic;
    Bout_dcb = QMcoder(Seq_dcb);
    toc
    
    %%%% Encoding the AC Coefficients:
    %Generating binary sequence from the AC coefficients:
    tic;
    Seq_acb = ac_seq(Ib);
    toc
    
    %%% Encoding using Arithmetic QM Coding:
    tic;
    Bout_acb = QMcoder(Seq_acb);
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Coding Cr Component %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Subsampling ratio used 4:2:0.
    clear CCrs;
    CCrs = zeros(H/2,W/2);
    k = 1;
    for i = 1:2:H
        l = 1;
        for j = 1:2:W
            G = CCr(i:i+1,j:j+1);
            Avg = (sum(sum(G)))/4;
            CCrs(k,l) = Avg;
            l = l + 1;
        end
        k = k + 1;
    end
    
    tic;
    s = 1;
    %I = zeros(m*n,64);
    for i = 1:8:(H/2)-7
        for j = 1:8:(W/2)-7
            b = CCrs(i:i+7,j:j+7);
            b = double(b);
            % Level shifting by -128
            bshifted = b-128;
            % Calculating 2D-DCT:
            bdct = dct2(bshifted);
            %Edct(s) = sum(sum(bdct.^2))/64;
            % Quantizing using the normalization matrix defined above:
            bq = round(bdct./Cnmat);
            %Eq(s) = sum(sum(bq.^2))/64;
            Ir(s,:) = zigzag(bq);
            s = s + 1;
        end
    end
    toc
    
    % DC Coefficient Sequence Generation:
    tic;
    Seq_dcr = dc_seq(Ir(:,1)');
    toc
    
    % Encoding using Arithmetic QM coding:
    tic;
    Bout_dcr = QMcoder(Seq_dcr);
    toc
    
    %%%% Encoding the AC Coefficients:
    %Generating binary sequence from the AC coefficients:
    tic;
    Seq_acr = ac_seq(Ir);
    toc
    
    %%% Encoding using Arithmetic QM Coding:
    tic;
    Bout_acr = QMcoder(Seq_acr);
    toc
    Tcomp(k10) = toc(T1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T2 = tic;
    %%%%%%%%%%%%%%%%%%%%%%% Decoding Y component %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Decoding DC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_dcyr = QMdecoder(Bout_dcy, length(Seq_dcy));
    toc
    
    % Reconstructing DC coefficients using Seq_dcr:
    tic;
    dc_y = dc_seq_re(Seq_dcyr, H, W);
    toc
    
    %%% Decoding AC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_acyr = QMdecoder(Bout_acy, length(Seq_acy));
    toc
    
    % Decoding AC coefficients using Seq_ACr:
    tic;
    Ir_y = ac_seq_re(Seq_acyr, H, W);
    toc
    
    tic;
    Irev_y = [dc_y' Ir_y];
    k = 1;
    for p = 1:8:H-7
        for q = 1:8:W-7
            Z = zigzag_re(Irev_y(k,:));
            % Reverse Quantizing using the same normalization matrix:
            brq = Z.*Nmat;
            % Calculating 2D-IDCT:
            bidct = idct2(brq);
            % Level shifting by +128
            b = bidct+128;
            CYrec(p:p+7,q:q+7) = b;
            k = k + 1;
        end
    end
    %figure, imshow(uint8(CYrec));
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%% Decoding Cb Component %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Decoding DC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_dcbr = QMdecoder(Bout_dcb, length(Seq_dcb));
    toc
    
    % Reconstructing DC coefficients using Seq_dcr:
    tic;
    dc_b = dc_seq_re(Seq_dcbr, H/2, W/2);
    toc
    
    %%% Decoding AC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_acbr = QMdecoder(Bout_acb, length(Seq_acb));
    toc
    
    % Decoding AC coefficients using Seq_ACr:
    tic;
    Ir_b = ac_seq_re(Seq_acbr, H/2, W/2);
    toc
    
    tic;
    Irev_b = [dc_b' Ir_b];
    k = 1;
    for p = 1:8:(H/2)-7
        for q = 1:8:(W/2)-7
            Z = zigzag_re(Irev_b(k,:));
            % Reverse Quantizing using the same normalization matrix:
            brq = Z.*Cnmat;
            % Calculating 2D-IDCT:
            bidct = idct2(brq);
            % Level shifting by +128
            b = bidct+128;
            CCbrec(p:p+7,q:q+7) = b;
            k = k + 1;
        end
    end
    toc
    
    % Resizing the Cb component:-
    % Subsampling ratio used: 4:2:0.
    tic;
    %CCbs = zeros(m/2,n/2);
    clear CCbrs;
    k = 1;
    for i = 1:2:H
        l = 1;
        for j = 1:2:W
            CCbrs(i:i+1,j:j+1) = CCbrec(k,l);
            l = l + 1;
        end
        k = k + 1;
    end
    
    %figure, imshow(uint8(CCbrs));
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%% Decoding Cr Component %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Decoding DC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_dcrr = QMdecoder(Bout_dcr, length(Seq_dcr));
    toc
    
    % Reconstructing DC coefficients using Seq_dcr:
    tic;
    dc_r = dc_seq_re(Seq_dcrr, H/2, W/2);
    toc
    
    %%% Decoding AC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_acrr = QMdecoder(Bout_acr, length(Seq_acr));
    toc
    
    % Decoding AC coefficients using Seq_ACr:
    tic;
    Ir_r = ac_seq_re(Seq_acrr, H/2, W/2);
    toc
    
    tic;
    Irev_r = [dc_r' Ir_r];
    k = 1;
    for p = 1:8:(H/2)-7
        for q = 1:8:(W/2)-7
            Z = zigzag_re(Irev_r(k,:));
            % Reverse Quantizing using the same normalization matrix:
            brq = Z.*Cnmat;
            % Calculating 2D-IDCT:
            bidct = idct2(brq);
            % Level shifting by +128
            b = bidct+128;
            CCrrec(p:p+7,q:q+7) = b;
            k = k + 1;
        end
    end
    toc
    
    % Resizing the Cb component:-
    % Subsampling ratio used: 4:2:0.
    tic;
    %CCbs = zeros(m/2,n/2);
    clear CCrrs;
    k = 1;
    for i = 1:2:H
        l = 1;
        for j = 1:2:W
            CCrrs(i:i+1,j:j+1) = CCrrec(k,l);
            l = l + 1;
        end
        k = k + 1;
    end
    
    %figure, imshow(uint8(CCrrs));
    toc
    
    Cres(:,:,1) = CYrec;
    Cres(:,:,2) = CCbrs;
    Cres(:,:,3) = CCrrs;
    
    Cres = Cres./255;
    figure, imshow(uint8(Cres));
    
    %Crgb = ycbcr2rgb(Cres);
    Crgb = colorspace('YPbPr->RGB',Cres);
    Crgb = uint8(Crgb.*255);
    figure, imshow(Crgb);
    title('Reconstructed Image')
    OPBaseFileName = sprintf('OPImage%d.bmp', k10);
    OPFullFileName = fullfile(OPFolder,OPBaseFileName);
    imwrite(Crgb,OPFullFileName,'bmp');
    Tdecomp(k10) = toc(T2)
    
    %%%%%%% Calculating CR, MSE, PSNR, SSIM:
    
    tic;
    % Calculating CR:
    no_bits = (length(Bout_dcy) + length(Bout_acy) + length(Bout_dcb) + length(Bout_acb) + length(Bout_dcr) + length(Bout_acr))*8;
    CR(k10) = (H*W*8*3)/no_bits;
    
    % Calculating MSE, PSNR, SSIM:
    
    clear Er;
    Er = (A - Crgb).^2;
    mse(k10) = sum(sum(sum(Er)))/(H*W*3);
    psnr(k10) = 10*log10(65025/mse(k10));
    
    Ssim(k10) = 0;
    for k = 1:3
        Ssim(k10) = Ssim(k10) + ssim(A(:,:,k),Crgb(:,:,k));
    end
    Ssim(k10) = Ssim(k10)/3;
    
    % Compressed Size in kB:
    Siz(k10) = no_bits/(8*1024);
    
    % Size of image:
    Szfull(k10,1:2) = [H W];
    toc
end

xlswrite(OPXFileName,Szfull);
sheet=2;
xlswrite(OPXFileName,Siz',sheet);
sheet=3;
xlswrite(OPXFileName,Tcomp',sheet);
sheet=4;
xlswrite(OPXFileName,Tdecomp',sheet);
% sheet=5;
% xlswrite(OPXFileName,CRe',sheet);
sheet=5;
xlswrite(OPXFileName,CR',sheet);
sheet=6;
xlswrite(OPXFileName,mse',sheet);
sheet=7;
xlswrite(OPXFileName,psnr',sheet);
sheet=8;
xlswrite(OPXFileName,Ssim',sheet);