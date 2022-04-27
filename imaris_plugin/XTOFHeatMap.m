%
%  Computes motility heatmap based on optical flow
%  
%  Revision: 3RC Build: 201604010906
%
%
%  Author: Pizzagalli Diego Ulisse (1,2)
%          1. Institute for Research in Biomedicine - Bellinzona (CH)
%          2. Institute of Computational Science,
%             Universita della Svizzera Italiana - Lugano (CH)
%
%
%  INSTALLATION ON COMPUTERS WITH MATLAB ALREADY INSTALLED:
%                1. Copy this file into an XTensions folder
%                   (e.g C:\Program Files\Bitplane\Imaris [ver]\XT\matlab).
%                2. Restart Imaris and you can find this function in the 
%                   Image Processing menu with the "SVMColoc" name.
%                NOTE: Tested with MATLAB r2012b, r2013a
%
%
%    <CustomTools>
%      <Menu name="IRB - immi toolbox">
%        <Item name="Heatmap motility intensity" icon="I" tooltip="computes the maximum opticval flow projection in time">
%          <Command>MatlabXT::XTOFHeatMap(%i)</Command>
%        </Item>
%      </Menu>
%    </CustomTools>
% 
% 
%  Copyright 2016 - Diego Ulisse Pizzagalli. 
%  Free to use and distribute under Creative Commmons Licence 3.0
%  http://creativecommons.org/licenses/by/3.0/


function XTOFHeatMap(aImarisApplicationID)
    if isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
        vImarisApplication = aImarisApplicationID;
    else
        javaaddpath ImarisLib.jar;
        vImarisLib = ImarisLib;
        if ischar(aImarisApplicationID)
            aImarisApplicationID = round(str2double(aImarisApplicationID));
        end
        vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
    end
    [dir, fn, ext] = fileparts(mfilename('fullpath'));
    addpath(strcat(dir, filesep, 'libs', filesep, 'opticalflow'));
    
    aDataSet = vImarisApplication.GetDataSet.Clone;
    dataset_size = [aDataSet.GetSizeX, aDataSet.GetSizeY, aDataSet.GetSizeZ, aDataSet.GetSizeC, aDataSet.GetSizeT];
    
    W = dataset_size(1); %556
    H = dataset_size(2); %556
    Z = dataset_size(3); %18
    C = dataset_size(4); %3
    T = dataset_size(5); %60
    
    DS_MAX_X = aDataSet.GetExtendMaxX(); 
    DS_MAX_Y = aDataSet.GetExtendMaxY(); 
    DS_MAX_Z = aDataSet.GetExtendMaxZ(); 
    DS_MIN_X = aDataSet.GetExtendMinX(); 
    DS_MIN_Y = aDataSet.GetExtendMinY(); 
    DS_MIN_Z = aDataSet.GetExtendMinZ();
    
    DS_dx = (DS_MAX_X - DS_MIN_X) / (W-1);
    DS_dy = (DS_MAX_Y - DS_MIN_Y) / (H-1);
    DS_dz = (DS_MAX_Z - DS_MIN_Z) / (Z-1);
    
    max_distance = W*DS_dx*H*DS_dy;
    
    TF = -1;
    
    prompt = {['Enter channel to use [1:',num2str(C),']:']};
    dlgtitle = 'Input';
    dims = [1];
    definput = {'1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput)
    chIdx = str2num(answer{1});
    
    %% Creating output RGB channel
    aDataSet.SetSizeC(C + 1);
    aDataSet.SetChannelName(C, ['Optical Flow intensity - ch',num2str(chIdx)]);
    aDataSet.SetChannelColorRGBA(C, 16711680);
    aDataSet.SetChannelRange(C, 0, max_distance);
    
    %% Processing
    z_stack = zeros(W,H,Z);
    z_stack_old = zeros(W,H,Z);
    
    %multiWaitbar( 'CloseAll' );
    
    v_stack = zeros(W,H,Z);
    u_stack = zeros(W,H,Z);
    a_stack = zeros(W,H,Z);
    Mv_stack = zeros(W,H,Z);
    Mu_stack = zeros(W,H,Z);
    a_stack = zeros(W,H,Z);
    Ma_stack = zeros(W,H,Z);
    
    cc = chIdx - 1;
    for TF=0:T-1
        %multiWaitbar('Frames',TF/(T-1));
        z_stack = medfilt3(aDataSet.GetDataVolumeShorts(cc, TF));
        if(TF == 0)
            z_stack_old = z_stack;
        end
        
        ofInt3d = zeros(W,H,Z);
        tic;
        for cz=1:Z
            %multiWaitbar('Planes',cz/(Z));
            [~, ~, ~] = grad2Dm(z_stack_old(:,:,cz), 1, 1);
            [dy, dx, dt] = grad2Dm(z_stack(:,:,cz), 1);
            [vt1, ut1] = DoFlowT1(dx, dy, dt);
            ofInt3d(:,:,cz) = sqrt(vt1.^2 + ut1.^2).*DS_dx; %magnitude
            a_stack(:,:,cz) = a_stack(:,:,cz)+(ofInt3d(:,:,cz)*(1/T-1));
            v_stack(:,:,cz) = v_stack(:,:,cz)+(vt1*(1/T-1));
            u_stack(:,:,cz) = u_stack(:,:,cz)+(ut1*(1/T-1));
            Mv_stack(:,:,cz) = max(Mv_stack(:,:,cz), vt1); %solo se maggiore della media? update both the components.
            Mu_stack(:,:,cz) = max(Mu_stack(:,:,cz), ut1);
            Ma_stack(:,:,cz) = max(Ma_stack(:,:,cz), ofInt3d(:,:,cz));
        end
        toc;
        ofInt3d(ofInt3d > max_distance-1) = max_distance-1;
        aDataSet.SetDataVolumeFloats(ofInt3d,  C,  TF);
        z_stack_old = z_stack;
    end
    %multiWaitbar('CloseAll');
    
    %% Visualization
    ranges = zeros(C, 2);
    colors = zeros(C, 1, 'uint32');
    rgb_colors = zeros(C, 3, 'uint8');
    for cc = 0:C - 1
        ranges(cc+1, 1) = aDataSet.GetChannelRangeMin(cc);
        ranges(cc+1, 2) = aDataSet.GetChannelRangeMax(cc);
        colors(cc+1) = aDataSet.GetChannelColorRGBA(cc);
        rgb_colors(cc+1,1) = (bitand(colors(cc+1), hex2dec('FF'))/hex2dec('FF'))*255;
        rgb_colors(cc+1,2) = (bitand(colors(cc+1), hex2dec('FF00'))/hex2dec('FF00'))*255;
        rgb_colors(cc+1,3) = (bitand(colors(cc+1), hex2dec('FF0000'))/hex2dec('FF0000'))*255;
    end
    
    z_stack = zeros(W,H,Z,C);
    TF = 0;
    for cc=0:C-1
        z_stack(:,:,:,cc+1) = aDataSet.GetDataVolumeFloats(cc, TF);
    end

    R = zeros(H,W,3);
    for cc=1:C
        min_int = ranges(cc,1);
        max_int = ranges(cc,2);
        m = 255 / (max_int -  min_int);
        q = (-min_int)*m;
        c_temp = (max(z_stack(:,:,:,cc),[],3)*m) + q;
        c_temp(c_temp <= 0) = 0;
        c_temp(c_temp >= 255) = 255;
        c_temp = c_temp';
        
        R(:,:,1) = max(R(:,:,1), (c_temp.*double(rgb_colors(cc, 1)/255)));
        R(:,:,2) = max(R(:,:,2), (c_temp.*double(rgb_colors(cc, 2)/255)));
        R(:,:,3) = max(R(:,:,3), (c_temp.*double(rgb_colors(cc, 3)/255)));
    end

    ds = 4;
    dsx = 1:ds:W;
    dsy = 1:ds:H;
    [DSX, DSY] = meshgrid(dsx, dsy);
    hold on;
    
    figure;
    imshow(flipud((uint8(R)))); %OK, solo flupud perche' gia' trasposto in c_temp
    hold on;
    quiver(DSY, DSX, -flipud((max(Mu_stack(1:ds:end,1:ds:end,:),[],3))'), -flipud((max(Mv_stack(1:ds:end,1:ds:end,:),[],3))'), 0, 'w', 'LineWidth',2);
    title('Maximum optical Flow direction');
    
    figure;
    HM = max(Ma_stack,[],3);
    HM(1:7, :) = 0;
    HM(end-7:end, :) = 0;
    HM(:,1:7) = 0;
    HM(:,end-7:end) = 0;
    HM = imgaussfilt(mat2gray(flipud(HM')),7);
    imagesc(imadjust(HM));
    hold on;
    title('Optical Flow heat map (\sigma = 7)');
    save('hm.mat', 'HM');
    %OK
    
    figure;
    imagesc(imadjust(HM));
    hold on;
    quiver(-flipud(Mu_stack(:,:,3)'), -flipud(Mv_stack(:,:,3)'), 0, 'w', 'LineWidth',2);
    title('max of displacements');

    pause;
    vImarisApplication.SetDataSet(aDataSet);
    close all;
    clear all;
    
    function colIm = flow2colIm(u,v)
        H  = (atan2(v,u)+ pi+0.000001)/(2*pi+0.00001);   
        V = min(0.99999999,sqrt(u.^2 + v.^2)*8/10);
        colIm = hsv2rgb(H,  ones(size(H),'single'),V);
    end
end