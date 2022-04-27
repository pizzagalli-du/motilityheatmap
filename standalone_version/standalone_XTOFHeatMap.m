%  Optical Flow viewm (color)
%  Revision: 4RC Build: 20210602

close all;
clear all;
clc;

[curr_fn, curr_path] = uigetfile('*.ims', 'Select IMS file');

try
    curr_ims = ImarisReader([curr_path, '\', curr_fn]);
catch
    disp(['ims read error ', curr_path, '\', curr_fn]);
    return;
end
min_x = curr_ims.DataSet.ExtendMinX;
min_y = curr_ims.DataSet.ExtendMinY;
min_z = curr_ims.DataSet.ExtendMinZ;

DS_dx = (curr_ims.DataSet.ExtendMaxX - curr_ims.DataSet.ExtendMinX) / curr_ims.DataSet.SizeX;
DS_dy = (curr_ims.DataSet.ExtendMaxY - curr_ims.DataSet.ExtendMinY) / curr_ims.DataSet.SizeY;
DS_dz = (curr_ims.DataSet.ExtendMaxZ - curr_ims.DataSet.ExtendMinZ) / curr_ims.DataSet.SizeZ;
    
W = curr_ims.DataSet.SizeX; %556
H = curr_ims.DataSet.SizeY; %556
Z = curr_ims.DataSet.SizeZ; %18
C = curr_ims.DataSet.SizeC; %3
T = curr_ims.DataSet.SizeT; %60
    
max_distance = W*DS_dx*H*DS_dy;
    
TF = -1;

prompt = {['Enter channel to use [1:',num2str(C),']:'], 'Sigma:'};
dlgtitle = 'Input';
dims = [1];
definput = {'1', '7'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
chIdx = str2num(answer{1});
sigma = str2num(answer{2});


%% Processing
z_stack = zeros(W,H,Z);
z_stack_old = zeros(W,H,Z);

v_stack = zeros(W,H,Z);
u_stack = zeros(W,H,Z);
a_stack = zeros(W,H,Z);
Mv_stack = zeros(W,H,Z);
Mu_stack = zeros(W,H,Z);
a_stack = zeros(W,H,Z);
Ma_stack = zeros(W,H,Z);

cc = chIdx - 1;
for TF=0:T-1
    z_stack = medfilt3(curr_ims.DataSet.GetDataVolume(cc, TF));
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
    z_stack_old = z_stack;
end

%% Visualization
ds = 4;
dsx = 1:ds:W;
dsy = 1:ds:H;
[DSX, DSY] = meshgrid(dsx, dsy);
hold on;

figure;
HM = max(Ma_stack,[],3);
HM(1:7, :) = 0;
HM(end-7:end, :) = 0;
HM(:,1:7) = 0;
HM(:,end-7:end) = 0;
HM = imgaussfilt(mat2gray(flipud(HM')),sigma);
imagesc(imadjust(HM));
axis equal; axis tight;
colorbar;
hold on;
title(['Optical Flow heat map (\sigma = ',num2str(sigma),')']);

function colIm = flow2colIm(u,v)
    H  = (atan2(v,u)+ pi+0.000001)/(2*pi+0.00001);   
    V = min(0.99999999,sqrt(u.^2 + v.^2)*8/10);
    colIm = hsv2rgb(H,  ones(size(H),'single'),V);
end