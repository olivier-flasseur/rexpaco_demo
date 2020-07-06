%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ##GOAL##
% Main file to reconstruct a flux image from an ADI stack with REXPACO.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% flags 
flag_display = 0;
font_size = 25;

%% fixed parameters
patch_halfsize = 4;

%% load data
% dataset 
lambda_vector = fitsread('./add_on_GlobalBioIm/data/dataset/rylup_lambda_vector.fits');
rotation_vector = fitsread('./add_on_GlobalBioIm/data/dataset/rylup_rotation_vector.fits');
master_cube = fitsread('./add_on_GlobalBioIm/data/dataset/rylup_master_cube.fits');
master_psf = fitsread('./add_on_GlobalBioIm/data/dataset/rylup_psf_cube.fits');
[R,C,T,L] = size(master_cube);

% instrument
mask_cube = fitsread('./add_on_GlobalBioIm/data/instrument/sphere_mask_cube.fits');
transmission_cube = fitsread('./add_on_GlobalBioIm/data/instrument/sphere_coronagraph_h2_transmission_cube.fits');

if(flag_display)
    % parallactic rotation
    figure; plot(1:T,rotation_vector(1)-rotation_vector); xlabel('temporal index'); ylabel('relative rotation'); set(gca,'Fontsize',font_size); grid on;

    % master_cube
    figure;
    for tt=1:T
        imshow(mask_cube.*squeeze(master_cube(:,:,tt)),[0 200]); colormap(jet(1024)); colorbar; title(sprintf('$t = %d$',tt),'Interpreter','latex'); set(gca,'Fontsize',font_size); pause();
    end;
    
    % master_psf
    figure;
    for ll=1:L
        imshow(squeeze(master_cube(:,:,1,ll)),[0 200]); colormap(jet); colorbar;
        title(sprintf('$\\lambda = %d$',ll),'Interpreter','latex');
        set(gca,'Fontsize',font_size);
        pause();
    end;
    
    % master_psf
    figure; imshow(master_psf,[]); colormap(jet(1024)); colorbar;

    % transmission_cube
    figure; imshow(transmission_cube,[]); colormap(jet(1024)); colorbar;
end;

%% patch definition
patch_size = 2*patch_halfsize+1;
patch_mask = double(ones(patch_size,patch_size));

mask_cube_holes = zeros(R,C);
mask_cube_holes(1:patch_size:end,1:patch_size:end) = 1;
mask_cube_holes = mask_cube & mask_cube_holes;

mask_cube_holes(1:patch_halfsize,:) = 0; mask_cube_holes(:,1:patch_halfsize) = 0; mask_cube_holes(R-patch_halfsize+1:R,:) = 0; mask_cube_holes(:,C-patch_halfsize+1:C) = 0;

%% compute pad
% define indexes
R_pad = round(sqrt(2)*R);
C_pad = round(sqrt(2)*C);

if(mod(R_pad,2))
    R_pad = R_pad + 1;
end;
if(mod(C_pad,2))
    C_pad = C_pad + 1;
end;
R_delta_pad = (R_pad-R)/2;
C_delta_pad = (C_pad-C)/2;

R_pad_idx_start = R_delta_pad+1;
C_pad_idx_start = C_delta_pad+1;
R_pad_idx_end = R_pad-R_delta_pad;
C_pad_idx_end = C_pad-C_delta_pad;

%% define operators and forward model
% -- To run on GPU (0: CPU / 1: Matlab Parrallel Computing Toolbox / 2: CudaMat)
useGPU(0)

% -- Convolution Operator definition
H = LinOpConv(repmat(fft2(ifftshift(master_psf)),[1,1,T]));
H.doPrecomputation = 1;
checkLinOp(H);

% -- Adi Motion Operator definition
Q = LinOpAdiMotion(rotation_vector,[R_pad,C_pad]);
Q.doPrecomputation = 1;
checkLinOp(Q);

% -- Adi Transmission Operator definition
G = LinOpAdiTransmission(transmission_cube,T);
G.doPrecomputation = 1;
checkLinOp(G);

% -- Adi Field of View Selection
V = LinOpSelectorPatch([R_pad,C_pad,T], [R_pad_idx_start,C_pad_idx_start,1], [R_pad_idx_end,C_pad_idx_end,T], false);
V.doPrecomputation = 1;
checkLinOp(V);

% -- Total forward model
DIRECT_MODEL = V*H*G*Q;
checkLinOp(DIRECT_MODEL);

% -- Regularization definition: hyperbolicTV
mu_hyperbolicTV = 1E5;
G = LinOpGrad([R_pad C_pad]);
hyperB = CostHyperBolic(G.sizeout,1E-7,3)*G;

% -- Regularization definition: CostL1
mu_L1 = 1E5;
costL1 = RegCostL1([R_pad C_pad],true);

% -- Save regularization terms
reg_terms{1} = mu_hyperbolicTV*hyperB;
reg_terms{2} = mu_L1*costL1;

%% compute iteratively the background statistics and the flux image
x_init = zeros(R_pad,C_pad);
opti_tol = [1E-8 1E-8];
opti_maxiter = 100;
bckgd_tol = 0.1;
bckgd_maxiter = 5;

[xopt_save,A_xopt_save,m_save,criterion_save] = background_object_iterative(master_cube,mask_cube_holes,x_init,DIRECT_MODEL,reg_terms,patch_halfsize,opti_tol,opti_maxiter,bckgd_tol,bckgd_maxiter);

%% display
figure; imshow(xopt_save(:,:,end),[]); colormap(hot(1024)); colorbar; hold on; plot((R_pad+1)/2,(C_pad+1)/2,'*','Color',[1,1,1],'MarkerSize',15,'LineWidth',3); title('$\widehat{x}$','Interpreter','latex'); set(gca,'Fontsize',font_size);
figure; imshow(conv2(xopt_save(:,:,end),master_psf,'same'),[]); colormap(hot(1024)); colorbar; hold on; plot((R_pad+1)/2,(C_pad+1)/2,'*','Color',[1,1,1],'MarkerSize',15,'LineWidth',3); title('H$\widehat{x}$','Interpreter','latex'); set(gca,'Fontsize',font_size);

