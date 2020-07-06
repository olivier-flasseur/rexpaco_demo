function [m,Cov_inv] = compute_background_statistics(master_cube,mask_cube,patch_halfsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ##GOAL##
% Function computing the background statistics (local estimation of the mean
% and covariances from multi-variate Gaussian models).
%
% ##INPUTS##
% master_cube: ADI stack.
% mask_cube: binary mask defining where to extract the REXPACO patches.
% patch_halfsize: half-size of the REXPACO patches in pixels.
%
% ##OUTPUTS##
% m: mean component.
% Cov_inv: precision matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % warning off
    warning off;
    
    % compute mean
    m = mean(master_cube,3);
    
    % compute inverse of shrunk covariances
    for r=1:size(master_cube,1)
        for c=1:size(master_cube,2)
            if(mask_cube(r,c))
                r_patch_extractor = r-patch_halfsize:r+patch_halfsize;
                c_patch_extractor = c-patch_halfsize:c+patch_halfsize;
                current_samples = squeeze(master_cube(r_patch_extractor,c_patch_extractor,:));
                current_mean = squeeze(m(r_patch_extractor,c_patch_extractor));
                current_samples = reshape(current_samples - repmat(current_mean,[1,1,size(current_samples,3)]),[size(current_samples,1)*size(current_samples,2),size(current_samples,3)]);
                S = (current_samples*current_samples')./size(current_samples,2);
                [current_Cov,~] = shrinkage_covariance(S,size(current_samples,2));
                Cov_inv{r}{c} = current_Cov^(-1);
            end;            
        end;
    end;

    % warning on
    warning on;
    
end

