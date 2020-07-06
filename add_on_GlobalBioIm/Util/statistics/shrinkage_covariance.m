function [shrinked_Cov,rho] = shrinkage_covariance(Cov,nb_samples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ##GOAL##
% Function regularizing covariance matrices by shrinkage.
%
% ##INPUTS##
% Cov: sample covariance matrix.
% nb_samples: number of samples used to compute Cov.
%
% ##OUTPUTS##
% skrinked_Cov: shrunk covariance matrix.
% rho: shrinkage factor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % display shrunk covariance
    flag_display = 0;
    
    % shrinkage procedure
    F = diag(diag(Cov));
    
    trace_Cov_Cov = trace(Cov*Cov);
    sum_square_diag_Cov = sum(diag(Cov).^2);
    
    rho = (trace_Cov_Cov + (trace(Cov))^2 -2*sum_square_diag_Cov) / ((nb_samples+1)*(trace_Cov_Cov-sum_square_diag_Cov));
    rho = min(max(rho,0),1);
    
    shrinked_Cov = (1 - rho)*Cov + rho*F;

    if(flag_display)
        figure; 
        subplot(1,3,1); imshow(Cov,[]); colormap(jet(1024)); title('$\widehat{\textbf{S}}$','Interpreter','latex');
        subplot(1,3,2); imshow(F,[min(Cov(:)) max(Cov(:))]); colormap(jet(1024)); title('$\widehat{\textbf{F}}$','Interpreter','latex');
        subplot(1,3,3); imshow(shrinked_Cov,[min(Cov(:)) max(Cov(:))]); colormap(jet(1024)); title('$\widehat{\textbf{C}}$','Interpreter','latex');
    end;
    
end
