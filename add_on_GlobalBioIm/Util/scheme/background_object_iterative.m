function [xopt_save,A_xopt_save,m_save,criterion_save] = background_object_iterative(master_cube,mask_cube,x_init,DIRECT_MODEL,reg_terms,patch_halfsize,opti_tol,opti_maxiter,bckgd_tol,bckgd_maxiter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ##GOAL##
% Function reconstructing a flux image of the sources from an ADI stack.
% To prevent photometric bias due to the contribution of the sources in the
% background statistics, the function alternates between the estimation of 
% the background statistics and the reconstruction of the source flux.
%
% ##INPUTS##
% master_cube: ADI stack.
% mask_cube: binary mask defining where to extract the REXPACO patches.
% x_init: initialization of the reconstructed flux image at first iteration. 
% DIRECT_MODEL: REXPACO forward model.
% reg_terms: regularization terms.
% patch_halfsize: half-size of the REXPACO patches in pixels.
% opti_tol: convergence criterion for optimization.
% opti_maxiter: maximum number of iterations for optimization.
% bckgd_tol: convergence criterion for updating background statistics / source flux.
% bckgd_maxiter: maximum number of iterations for updating background statistics / source flux.
%
% ##OUTPUTS##
% The following outputs are saved for each iteration background statistics / source flux;
% xopt_save: reconstructed source flux image ($\widehat{x}$).
% A_xopt_save: convolded version of xopt_save ($H \widehat{x}$).
% m_save: mean component of the background.
% criterion_save: convergence criterion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % parameters
    [R_pad,C_pad] = size(x_init);
    [~,~,T] = size(master_cube);

    % initialization
    current_cleaned_master_cube = master_cube;
    xopt_old = x_init;
    flag_stop = 0;
    current_iter = 0;
    
    % outputs
    xopt_save = [];
    A_xopt_save = [];
    m_save = [];
    criterion_save = [];
    
    % main loop
    while(~flag_stop && current_iter <= bckgd_maxiter)
        % current iteration
        current_iter = current_iter + 1;
        
        % compute background statistics and data centering
        [current_m,current_Cov_inv] = compute_background_statistics(current_cleaned_master_cube,mask_cube,patch_halfsize);
        m_save = cat(3,m_save,current_m);
        
        % -- Cost definition: Patch Weighted L2 norm cost function
        PWL2 = CostPWL2(current_Cov_inv,mask_cube,patch_halfsize,master_cube-repmat(current_m,[1,1,T]));            
        F = PWL2*DIRECT_MODEL;
        F.doPrecomputation=1;
        
        % -- VMLMB PWL2 + NonNeg + hyperbolicTV + L1norm
        ObjectiveFunction = F;
        for rr=1:numel(reg_terms)
            ObjectiveFunction = ObjectiveFunction + reg_terms{rr};
        end;
        VMLMB = OptiVMLMB(ObjectiveFunction,zeros(R_pad,C_pad),[]);
        VMLMB.verbose = true;
        VMLMB.ItUpOut = 1;
        VMLMB.maxiter = opti_maxiter;
        VMLMB.CvOp=TestCvgCombine('CostRelative',opti_tol(1),'StepRelative',opti_tol(2));
        VMLMB.m = 3;
        VMLMB.OutOp.computecost = true;
        
        VMLMB.run(xopt_old);
        
        % current results
        xopt = VMLMB.xopt;
        A_xopt = DIRECT_MODEL*xopt;
        
        % save
        xopt_save = cat(3,xopt_save,xopt);
        A_xopt_save = cat(3,A_xopt_save,A_xopt(:,:,1));

        % testing convergence
        criterion_save(current_iter) = max(abs(xopt_old(:)-xopt(:))./abs(xopt_old(:)+eps))*100;
        if((criterion_save(current_iter)<bckgd_tol) || ((current_iter>=3)&&(criterion_save(current_iter)-criterion_save(current_iter-1))>(criterion_save(current_iter-2)-criterion_save(current_iter-1))) )
            flag_stop = 1;
        end;
        
        % update
        current_cleaned_master_cube = master_cube - A_xopt;
        xopt_old = xopt;
    end;

end

