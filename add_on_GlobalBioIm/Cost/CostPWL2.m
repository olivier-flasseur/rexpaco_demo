classdef CostPWL2 < Cost
    % CostPWL2: Patch Weighted L2 norm cost function
    % $$C(\\mathrm{x}) := \sum_n P_n \\frac12\\|\\mathrm{x} - \\mathrm{y}\\|^2_W = \\frac12 (\\mathrm{x} - \\mathrm{y})^T W (\\mathrm{x} - \\mathrm{y}) $$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param W: precision matrix :class:`LinOp` object or scalar
    % :param Mask: binary mask defining where to extract the REXPACO patches
    % :param patch_halfsize: half-size of the REXPACO patches in pixels
    % :param y: data vector
    %
    % **Example** C=CostPWL2(W,Mask,patch_halfsize,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOp`
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        W; % precision matrix
        Mask; % binary mask defining where to extract the PACO patches
        patch_halfsize; % half-size of the PACO patches in pixels
        R; % number of lines of the reconstructed object
        C; % number of columns of the reconstructed object
        T; % number of temporal frames
    end;
    
    %% Constructor
    methods
        function this = CostPWL2(varargin)
            this.W = varargin{1};
            this.Mask = varargin{2};
            this.patch_halfsize = varargin{3};
            this.y = varargin{4};
            this.R = size(this.y,1);
            this.C = size(this.y,2);
            this.T = size(this.y,3);
            this.sizein = [this.R,this.C,this.T];
            this.sizeout = 1;
            this.name = 'CostPWL2';
            this.isConvex = true;
            this.isSeparable = false;
            this.isDifferentiable = true;
        end;
    end;
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
	methods (Access = protected)
        
        function cost=apply_(this,x)
            cost = 0;
            res = this.y - x;
            for r=1:this.R
                for c=1:this.C
                    if(this.Mask(r,c))
                        r_patch_extractor = r-this.patch_halfsize:r+this.patch_halfsize;
                        c_patch_extractor = c-this.patch_halfsize:c+this.patch_halfsize;       
                        current_res = reshape(res(r_patch_extractor,c_patch_extractor,:),[(2*this.patch_halfsize+1)^2,this.T]);
                        w_res = this.W{r}{c} * current_res;
                        cost = cost + 0.5*current_res(:)'*w_res(:);
                    end;
                end;
            end;
        end;
        
        function grad=applyGrad_(this,x)
            grad = zeros(this.sizein);
            res = x - this.y;
            for r=1:this.R
                for c=1:this.C
                    if(this.Mask(r,c))
                        r_patch_extractor = r-this.patch_halfsize:r+this.patch_halfsize;
                        c_patch_extractor = c-this.patch_halfsize:c+this.patch_halfsize;
                        current_res = reshape(res(r_patch_extractor,c_patch_extractor,:),[(2*this.patch_halfsize+1)^2,this.T]);
                        w_res = this.W{r}{c} * current_res;
                        grad(r_patch_extractor,c_patch_extractor,:) = grad(r_patch_extractor,c_patch_extractor,:) + reshape(w_res,[2*this.patch_halfsize+1,2*this.patch_halfsize+1,this.T]);
                    end;
                end;
            end;
        end;

    end
    
    methods (Access = protected)
         %% Copy
         function this = copyElement(obj)
             this = copyElement@Cost(obj);
             this.W = copy(obj.W);          
         end
    end  
end
