classdef LinOpAdiMotion <  LinOp
    % LinOpAdiMotion: Operator modeling the apparent motion of the field of
    % view of ADI series along time. It accounts for the interpolation in
    % the adjoint operator.
    %
    % :param rotation_vector: ADI rotation vector of parallactic angles
    % :param R_C_padded: vector of spatial dimensions of the images
    %
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    % **Example** H=LinOpAdiMotion(rotation_vector,R_C_padded)
    %
    % See also :class:`LinOp`, :class:`Map`

    properties (SetAccess = protected,GetAccess = public)
        rotation_vector;       % Parallactic angles
        T;                     % Number of temporal frames
        Zrot;                  % Interpolation operator for each temporal frame
    end
    
    %% Constructor
    methods
        function this = LinOpAdiMotion(varargin)

            %====  Read arguments         
            rotation_vector = varargin{1};
            R_C_padded = varargin{2};

            %===== Set name, dimensions etc...
            if(~isempty(rotation_vector))
                assert(isvector(rotation_vector),'Parameter rotation_vector should be a vector');
            end;
            if(~isempty(R_C_padded))
                assert(isvector(R_C_padded),'Parameter R_C_padded should be a vector');
            end;
            this.rotation_vector = rotation_vector;
            this.T = numel(rotation_vector);
            this.name = 'LinOpAdiMotion';      
            this.isInvertible = false;
            this.isDifferentiable = true;            
            this.sizein = R_C_padded;
            this.sizeout = [R_C_padded, this.T];
            %===== Build interpolation operators
            [X,Y] = meshgrid(1:this.sizein(1),1:this.sizein(2));
            this.Zrot = {};
            for t=1:this.T
                theta = (this.rotation_vector(1)-this.rotation_vector(t))/180*pi;
                rotation_matrix = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
                XYq = rotation_matrix*([X(:)-this.sizein(2)/2 Y(:)-this.sizein(1)/2]');
                Xq = reshape(XYq(1,:),[this.sizein(1),this.sizein(2)])+this.sizein(2)/2;
                Yq = reshape(XYq(2,:),[this.sizein(1),this.sizein(2)])+this.sizein(1)/2;
                this.Zrot{t} = LinOpInterp2(X,Y,Xq,Yq);
            end;
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        
        function y = apply_(this,x)
            y = zeros(this.sizeout);
            for t=1:this.T
                y(:,:,t) = this.Zrot{t}*x;
            end;
        end;
        
        function y = applyAdjoint_(this,x)            
            y = zeros(this.sizein);
            for t=1:this.T
                y = y + this.Zrot{t}'*x(:,:,t);
            end;
        end;
        
    end
    
end
