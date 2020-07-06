classdef LinOpAdiTransmission <  LinOp
    % LinOpAdiTransmission: Operator modeling the intensity transmission of
    % the coronagraph
    %
    % :param transmission_map: ADI transmission map of the coronagraph
    % :param T: number of temporal dates
    %
    % All attributes of parent class :class:`LinOp` are inherited.
    %
    % **Example** H=LinOpAdiTransmission(transmission maps)
    %
    % See also :class:`LinOp`, :class:`Map`

    properties (SetAccess = protected,GetAccess = public)
        transmission_map;      % ADI transmission map of the coronagraph
        T;                     % Number of temporal frames
    end
    
    %% Constructor
    methods
        function this = LinOpAdiTransmission(varargin)

            %====  Read arguments         
            transmission_map = varargin{1};
            T = varargin{2};

            %===== Set name, dimensions etc...
            if(~isempty(transmission_map))
                assert(ismatrix(transmission_map),'Parameter transmission_map should be a vector');
            end;
            if(~isempty(T))
                assert(isfinite(T) & T==floor(T),'Parameter T should be an integer');
            end;
            this.transmission_map = transmission_map;
            this.T = T;
            this.name = 'LinOpAdiTransmission';      
            this.isInvertible = false;
            this.isDifferentiable = true;            
            this.sizein = [size(transmission_map,1),size(transmission_map,2),T];
            this.sizeout = this.sizein;
           
        end
    end
    
    %% Core Methods containing implementations (Protected)
    methods (Access = protected)
        
        function y = apply_(this,x)
            y = x.*repmat(this.transmission_map,[1,1,this.T]);
        end;
        
        function y = applyAdjoint_(this,x)
            y = x.*repmat(conj(this.transmission_map),[1,1,this.T]);
        end;
        
    end
    
end
