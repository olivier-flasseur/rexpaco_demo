classdef LinOpInterp2 < LinOp
    % LinOpInterp: Linear interpolation operator on a cartesian grid for
    % dimensions lower than 3
    % 
    % :param coord_in: coordinates of the sample points.
    %
    % :param coord_out: coordinates of the query points.
    %
    % :param extrapVal: extrapolation value for the query points which fall
    % outside the intial grid (default = 0, see notes).
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**:
    %   1. it is possible to extrapolate with values different from 0 but
    %   in that case the interpolation operation is not linear.
    %   2. this operator is implemented for a number of dimensions
    % lower than 3
    %   3. for 2D and 3D matrices, by convention the coordinates order is 
    %   y,x and y,x,z
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out)
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out, extrapVal)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Based on the following code:
    %     https://zenodo.org/record/3585632#.XwMaNXUzY5k
    %     Created: 03/27/2018 (mm/dd/yyyy)
    %     Anthony Berdeu (Laboratoire Hubert Curien)
    %
    %     Updated: 07/05/2020 (mm/dd/yyyy) for REXPACO: 
    %     Olivier Flasseur (Observatory of Lyon)
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.

    properties
        ndms ;      % number of dimensions
        weight ;    % weighting coefficients of the linear interpolation
        V_ind       % index of the interpolated values
        inner_ind   % index of the query positions which lie inside the
            % grid
        extrapVal ; % extrapolation value for the query points which fall
            % outside the intial grid.
    end
    
    %% Constructor
    methods
        function this = LinOpInterp2(varargin)
            %% Initialization
            this.name = 'LinOpInterp2' ;             
            this.isInvertible = false ;
            
            % Coordinates of the sample points.
            coord_in{1} = varargin{2}(:,1);
            coord_in{2} = varargin{1}(1,:);
            this.ndms = length(coord_in) ;
            
            % Determination of sizein and sizeout
            this.sizein = zeros(1, this.ndms) ;
            for dim = 1:this.ndms
                this.sizein(dim) = length(coord_in{dim}) ;
            end

            % Query points
            coord_out = [varargin{4}(:) varargin{3}(:)];
            nb_p = size(coord_out,1) ;
            this.sizeout = [size(varargin{3},1) size(varargin{4},2)];

            % Extrapolation value for the query points which fall
            % outside the intial grid.
            if nargin>4
                this.extrapVal = varargin{5} ;
            else
                this.extrapVal = 0 ;
            end
            if this.extrapVal ~= 0
                warning(['Extrapolation value different from 0. ', ...
                    'The operator is not linear...']) ;
            end
            
            %% Determination of the interpolated indexes and weighting
            % coefficients
            
            % Initialization
            this.weight = ones(nb_p, 2*this.ndms) ;
            this.V_ind = ones(nb_p, 2^this.ndms, 'uint32') ;
            this.inner_ind = uint32(1:nb_p) ;
            
            % Loop on the query points
            nb_y = this.sizein(1) ;
            y = coord_in{1} ;
            x = coord_in{2} ;
            for p = 1:nb_p
                % Differential positions
                dy = y-coord_out(p,1) ;
                dx = x-coord_out(p,2) ;
                ind_y_p = find(dy>=0 & coord_out(p,1)>=0,1) ;
                ind_x_p = find(dx>=0 & coord_out(p,2)>=0,1) ;
                
                % Is the query point lying in the grid on the
                % y-axis?
                if ind_y_p > 1
                    ind_y_m = ind_y_p-1 ;
                elseif(isempty(ind_y_p))  % This query point falls
                    % outside the grid
                    this.inner_ind(p) = 0 ;
                else % (ind_y_p==1) && (dx(ind_y_p)==0) -> point on
                    % the left edge
                    ind_y_p = 2 ;
                    ind_y_m = 1 ;
                end
                
                % Is the query point lying in the grid on the
                % x-axis?
                if ind_x_p > 1
                    ind_x_m = ind_x_p-1 ;
                elseif isempty(ind_x_p) % This query point falls
                    % outside the grid
                    this.inner_ind(p) = 0 ;
                else % (ind_x_p==1) && (dx(ind_x_p)==0) -> point on
                    % the left edge
                    ind_x_p = 2 ;
                    ind_x_m = 1 ;
                end
                
                % Determination of the interpolated indexes and
                % weighting coefficients
                if this.inner_ind(p)
                    % weight
                    w_x = dx(ind_x_p)/(dx(ind_x_p)-dx(ind_x_m)) ;
                    w_y = dy(ind_y_p)/(dy(ind_y_p)-dy(ind_y_m)) ;
                    this.weight(p,:) = [ ...
                        w_y, 1-w_y, ...         % y_0 y_1
                        w_x, 1-w_x ] ;          % x_0 x_1
                    
                    % index
                    ind = ind_y_m + (ind_x_m-1)*nb_y ;
                    this.V_ind(p,:) = [ ...
                        ind, ...                % y_0 x_0
                        ind+1, ...              % y_1 x_0
                        ind+nb_y, ...           % y_0 x_1
                        ind+1+nb_y ] ;          % y_1 x_1
                end
            end
            
            % Exclusion of the points lying outside the grid
            this.inner_ind = this.inner_ind(this.inner_ind>0) ;
            this.weight = this.weight(this.inner_ind,:) ;
            this.V_ind = this.V_ind(this.inner_ind,:) ;
        end
    end
	
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
	methods (Access = protected)
        %% apply_
        function y = apply_(this,x)
            y = this.extrapVal*ones(this.sizeout(1)*this.sizeout(2),1) ;
            y(this.inner_ind) = ...
                this.weight(:,1) .* this.weight(:,3) .* x(this.V_ind(:,1)) + ...   % y_0 x_0
                this.weight(:,2) .* this.weight(:,3) .* x(this.V_ind(:,2)) + ...   % y_1 x_0
                this.weight(:,1) .* this.weight(:,4) .* x(this.V_ind(:,3)) + ...   % y_0 x_1
                this.weight(:,2) .* this.weight(:,4) .* x(this.V_ind(:,4));        % y_1 x_1
            y = reshape(y,[this.sizeout(1),this.sizeout(2)]);
        end
        
        %% applyAdjoint_
        function y = applyAdjoint_(this,x)
            x_inner = x(this.inner_ind)';
            
            aux = [this.weight(:,1).*this.weight(:,3).*x_inner ; ...             % y_0 x_0
                this.weight(:,2).*this.weight(:,3).*x_inner ; ...       	     % y_1 x_0
                this.weight(:,1).*this.weight(:,4).*x_inner ; ...          	     % y_0 x_1
                this.weight(:,2).*this.weight(:,4).*x_inner];                  	 % y_1 x_1
            
            y = accumarray(this.V_ind(:), aux, [prod(this.sizein), 1]);
            y = reshape(y, this.sizein);            
        end
        
    end
end
