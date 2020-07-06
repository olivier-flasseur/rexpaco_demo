classdef RegCostL1 < Cost
    % RegCostL1: L1 norm cost function
    % $$C(x) := \\|\\mathrm{x}\\|_1 $$
    % 
    % :param sz: size of the input
    %
    % :param flag_constraint: flag for the flag_input (can be a matrix of
    %  size sz to apply the constraint on each element of x)
    %   ->  1: positivity constraint on x
    %   -> -1: negativity constraint on x
    %   ->  0: no constraint on x
    %
    % :param W: the weighting matrix (default 1)
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Example** C=RegCostL1(sz)
    %
    % **Example** C=RegCostL1(sz, flag_constraint)
    %
    % **Example** C=RegCostL1(sz, flag_constraint, W)
    %
    % See also :class:`Map` :class:`Cost`, :class:`LinOp`
    
    %%    Copyright (C) 2015
    %     F. Soulez ferreol.soulez@epfl.ch
    %     Modified: 06/05/2018 (mm/dd/yyyy)
    %     Anthony Berdeu (Laboratoire Hubert Curien)
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
    
    %% Constructor
    properties
        flag_constraint ;
        W ;
    end
    methods
        function this = RegCostL1(sz, flag_constraint, W)
            this.name='RegCostL1';
            this.isConvex=true;
            
            if nargin < 3 || isempty(W)
                this.W = 1 ;
            else
                this.W = W ;
            end
            
            if nargin < 2 || isempty(flag_constraint)
                flag_constraint = 0 ;
            elseif isempty(flag_constraint)
                flag_constraint = 0 ;
            else
                flag_constraint = sign(double(flag_constraint)) ;
            end
            
            this.sizein = sz ;
            this.flag_constraint = flag_constraint ;
            
            if any(flag_constraint(:)==0)
                this.isDifferentiable = false ;
            else
                this.isDifferentiable = true ; % The cost is differentiable
                    % on its definition domain
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)
% % % %         TODO: put W
% % % %         function y=applyProx_(~,x,alpha)
% % % %             % Reimplemented from parent class :class:`Cost`.
% % % %             % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) = \\mathrm{sign(x)} \\mathrm{max}(\\vert x \\vert- \\alpha,0)$$
% % % %             y =  sign(x) .* max( abs( x) - alpha,0);
% % % %         end
        
        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            if this.isDifferentiable
                if any(x(:).*this.flag_constraint(:)<0)
                    y = + inf;
                else
                    y = x.*this.W ;
                    y=sum(abs(y(:)));
                end
            else
                y = x.*this.W ;
                y=sum(abs(y(:)));
            end
        end

        function g=applyGrad_(this,~)
            % Reimplemented from parent class :class:`Cost`.
            if ~this.isDifferentiable
                error(['L_1 norm not differentiable if no positivity', ...
                    ' or negativity constraint is applied...']) ;
            end
            g = this.W.*ones(this.sizein).*this.flag_constraint ;
        end
    end
end
