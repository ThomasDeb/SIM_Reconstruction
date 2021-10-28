classdef CostTV1D < Cost
    % CostMixNormSchatt1 Mixed Schatten-l1 Norm [1]
    % $$C(\\mathrm{x}) :=   \\sum_n  \\| \\mathrm{x}_{n\\cdot} \\|_{Sp}, $$
    % for \\(p \\in [1,+\\infty]\\). Here, \\(\\|\\cdot\\|_{Sp}\\)  denotes the p-order Shatten norm
    % defined by
    % $$ \\|\\mathrm{X}\\|_{Sp} = \\left[\\sum_k (\\sigma_k(\\mathrm{X}))^p\\right]^{1/p},$$
    % where \\(\\sigma_k(\\mathrm{X})\\) is the k-th singular value of \\(\\mathrm{X}\\). In other words it is the lp-norm
    % of the signular values of  \\(\\mathrm{X}\\).
    %
    % :param p: order of the Shatten norm (default 1)
    %
    % All attributes of parent class :class:`Cost` are inherited.
    %
    % **Note** The actual implementation works for size (sz) having one of the two following forms:
    %
    %   * (NxMx3) such that the Sp norm will be applied on each symetric 2x2
    %     $$ \\begin{bmatrix} \\mathrm{x}_{n m 1} & \\mathrm{x}_{n m 2} \\newline
    %     \\mathrm{x}_{n m 2} & \\mathrm{x}_{n m 3} \\end{bmatrix}$$
    %     and then the \\(\\ell_1\\) norm on the two other dimensions.
    %   * (NxMxKx6) such that the Sp norm will be applied on each symetric 3x3
    %     $$ \\begin{bmatrix} \\mathrm{x}_{n m k 1} & \\mathrm{x}_{n m k 2} & \\mathrm{x}_{n m k 3} \\newline
    %     \\mathrm{x}_{n m k 2} & \\mathrm{x}_{n m k 4} & \\mathrm{x}_{n m k 5} \\newline
    %     \\mathrm{x}_{n m k 3} & \\mathrm{x}_{n m k 5} & \\mathrm{x}_{n m k 6} \\newline  \\end{bmatrix}$$
    %     and then the \\(\\ell_1\\) norm on the three other dimensions.
    %
    % **References**
    % [1] Lefkimmiatis, S., Ward, J. P., & Unser, M. (2013). Hessian Schatten-norm regularization
    % for linear inverse problems. IEEE transactions on image processing, 22(5), 1873-1888.
    %
    % **Example** C=CostMixNormSchatt1(sz,p,y)
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOp`

    %%    Copyright (C) 2017
    %     E. Soubies emmanuel.soubies@epfl.ch
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

    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        index; % Index along which T1 TV is performed
    end

    %% Constructor
    methods
        function this = CostTV1D(sz, index, y)
            % Verify if the mexgl files exist
            if (exist('TV1D_denoise_mex')~=3)
                [mpath,~,~] = fileparts(which('TV1D_denoise_mex.c'));
                pth = cd;
                cd(mpath);
                mex -v TV1D_denoise_mex.c -largeArrayDims;
                cd(pth);
            end

            if nargin<3, y=0; end
            this@Cost(sz,y);
            this.name='CostTV1D';
            this.isConvex=true;
            this.isDifferentiable=false;
            assert(index <= length(sz),'Index must be smaller than the dimension of the input');
            this.index = index;
        end
    end

    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)

        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            S1.type = '()'; S2.type = '()';
            subs1 = repmat({':'}, 1, length(this.sizein));
            subs2 = subs1; subs1{this.index} = 2:this.sizein(this.index);
            subs2{this.index} = 1:this.sizein(this.index)-1;
            S1.subs = subs1; S2.subs = subs2;
            diff = subsref(x, S1) - subsref(x, S2);
            y = sum(abs(diff(:)));
        end
        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            sz = this.sizein;
            non_index = 1:length(sz); non_index(this.index) = [];
            y = zeros(sz);
            S.type = '()';
            subs = cell(1, length(sz));
            subs{this.index} = ':'; S.subs = subs;
            indices = cell(1, length(non_index));
            for i = 1 : prod(sz(non_index)) % Could be parallelized
                [indices{:}] = ind2sub(sz(non_index), i);
                S.subs(non_index) = indices;
                y = subsasgn(y, S, TV1D_denoise_mex(subsref(x, S), alpha));
            end
        end
    end
end
