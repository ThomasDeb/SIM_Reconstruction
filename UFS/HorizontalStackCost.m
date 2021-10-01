classdef HorizontalStackCost < Cost
    properties (SetAccess = protected,GetAccess = public)
        costs;      % Cell array of costs
        idx;        % Dimensions corresponding to input of the costs
        costsSz;    % Shape of costs array (in case of repetition of the same cost)
    end

    %% Constructor
    methods
        function this = HorizontalStackCost(sz, costs, idx, costsSz)
            if nargin < 4
                costsSz = size(costs);
            end
            assert(isequal(sz(idx), costs{1}.sizein), ...
                'idx should be compatible with input size of costs');
            assert(isequal(sz(setdiff(1:length(sz), idx)), costsSz), ...
                'remaining indices should be compatible with array of costs');

            for i = 1 : length(costs)
                assert(isequal(costs{i}.sizein, costs{1}.sizein), 'All costs should have the same input size')
            end

            this@Cost(sz);
            this.costs = costs;
            this.idx = idx;
            this.costsSz = costsSz;

        end
    end

    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyProx_(this,x,alpha)
    methods (Access = protected)

        function y=apply_(this,x)
            % Reimplemented from parent class :class:`Cost`.
            y = 0;
            S.type = '()';
            subs = cell(this.sizein);
            for i = 1 : length(this.idx), subs{this.idx(i)} = ':'; end
            S.subs = subs;
            other_idx = setdiff(1:length(this.sizein), this.idx);
            if length(this.costs) == 1
                for i = 1 : prod(this.costsSz)
                    S.subs{other_idx} = ind2sub(this.costsSz, i);
                    y = y + this.costs{1}.apply(subsref(x, S));
                end
            else
                for i = 1 : prod(this.costsSz)
                    S.subs{other_idx} = ind2sub(this.costsSz, i);
                    y = y + this.costs{i}.apply(subsref(x, S));
                end
            end
        end

        function y=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            y = zeros(size(x));
            S.type = '()';
            subs = cell(this.sizein);
            for i = 1 : length(this.idx), subs{this.idx(i)} = ':'; end
            S.subs = subs;
            other_idx = setdiff(1:length(this.sizein), this.idx);
            if length(this.costs) == 1
                for i = 1 : prod(this.costsSz)
                    S.subs{other_idx} = ind2sub(this.costsSz, i);
                    y = subsasgn(y, S, this.costs{1}.applyProx(subsref(x, S)));
                end
            else
                for i = 1 : prod(this.costsSz)
                    S.subs{other_idx} = ind2sub(this.costsSz, i);
                    y = subsasgn(y, S, this.costs{i}.applyProx(subsref(x, S), alpha));
                end
            end
        end
    end
end
