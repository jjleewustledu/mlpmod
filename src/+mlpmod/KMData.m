classdef KMData < handle & mlsystem.IHandle
    %% Time ~ seconds; pTAC ~ Bq/cc
    %  
    %  Created 12-Jan-2024 17:30:33 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlpmod/src/+mlpmod.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        C
        fqfp
        header
        Hill_d_min
        input_func
        regions
        scanner_data
        timesStart
        timesEnd
    end

    properties (Dependent)
        filepath
        fqfn
        taus
        timesMid
    end

    methods %% GET
        function g = get.filepath(this)
            g = myfileparts(this.fqfp);
        end
        function g = get.fqfn(this)
            g = this.fqfp + ".nii.gz";
        end
        function g = get.taus(this)
            g = this.timesEnd - this.timesStart;
        end
        function g = get.timesMid(this)
            g = (this.timesStart + this.timesEnd)/2;
        end
    end

    methods
        function ic = imaging_metab_corr_ptac(this)
            T = this.table_metab_corr_ptac();
            ic = mlfourd.ImagingContext2(asrow(T.pTAC));
            ic.fqfp = this.fqfp + "_metab_corr_ptac";
            j = struct( ...
                'KMData_header', this.header, ...
                'timesMid', T.Time + 0.5, ...
                'taus', ones(size(T.Time)), ...
                'times', T.Time);
            ic.json_metadata = j;
        end
        function ic = imaging_total_ptac(this)
            T = this.table_total_ptac();
            ic = mlfourd.ImagingContext2(asrow(T.pTAC));
            ic.fqfp = this.fqfp + "_total_ptac";
            j = struct( ...
                'KMData_header', this.header, ...
                'timesMid', T.Time + 0.5, ...
                'taus', ones(size(T.Time)), ...
                'times', T.Time);
            ic.json_metadata = j;
        end
        function ic = imaging_tac(this, varargin)
            T = this.table_tac(varargin{:});
            ic = mlfourd.ImagingContext2(asrow(T.TAC));
            ic.fqfp = this.fqfp + "_total_ptac";
            j = struct( ...
                'KMData_header', this.header, ...
                'timesMid', this.timesMid, ...
                'taus', this.taus, ...
                'times', this.timesStart);
            ic.json_metadata = j;
        end
        function T1 = table_metab_corr_ptac(this)
            %% plasma time activity curves, w/ metab corrections; interpolated to 1 sec

            if ~isempty(this.T_metab_corr_ptac_)
                T1 = this.T_metab_corr_ptac_;
                return
            end

            if isempty(this.T_parent_frac_Hill_)
                T1 = this.table_total_ptac();
                return
            end

            T = this.table_total_ptac();
            Time = T.Time;
            PF = interp1(this.T_parent_frac_Hill_.Time, this.T_parent_frac_Hill_.ParentFraction, Time, ...
                "pchip", this.T_parent_frac_Hill_.ParentFraction(end));
            pTAC = T.pTAC .* PF;

            T1 = table(Time, pTAC);
            this.T_metab_corr_ptac_ = T1;
        end
        function T1 = table_tac(this, roi)
            arguments
                this mlpmod.KMData
                roi {mustBeTextScalar} = 'WholeBrain'
            end

            T = this.scanner_data;
            Time = ascol(this.timesMid);
            vnames = T.Properties.VariableNames;
            assert(any(contains(vnames, roi, IgnoreCase=true)))
            if all(contains(vnames, "[Bq/cc]"))
                TAC = T.(sprintf("tac.%s[Bq/cc]", lower(roi)));
                T1 = table(Time, TAC);
                return
            end
            if all(contains(vnames, "[kBq/cc]"))
                TAC = 1000*T.(sprintf("tac.%s[kBq/cc]", lower(roi)));
                T1 = table(Time, TAC);
                return
            end
            if all(contains(vnames, "[muCi/cc]"))
                TAC = 37000*T.(sprintf("tac.%s[muCi/cc]", lower(roi)));
                T1 = table(Time, TAC);
                return
            end
            if all(contains(vnames, "[uCi/cc]"))
                TAC = 37000*T.(sprintf("tac.%s[muCi/cc]", lower(roi)));
                T1 = table(Time, TAC);
                return
            end
            error("mlpmod:ValueError", stackstr())
        end
        function T1 = table_total_ptac(this)
            %% plasma time activity curves, w/o metab corrections; interpolated to 1 sec

            if ~isempty(this.T_total_ptac_)
                T1 = this.T_total_ptac_;
                return
            end

            T = this.input_func;
            T.Properties.VariableNames = ["Time", "pTAC"];
            switch this.input_func.Properties.VariableNames{1}
                case 'sample-time[minutes]'
                    T.Time = T.Time*60;
                case 'sample-time[seconds]'
                otherwise
                    error("mlpmod:ValueError", stackstr())
            end
            switch this.input_func.Properties.VariableNames{2}
                case 'plasma.1[Bq/cc]'
                case 'plasma.1[kBq/cc]'
                    T.pTAC = 1000*T.pTAC;
                case {'plasma.1[muCi/cc]', 'plasma.1[uCi/cc]'}
                    T.pTAC = 37000*T.pTAC;
                otherwise
                    error("mlpmod:ValueError", stackstr())
            end

            Time = ascol(0:this.timesMid(end)); % interpolated to 1 sec
            pTAC = ascol(interp1(T.Time, T.pTAC, Time, ...
                "pchip", T.pTAC(end)));
            T1 = table(Time, pTAC);
            this.T_total_ptac_ = T1;
        end
    end

    methods (Static)
        function this = create(fqfn, opts)
            %% CREATE is a factory method for KMData.
            %
            % Args:
            %     fqfn {mustBeFile} : e.g., *kmData.xls
            %     opts.T_parent_frac table = table()
            %     opts.Hill_d_min double = 0.95 % Hill's parameter d, minimum

            arguments
                fqfn {mustBeFile}
                opts.T_parent_frac table = table()
                opts.Hill_d_min double = 0.95 % Hill's parameter d, minimum
            end
            this = mlpmod.KMData;
            this.C = readcell(fqfn);
            this.fqfp = myfileprefix(fqfn);
            this.Hill_d_min = opts.Hill_d_min;

            % parse kmData structures
            data_rows = find(cellfun(@(x) strcmp(x, "# DATA"), this.C));
            comment_rows = find(cellfun(@(x) strcmp(x, "# NON_PMOD_COMMENTS"), this.C));
            regions_row = find(cellfun(@(x) strcmp(x, "// regions"), this.C));
            last_row = size(this.C, 1);
            last_col = size(this.C, 2);

            % extract kmData structures
            this.header = cell2table(this.C(1:data_rows(1)-1, 1:3), ...
                VariableNames={'Key', 'Value1', 'Value2'});
            this.input_func = cell2table(this.C(data_rows(1)+2:comment_rows(1)-1, 1:2), ...
                VariableNames=strrep(this.C(data_rows(1)+1, 1:2), "// ", ""));
            this.regions = this.C(regions_row, 3:last_col);
            this.scanner_data = cell2table(this.C(data_rows(2)+2:last_row, 3:last_col), ...
                VariableNames=this.C(data_rows(2)+1, 3:last_col));
            this.timesStart = asrow(cell2mat(this.C(data_rows(2)+2:last_row, 1)));
            this.timesEnd = asrow(cell2mat(this.C(data_rows(2)+2:last_row, 2)));

            % manage T_parent_frac_Hill_
            if ~isempty(opts.T_parent_frac)
                this.T_parent_frac_Hill_ = this.build_Hill(opts.T_parent_frac);
            end

            %are_missing = cellfun(@(x) isa(x, "missing"), this.input_func);
        end
    end

    %% PRIVATE

    properties (Access = private)
        T_metab_corr_ptac_
        T_parent_frac_Hill_
        T_total_ptac_        
    end

    methods (Access = private)
        function this = KMData()
        end

        function [h,T1] = build_Hill(this, T)
            T.Properties.VariableNames = ["Time", "ParentFraction"];
            if max(T.Time) <= 190 
                % ensure sec
                T.Time = T.Time*60; 
            end
            if max(T.ParentFraction) > 1
                % ensure frac \in [0,1]
                T.ParentFraction = 0.01*T.ParentFraction;
            end

            hill = mlwong.Hill(T, dt=1, d_min=this.Hill_d_min, visualize_anneal=false);
            T1 = hill.solve();
            T1.Properties.VariableNames = ["Time", "ParentFraction"];
            h = plot(hill);
            saveFigure2(h, this.fqfp+"_Hill");
            writetable(T1, this.fqfp+"_Hill.csv")
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
