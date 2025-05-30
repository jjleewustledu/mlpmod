classdef Pmod < handle
    %% line1
    %  line2
    %  
    %  Created 21-May-2024 17:35:31 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlpmod/src/+mlpmod.
    %  Developed on Matlab 24.1.0.2603908 (R2024a) Update 3 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        xls_filename
    end

    properties (Dependent)
        age_acquisition
        date_acquisition
        date_birth
        ifc_template
        kmdata_table
        weight  % kg
        dose  % mCi
    end

    methods %% GET
        function g = get.age_acquisition(this)
            g = years(this.date_acquisition - this.date_birth);
        end
        function g = get.date_acquisition(this)
            row1 = contains(this.kmdata_table.ColumnA, "PATIENT_ACQUISITION_DATE");
            g = datetime(this.kmdata_table{row1, 2}, InputFormat="yyyyMMdd");
        end
        function g = get.date_birth(this)
            row1 = contains(this.kmdata_table.ColumnA, "PATIENT_BIRTH_DATE");
            g = datetime(this.kmdata_table{row1, 2}, InputFormat="yyyyMMdd");
        end
        function g = get.dose(this)
            row1 = contains(this.kmdata_table.ColumnA, "PATIENT_ACTIVITY_PRE");
            dose1 = double(this.kmdata_table{row1, 2});
            row2 = contains(this.kmdata_table.ColumnA, "PATIENT_ACTIVITY_POST");
            dose2 = double(this.kmdata_table{row2, 2});
            g = dose1 - dose2;
        end
        function g = get.ifc_template(~)
            g = mlfourd.ImagingFormatContext2( ...
                fullfile(mlwong.TZ3108.home, "ifc_template.nii.gz"));
        end
        function g = get.kmdata_table(this)
            if isempty(this.kmdata_table_)
                this.kmdata_table_ = this.readtable();
            end
            g = this.kmdata_table_;
        end
        function g = get.weight(this)
            row = contains(this.kmdata_table.ColumnA, "PATIENT_WEIGHT");
            g = double(this.kmdata_table{row, 2});
        end
    end

    methods
        function this = Pmod(opts)
            arguments
                opts.filename {mustBeTextScalar} = ""
            end
            this.xls_filename = opts.filename;
            warning("off", "MATLAB:table:ModifiedAndSavedVarnames");
        end

        function u = aif(this)
            t = this.kmdata_table;
            top = find(t.ColumnA == "// sample-time[minutes]" & t.ColumnB == "plasma.1[kBq/cc]");
            bottom = find(t.ColumnA == "# NON_PMOD_COMMENTS");
            u = t(top+1:bottom-1, 1:2);
            u.Properties.VariableNames = ["time_min", "activity_kBq_cc"];
            u.time_min = double(u.time_min);
            u.activity_kBq_cc = double(u.activity_kBq_cc);

            % plot(u, "time_min", "activity_kBq_cc")
        end

        function u = aif_over_auc(this)
            %% activity units ~ (kBq/cm^3)/(kBq*min/cm^3) ~ 1/min
            u = this.aif();
            auc = trapz(u.time_min, u.activity_kBq_cc);
            activity_inv_min = u.activity_kBq_cc/auc;
            time_min = u.time_min;
            u = table(time_min, activity_inv_min);
        end

        function v = aif_as_suv(this, opts)
            %% activity units ~ (kBq/cm^3)/(kBq/g) ~ g/cm^3

            arguments
                this mlpmod.Pmod
                opts.interp_method = []
                opts.dt = 1/60  % min
                opts.T = []
            end

            u = this.aif();
            dose_kBq = 37e3*this.dose;
            weight_g = 1e3*this.weight;
            time_min = u.time_min;
            suv = u.activity_kBq_cc*weight_g/dose_kBq;
            v = table(time_min, suv);

            if isempty(opts.T)
                opts.T = max(u.time_min);
            end
            if ~isempty(opts.interp_method)
                time_min = ascol(0:opts.dt:opts.T);
                suv = ascol(interp1(v.time_min, v.suv, time_min, opts.interp_method, v.suv(end)));
                v = table(time_min, suv);
            end
        end

        function v = aif_as_suvr(this, opts)
            %% activity units ~ ptac / tac

            arguments
                this mlpmod.Pmod
                opts.interp_method = []
                opts.dt = 1/60  % min
                opts.T = []
            end

            u = this.aif();
            time_min = u.time_min;

            T_wb = sum(this.tac_wholebrain.taus);
            timesMid = this.tac_wholebrain.timesMid;
            tac = this.tac_wholebrain.tac_wholebrain_kBq_cc_;
            wholebrain_time_avg = trapz(timesMid, tac)/T_wb;

            suvr = u.activity_kBq_cc/wholebrain_time_avg;
            v = table(time_min, suvr);

            if isempty(opts.T)
                opts.T = max(u.time_min);
            end
            if ~isempty(opts.interp_method)
                time_min = ascol(0:opts.dt:opts.T);
                suvr = ascol(interp1(v.time_min, v.suvr, time_min, opts.interp_method, v.suvr(end)));
                v = table(time_min, suvr);
            end
        end

        function v = aif_rescaled(this, opts)
            %% activity units ~ kBq/cm^3
            %  rescaled by new dose (mCi) and new weight (kg) for repurposing an animal's AIF for closely related studies

            arguments
                this mlpmod.Pmod
                opts.new_dose double = this.dose
                opts.new_weight double = this.weight
                opts.interp_method = []
                opts.dt = 1/60  % min
                opts.T = []
            end

            u = this.aif();
            dose_adj = opts.new_dose/this.dose;
            weight_adj = opts.new_weight/this.weight;
            time_min = u.time_min;
            activity_kBq_cc = u.activity_kBq_cc*weight_adj/dose_adj;
            v = table(time_min, activity_kBq_cc);

            if isempty(opts.T)
                opts.T = max(u.time_min);
            end
            if ~isempty(opts.interp_method)
                time_min = ascol(0:opts.dt:opts.T);
                suvr = ascol(interp1(v.time_min, v.suvr, time_min, opts.interp_method, v.suvr(end)));
                v = table(time_min, suvr);
            end
        end

        function [aif,tacs] = build_nifti(this)
            ifc = this.ifc_template;
            aif = copy(ifc);
            aif.img = asrow(this.aif.activity_kBq_cc);
            aif.fqfp = strrep(myfileprefix(this.xls_filename), "-pkin.kmData", "-aif");
            aif.fileprefix = strrep(aif.fileprefix, "_proc-", "_trc-tz3108_proc-");
            aif.filepath = strrep(aif.filepath, mlwong.TZ3108.home, mlwong.TZ3108.bids_home);
            aif.filepath = strrep(aif.filepath, "chemistry", "pet");
            Nt = size(this.aif.time_min, 1);
            aif.json_metadata = struct( ...
                "times", 60*this.aif.time_min, ...
                "timesMid", 60*this.aif.time_min, ...
                "taus", ones(Nt, 1));
            aif = mlfourd.ImagingContext2(aif);

            tacs = copy(ifc);
            arr = table2array(this.tacs);
            tacs.img = arr(:, 4:end)';
            tacs.fqfp = strrep(myfileprefix(this.xls_filename), "-pkin.kmData", "-tacs");
            tacs.fileprefix = strrep(tacs.fileprefix, "_proc-", "_trc-tz3108_proc-");
            tacs.filepath = strrep(tacs.filepath, mlwong.TZ3108.home, mlwong.TZ3108.bids_home);
            tacs.filepath = strrep(tacs.filepath, "chemistry", "pet");
            tacs.json_metadata = struct( ...
                "times", this.tacs.times, ...
                "timesMid", this.tacs.timesMid, ...
                "taus", this.tacs.taus);
            tacs = mlfourd.ImagingContext2(tacs);
        end

        function r = regions(this)
            u = this.tacs();
            r = string(u.Properties.VariableNames);
            r = r(4:end);
        end

        function u = tac_wholebrain(this)
            t = this.tacs();
            times = t.times;
            timesMid = t.timesMid;
            taus = t.taus;
            try
                tac_wholebrain_kBq_cc_ = t.tac_wholebrain_kBq_cc_;
                u = table(times, timesMid, taus, tac_wholebrain_kBq_cc_);
            catch ME
                handexcept(ME)
            end
        end

        function u = tacs(this)
            t = this.kmdata_table;
            top = find(t.ColumnA == "// start[seconds]" & t.ColumnB == "end");
            u = t(top+1:end, 3:end);  % already double
            select = contains(u.Properties.VariableNames, "kBq_cc");
            u = u(:, select);

            timings = double(t{top+1:end, 1:2});  % cast to double
            times = timings(:,1);  % a.k.a. starts
            timesMid = mean(timings(:,1:2), 2);
            taus = timings(:,2) - timings(:,1);
            u = addvars(u, times, timesMid, taus, Before=1, NewVariableNames=["times", "timesMid", "taus"]);
        end

        function t = readtable(this, opts)
            %% https://chatgpt.com/c/b95ff87f-23da-4aae-9cc9-26f2cc58809c

            arguments
                this mlpmod.Pmod
                opts.filename {mustBeTextScalar} = this.xls_filename
            end
            this.xls_filename = opts.filename;

            % Define import options
            xls_opts = detectImportOptions(opts.filename);

            % Set the variable types for the first two columns to 'text'
            xls_opts.VariableTypes(1:2) = {'char', 'char'};

            % Set the variable names to be read from row row_vars
            xls_opts.DataRange = 'A1';  % MANUAL
            row_vars = this.row_vars(filename=opts.filename);
            xls_opts.VariableNamesRange = sprintf('%i:%i', row_vars, row_vars);  % '69:69';

            % Set the rest of the columns to be numeric
            numCols = width(xls_opts.VariableTypes);
            xls_opts.VariableTypes(3:numCols) = repmat({'double'}, 1, numCols - 2);

            % Read the table
            t = readtable(opts.filename, xls_opts);

            % Adjustments MANUAL
            t.Properties.VariableNames{1} = 'ColumnA';
            t.Properties.VariableNames{2} = 'ColumnB';
            t.ColumnA = string(t.ColumnA);
            t.ColumnB = string(t.ColumnB);            

            % Display the first and last rows of the table to verify
            this.kmdata_table_ = t %#ok<NOPRT>
        end

        function r = row_vars(this, opts)
            arguments
                this mlpmod.Pmod
                opts.filename {mustBeTextScalar} = this.xls_filename
            end
            this.xls_filename = opts.filename;

            c = readcell(this.xls_filename);
            c = c(:,1);
            r = find(strcmp(c, "// start[seconds]"));
        end
    end

    methods (Static)
    end

    %% PRIVATE

    properties (Access = private)
        kmdata_table_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
