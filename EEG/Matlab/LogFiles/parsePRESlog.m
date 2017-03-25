function stats=parsePRESlog(data,conv,savit)
% PARSEPRESLOG parses data read from a Presentation log and returns stats

if nargin<3, savit=false;   end

stats=struct;
dnames=fieldnames(data);
switch conv
    case 'gestalt'
        conds={'letter' 'street' 'ykanji'};
        bloks={'b1','b2','b3','b4'};
        for i=1:numel(dnames)
            data.(dnames{i}).type;
        end
end

if savit,   save(stats);    end
end % parsePRESlog