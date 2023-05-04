function structB = addToStructure(structA, structB)
% Add fields in structA to structB, if they don't exist already.
f = fieldnames(structA);
for i = 1:length(f)
    
    if ~isfield(structB, f{i})
        structB.(f{i}) = structA.(f{i});
    end
end
end