function structB = mergeStructures(structA, structB)
% Add & replace fields from structA to structB.
f = fieldnames(structA);
for i = 1:length(f)
    structB.(f{i}) = structA.(f{i});
end
end