function [stoich,bofMets] = getBOFstoich(model)
biomassReactionID = find(model.c == 1);

bofMetIDs = find(model.S(:,biomassReactionID)<0);

bofMets = model.mets(bofMetIDs);

stoich = model.S(bofMetIDs,biomassReactionID);

end