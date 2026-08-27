function ppc_val = ppc(inputArg1)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
C           = nchoosek(1:size(inputArg1,1),2);
cos_matrix  = cos(inputArg1(C(:,1),:)-inputArg1(C(:,2),:));
ppc_val     = mean(cos_matrix);

end