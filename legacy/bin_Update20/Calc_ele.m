function [e_mat] = Calc_ele(A,c,p)
%Calculate relative concentration of element for each phase
ny    = size(c{1}{1},1);
nx    = size(c{1}{1},2);
for ip = 1:size(p,3)
    c_ary                 = zeros(length(c{ip}),nx*ny);
    for ic = 1:length(c{ip})
        c_ary(ic,:)       = c{ip}{ic}(:)';
    end
    %Calculate element concentration
    e_ary                 = A{ip}'*[c_ary;1-sum(c_ary,1)];
    %Normalize element concentration
    e_ary                 = e_ary./repmat(sum(e_ary,1),size(e_ary,1),1);
    %Dependent element
    e_ary                 = e_ary(1:end-1,:);
    for ie = 1:size(e_ary,1)
        e_mat{ip}{ie}     = reshape(e_ary(ie,:),ny,nx);
    end
end
end