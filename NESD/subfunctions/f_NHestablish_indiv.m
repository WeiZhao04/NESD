function [NH] = f_NHestablish_indiv(distbook,TP)
dist = 0.1;
distbook = cellfun(@abs,distbook,'UniformOutput',false);
Gs = size(distbook,2);
%%
for ks = 1:Gs
    tempbook = cell(Gs,1);
    for kk = 1:Gs
        if isempty(tempbook{kk})
            tempbook{kk} = double(distbook{kk,ks}); %tril
        end
    end
    for kk = 1:Gs
        if isempty(tempbook{kk})
            tempbook{kk} = double(distbook{ks,kk}'); % triu 
        end
    end
    tempbook{ks} = zeros(TP(ks),TP(ks)); % add up zeros

%     NH.SZ{ks} = cellfun(@size,tempbook,'UniformOutput',false);
    tempbook = cell2mat(tempbook);
    tempbook(tempbook<dist) = 0;
    
    [tempbook,sInd] = sort(tempbook,'descend');
    nh_index = cell(TP(ks),1);
    [nh_weight,sind_w] = sort(sum(tempbook));
    temp = smooth(nh_weight);
    grad = (temp(end)-temp(1))/TP(ks);
    th = find((diff(temp)-grad)<0,1,'first');
%     sort(sind_w(th:end))
    for kt = sort(sind_w(th:end))
        % neigborhood ruling out 
        NHtemp = tempbook(1:10*Gs,kt);
        cutoff = find(NHtemp>0,1,'last');
        NumNH = sInd(1:cutoff,kt);
        nh_index{kt} = sort(NumNH);
    end
    NH.INX{:,ks} = nh_index;
end

NH.TP = TP;



% reshape and reorder 

for kr = 1:Gs
    NH.KI{kr} = find(cellfun(@isempty,NH.INX{:,kr})==0);
    NH.Len(kr) = length(NH.KI{kr});
end
% [NH.Len,rind] = sort(NH.Len,'descend');
% NH.INX = NH.INX(:,rind);
% NH.KI = NH.KI(:,rind);

end




%%
% search NH size with criterion of group size and sum-of-dist

function sumK = f_assemNH(D,Gs)
sumK = 0;
cutoff = find(D>0,1,'last');
D = D(1:cutoff);
if ~isempty(cutoff)
    if cutoff>floor(Gs*0.5) && sum(D)>1
        sumK = cutoff;
    else
        return
    end
end



end











