function [Y] = f_GSL_appr_indiv(neighbors,Esvd)

%%

N = size(Esvd.vec,2);
P = size(Esvd.vec{1},1);
T = sum(neighbors.TP);
%%
% computing weights with neighbors
for kn = 1:N
    NHs = neighbors.INX{:,kn};
    KI = neighbors.KI{kn};
    NoK = neighbors.Len(kn);
    Ytemp = zeros(P,NoK);
    Coef = zeros(neighbors.TP(kn),NoK);
    Lat = zeros(1,NoK);
    for kt = 1:NoK
        if size(NHs{KI(kt)},1) > ceil(N/2)
            Bin = false(T,1);
            Bin(NHs{KI(kt)}) = true;
            Bin = mat2cell(Bin,neighbors.TP);
            comps = [];
            for ks = 1:N %
                Ind = find(Bin{ks}>0);
                if ~isempty(Ind)
                    Rind = mod(Ind,neighbors.TP(ks))+1;
                    comps = [comps Esvd.vec{ks}(:,Rind)];
                end
            end
            centraComp = Esvd.vec{kn}(:,KI(kt));
            w = mldivide(comps,centraComp);
            Ytemp(:,kt) = comps*w; % projection with weights
            Coef(:,kt) = Esvd.cof{kn}(:,KI(kt));
            Lat(kt) = Esvd.lat{kn}(kt);
        else
            Ytemp(:,kt) = Esvd.vec{kn}(:,KI(kt));
            Coef(:,kt) = Esvd.cof{kn}(:,KI(kt));
            Lat(kt) = Esvd.lat{kn}(kt);
        end
    end
    % reform the subject-specific Emaps
    Y.vec{kn} = Ytemp;
    Y.cof{kn} = Coef;
    Y.lat{kn} = Lat;
end



end


