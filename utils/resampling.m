function i = resampling(q, type)
% 1 = multinomial, 2 = systematic, 3 = stratified

M = length(q);

if(type == 1) % -- Multinomial
    u = rand(M,1);
    qc = cumsum(q);
    qc = qc(:);
    qc=qc/qc(M);
    [~,ind1]=sort([u;qc]);
    ind2=find(ind1<=M);
    i=ind2'-(0:M-1);
elseif(type == 2) % -- Systematic
    qc = cumsum(q);
    u = ((0:M-1)+rand(1))/M;
    i = zeros(1,M); k = 1;
    for j = 1:M
        while(qc(k)<u(j))
            k = k + 1;
        end
        i(j) = k;
    end
elseif(type == 3) % -- Stratified
    u=([0:M-1]'+(rand(M,1)))/M;
    qc=cumsum(q);
    qc=qc(:);
    [~,ind1]=sort([u ; qc]);
    ind2=find(ind1<=M);
    i = ind2'-(0:M-1);
else
    error('No such resampling type');
end