function ind = stratified_resampling_MtoN(q,N)

qc = cumsum(q);
u = ((0:N-1)+rand(1,N))/N;

[~,ind1]=sort([u(:) ; qc(:)]);
ind2=find(ind1<=N);
ind = ind2'-(0:N-1);

end