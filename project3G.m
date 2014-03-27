N = 256;
a=-3; b=3;
h=(b-a)/N;
r = zeros(1,N);
r(1:5) = -[-205/72, 8/5, -1/5, 8/315, -1/560]./h^2;
D2 = toeplitz(r);
D2(1:4,:) = zeros(4, N);
D2(end-3:end, :) = zeros(4,N);

X = zeros(N,N);
for i=1:N
    X(i,i)=1 - (a+(i-1)*h)^2 - exp(-16*((a+(i-1)*h)^2));
    Vp(i) = exp(-16*((a+(i-1)*h)^2));
end

[V, L] = eig(D2, X);

l = diag(L, 0);

lambda = sort(real(sqrt(l)))

norm = trapz(V(:,2).^2)*h

dE = trapz(V(:,2)*Vp*V(:,2)/norm)*h