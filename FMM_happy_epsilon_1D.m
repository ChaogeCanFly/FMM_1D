function [eps] = FMM_happy_epsilon_1D(epsilon, N, M)

eps_m = fft(epsilon)/M;
e1 = eps_m(1:(2*N+1));
e2 = eps_m((M-2*N):M);

%we have to renumerate e2 for Toeplitz matrix, because matlab
%writes eps(m) with negative m in reverse order
for i=1:N
    b = e2(i+1);
    e2(i+1) = e2(2*N+1-i+1);
    e2(2*N+1-i+1)=b;
end
e2(1)=e1(1);

eps = toeplitz(e1,e2);