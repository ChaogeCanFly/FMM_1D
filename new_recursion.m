function [S_p] = new_recursion(S_p1, W_p, W_p1, pplus, pminus, N)

W11_p = W_p( 1:(2*N+1), 1:(2*N+1));
W21_p = W_p( ((2*N+2):(4*N+2)), (1:(2*N+1)) );
W12_p = W_p( (1:(2*N+1)), ((2*N+2):(4*N+2)) );
W22_p = W_p( ((2*N+2):(4*N+2)), ((2*N+2):(4*N+2)) );

W11_p1 = W_p1( (1:(2*N+1)), (1:(2*N+1)));
W21_p1 = W_p1( ((2*N+2):(4*N+2)), (1:(2*N+1)) );
W12_p1 = W_p1( (1:(2*N+1)), ((2*N+2):(4*N+2)) );
W22_p1 = W_p1( ((2*N+2):(4*N+2)), ((2*N+2):(4*N+2)) );


Tuu_p1 = S_p1( 1:(2*N+1), 1:(2*N+1));
Rud_p1 = S_p1( 1:(2*N+1), (2*N+2):(4*N+2));
Rdu_p1 = S_p1( (2*N+2):(4*N+2), 1:(2*N+1));
Tdd_p1 = S_p1( (2*N+2):(4*N+2), (2*N+2):(4*N+2));

Rud_p1_w = pplus*Rud_p1/pminus;
Tdd_p1_w = Tdd_p1/pminus;
Tuu_p1_w = pplus*Tuu_p1;

Z11 = W11_p1;
Z12 = -W11_p*Rud_p1_w - W12_p;
Z21 = W21_p1;
Z22 = -W21_p*Rud_p1_w - W22_p;

Z_up = cat(2, Z11, Z12);
Z_down = cat(2, Z21, Z22);
Z = cat(1, Z_up, Z_down);

X11 = W11_p*Tuu_p1_w;
X12 = -W12_p1;
X21 = W21_p*Tuu_p1_w;
X22 = - W22_p1;

X1 = cat(1, X11, X21);
X2 = cat(1, X12, X22);

ZX1 = Z\X1;
ZX2 = Z\X2;

Rud_p = ZX2(1:(2*N+1),:);
Tdd_p = Tdd_p1_w * ZX2((2*N+2):(4*N+2),:) ;
Tuu_p = ZX1(1:(2*N+1),:);
Rdu_p = Rdu_p1 + Tdd_p1_w * ZX1((2*N+2):(4*N+2),:) ;

S_p_up = cat(2, Tuu_p, Rud_p);
S_p_down = cat(2, Rdu_p, Tdd_p);
S_p = cat(1, S_p_up, S_p_down);

