
x = linspace(1,2,20);
y = 2*x.^2+1;
polyfit(x,y,2)

x1 = [1 ; 2; 3]
y1 = [1; 2; 3]
z1 = x1+y1;
surffit = fit([x1, y1],z1,'poly23','normalize','on')
%{
g_sqrt_11 = 1;
g_sqrt_12 = P1/2-(R1/R2)*(R2^2-(x2-P2/2)^2)^0.5/x1_t_minus;
g_sqrt_13 = 1;

g_sqrt_21 = P2/2-(R2/R1)*(R1^2-(x1-P1/2)^2)^0.5/x2_t_minus;
g_sqrt_22 = ( 1/(x1_t_plus-x1_t_minus)*(x2_t_plus-x2_t_minus) )*...
    ( 4*(R1^2-(x1-P1/2)^2)^0.5*(R2^2-(x2-P2/2)^2)^0.5+...
    (2*x1-x1_t_plus-x1_t_minus)*(2*x2-x2_t_plus-x2_t_minus)*...
    (x1-P1/2)*(x2-P2/2)/((R1^2-(x1-P1/2)^2)^0.5*(R2^2-(x2-P2/2)^2)^0.5) );
g_sqrt_23 = ( (R2/R1)*(R1^2-(x1-P1/2)^2)^0.5 - P2/2 )/(x2_t_plus - P2);

g_sqrt_31 = 1;
g_sqrt_32 = ( (R1/R2)*(R2^2-(x2-P2/2)^2)^0.5 - P1/2)/(x1_t_plus-P1);
g_sqrt_33 = 1;
%}