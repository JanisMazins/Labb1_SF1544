function J = Jacobian_transfer_functions(k1,k2)
%% Function that takes spring constants k1 and k2 as input parameters and returns 
% the Jacobian of the vector function F = [Fk; Fs]5 from the file transfer_functions.m.
% The Jacobian is a 2x2 matrix organized as J = [dFk/dk1 dFk/dk2; dFs/dk1 dFs/dk2], 
% where dFi/dkj is the partial derivative of Fi with respect to the spring constant kj. 

% Masses m1,m2, damping constants c1,c2.
% Feather constants k1,k2 come as input

m1=475; m2=53; c1=310; c2=1200;
k1_ref=5400; k2_ref=135000; 

vel = 65/3.6;

% renaming the constants
M = m1;
m = m2;
B = c1;
C = c2;
K = k1;
k = k2;
Kref = k1_ref;
kref = k2_ref;

% Frequency comes from the speed
omega = vel;


%% compute the reference values of F1 and F2
%D_ref = (m.*s.^2 + C.*s + kref).*(M.*s.^2 + B.*s + Kref) + M.*s.^2.*(B.*s + Kref);

% compute F1_ref (Fk_ref)
real_D_ref = (M .* m).*omega.^4 - (m.*Kref + C.*B + M.*kref + M.*Kref).*omega.^2 + Kref.*kref;
imag_D_ref = -omega.^3.*(B.*m + C.*M + M.*B) + omega.*(C.*Kref + B.*kref);

%real_f1_ref = -omega.^2.*C.*B + kref.*Kref;
%imag_f1_ref = omega.*(C.*Kref + B.*kref);

real_F1_bar_ref = M.*m.*kref.*Kref.*omega.^4 - M.*m.*C.*B.*omega.^6 + (m.*Kref + C.*B + M.*kref + M.*Kref).*C.*B.*omega.^4 ...
    - (m.*Kref + C.*B + M.*kref + M.*Kref).*kref.*Kref.*omega.^2 + (-omega.^2.*C.*B + kref.*Kref).*kref.*Kref + omega.^2.*(C.*Kref + B.*kref).^2 ...
    - omega.^4.*(C.*Kref + B.*kref).*(B.*m + C.*M + M.*B);

imag_F1_bar_ref = (C.*Kref + B.*kref).*M.*m.*omega.^5 - (C.*Kref + B.*kref).*(m.*Kref + M.*kref + M.*Kref).*omega.^3 + (-omega.^2.*C.*B + kref.*Kref).*(B.*m + C.*M + M.*B).*omega.^3;

F1_tilde_ref = real_F1_bar_ref.^2 + imag_F1_bar_ref.^2;
D_tilde_ref = real_D_ref.^2 + imag_D_ref.^2;
F1_ref = sqrt(F1_tilde_ref)./D_tilde_ref;


% compute F2_ref (Fs_ref)

real_F2_bar_ref = -m.^2.*M.^2.*omega.^8 + omega.^6.*M.*m.*(M + m).*Kref + (m.*Kref + C.*B + M.*kref + M.*Kref).*omega.^6.*m.*M - (m.*Kref + C.*B + M.*kref + M.*Kref).*omega.^4.*(M + m).*Kref ...
    - m.*M.*omega.^4.*Kref.*kref + omega.^2.*(M + m).*Kref.^2.*kref - omega.^6.*(M+m).*B.*(B.*m + C.*M + M.*B) + omega.^4.*(C.*Kref + B.*kref).*(M+m).*B;
    
imag_F2_bar_ref = omega.^7.*(M+m).*M.*m.*B - omega.^5.*(M+m).*B.*((M+m).*Kref + C.*B + M.*kref) + omega.^3.*(M+m).*B.*Kref.*kref ...
    - omega.^7.*m.*M.*(B.*m + C.*M + M.*B) + omega.^5.*m.*M.*(C.*Kref + B.*kref) - omega.^3.*(m+M).*Kref.*(C.*Kref + B.*kref) ...
    + omega.^5.*(B.*m + C.*M + M.*B).*(M + m).*Kref;

F2_tilde_ref = real_F2_bar_ref.^2 + imag_F2_bar_ref.^2;

F2_ref = sqrt(F2_tilde_ref)./D_tilde_ref;


% compute dF1_dk, dF1_dK
dreal_F1_bar_dk = M.*m.*K.*omega.^4 - (M+m).*K.^2.*omega.^2 - 2.*M.*k.*K.*omega.^2 + 2*K.^2.*k + 2*omega.^2.*k.*B.^2 - B.^2.*omega.^4.*(M+m);
dreal_F1_bar_dK = M.*m.*k.*omega.^4 - omega.^2.*M.*k.^2 - 2*omega.^2.*(M+m).*K.*k + 2*K.*k.^2 + 2*omega.^2.*C.^2.*K - omega.^4.*M.*C.^2;
dimag_F1_bar_dk = omega.^5.*B.*M.*m - 2*omega.^3.*B.*M.*k;
dimag_F1_bar_dK = omega.^5.*C.*M.*m - 2*omega.^3.*C.*(M+m).*K;

dreal_D_dk = - M.*omega.^2 + K;
dreal_D_dK = - (M + m).*omega.^2 + k;
dimag_D_dk = omega.*B;
dimag_D_dK = omega.*C;

real_F1_bar = M.*m.*k.*K.*omega.^4 - M.*m.*C.*B.*omega.^6 + (m.*K + C.*B + M.*k + M.*K).*C.*B.*omega.^4 ...
    - (m.*K + C.*B + M.*k + M.*K).*k.*K.*omega.^2 + (-omega.^2.*C.*B + k.*K).*k.*K + omega.^2.*(C.*K + B.*k).^2 ...
    - omega.^4.*(C.*K + B.*k).*(B.*m + C.*M + M.*B);

imag_F1_bar = (C.*K + B.*k).*M.*m.*omega.^5 - (C.*K + B.*k).*(m.*K + M.*k + M.*K).*omega.^3 + (-omega.^2.*C.*B + k.*K).*(B.*m + C.*M + M.*B).*omega.^3;

F1_tilde = real_F1_bar.^2 + imag_F1_bar.^2;

real_D = (M .* m).*omega.^4 - (m.*K + C.*B + M.*k + M.*K).*omega.^2 + K.*k;
imag_D = -omega.^3.*(B.*m + C.*M + M.*B) + omega.*(C.*K + B.*k);

D_tilde = real_D.^2 + imag_D.^2;

real_D = (M .* m).*omega.^4 - (m.*K + C.*B + M.*k + M.*K).*omega.^2 + K.*k;
imag_D = -omega.^3.*(B.*m + C.*M + M.*B) + omega.*(C.*K + B.*k);

dF1_tilde_dk = 2*real_F1_bar.*dreal_F1_bar_dk + 2*imag_F1_bar.*dimag_F1_bar_dk;
dF1_tilde_dK = 2*real_F1_bar.*dreal_F1_bar_dK + 2*imag_F1_bar.*dimag_F1_bar_dK;
dD_tilde_dk = 2*real_D.*dreal_D_dk + 2*imag_D.*dimag_D_dk;
dD_tilde_dK = 2*real_D.*dreal_D_dK + 2*imag_D.*dimag_D_dK;

dF1_dk = (1/2*(F1_tilde).^(-1/2).*dF1_tilde_dk.*D_tilde - sqrt(F1_tilde).*dD_tilde_dk)./(D_tilde.^2);

dF1_dK = (1/2*(F1_tilde).^(-1/2).*dF1_tilde_dK.*D_tilde - sqrt(F1_tilde).*dD_tilde_dK)./(D_tilde.^2);

% compute dF2_dk, dF2_dK
dreal_F2_bar_dk = m.*M.^2.*omega.^6 - m.*M.*omega.^4.*K - M.*(M + m).*omega.^4.*K + omega.^2.*(M+m).*K.^2 + omega.^4.*B.^2.*(M+m);
dreal_F2_bar_dK = 2*omega.^6.*M.*m.*(M + m) - omega.^4.*m.*M.*k - omega.^4.*(M+m).*(2*(M+m).*K + M.*k) + 2*(M+m).*K.*k.*omega.^2;
dimag_F2_bar_dk = - omega.^5.*B.*M.^2;
dimag_F2_bar_dK = - omega.^5.*(M+m).^2.*B + omega.^5.*m.*M.*C - 2*omega.^3.*(M+m).*C.*K + (M+m).*(B.*m + C.*M + M.*B).*omega.^5;

real_F2_bar = -m.^2.*M.^2.*omega.^8 + omega.^6.*M.*m.*(M + m).*K + (m.*K + C.*B + M.*k + M.*K).*omega.^6.*m.*M - (m.*K + C.*B + M.*k + M.*K).*omega.^4.*(M + m).*K ...
    - m.*M.*omega.^4.*K.*k + omega.^2.*(M + m).*K.^2.*k - omega.^6.*(M+m).*B.*(B.*m + C.*M + M.*B) + omega.^4.*(C.*K + B.*k).*(M+m).*B;
    
imag_F2_bar = omega.^7.*(M+m).*M.*m.*B - omega.^5.*(M+m).*B.*((M+m).*K + C.*B + M.*k) + omega.^3.*(M+m).*B.*K.*k ...
    - omega.^7.*m.*M.*(B.*m + C.*M + M.*B) + omega.^5.*m.*M.*(C.*K + B.*k) - omega.^3.*(m+M).*K.*(C.*K + B.*k) ...
    + omega.^5.*(B.*m + C.*M + M.*B).*(M + m).*K;

F2_tilde = real_F2_bar.^2 + imag_F2_bar.^2;

dF2_tilde_dk = 2*real_F2_bar.*dreal_F2_bar_dk + 2*imag_F2_bar.*dimag_F2_bar_dk;
dF2_tilde_dK = 2*real_F2_bar.*dreal_F2_bar_dK + 2*imag_F2_bar.*dimag_F2_bar_dK;

dF2_dk = (1/2*(F2_tilde).^(-1/2).*dF2_tilde_dk.*D_tilde - sqrt(F2_tilde).*dD_tilde_dk)./(D_tilde.^2);

dF2_dK = (1/2*(F2_tilde).^(-1/2).*dF2_tilde_dK.*D_tilde - sqrt(F2_tilde).*dD_tilde_dK)./(D_tilde.^2);

% Scale the partial derivatives with F1_ref and F2_ref
dF1_dK = dF1_dK/F1_ref;
dF1_dk = dF1_dk/F1_ref;
dF2_dK = dF2_dK/F2_ref;
dF2_dk = dF2_dk/F2_ref;

% Jacobian
J = [dF1_dK dF1_dk; dF2_dK dF2_dk];

