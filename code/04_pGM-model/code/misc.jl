# c_H,    q_H,    r_H    = 0.1, +1, 1.97e-10; # Levin (2009)
# c_ClO4, q_ClO4, r_ClO4 = 0.1, -1, 2.83e-10; # Levin (2009)

# function U_beta1(z, ah, kappa, elc, beta, epsilon_o, epsilon_w)
#     if z < 1e10ah # while inside hydrated ion radius
#         U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z) - 3.05
#     else # outside of radius
#         U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z)
#     end
# end;
# U_H(t) = U_beta1(t, r_H, kappa, elc, beta, epsilon_o, epsilon_w)
# W_H    = W(q_H, r_H);
# H      = (c_H, q_H, U_H);


# function U_beta2(z, ah, kappa, W)
#     if z < 1e10ah # while inside hydrated ion radius
#         U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah)) - 2.1
#     else # outside of radius
#         U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah))
#     end
# end;
# W_ClO4 = W(q_ClO4, r_ClO4);
# U_ClO4(t) = U_beta2(t, r_ClO4, kappa, W_ClO4)
# ClO4 = (c_ClO4, q_ClO4, U_ClO4);