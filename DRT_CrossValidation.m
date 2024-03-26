function [CVE,r_real,r_imag] = DRT_CrossValidation(z_exp,t,K_adj,D_adj,...
        para,type_kernel,M,type_weight)

    
    % v2
    l_now = para(1);
    Rshift_now = para(2);
    z_exp = z_exp - Rshift_now; % data shifted by the intercept resistance

    % weighting matx
    scale_vector = z_exp;
    [W_re,W_im] = DRT_Weight(scale_vector,M,type_weight); % take the vector to be the scale
    
    %% INVERSION ( y -> Q ) BASED ON REAL
    % prescribe the t-values on which the inverse distribution is evaluated
    [r_real] = DRT_Inversion(z_exp,t,K_adj,W_re,W_im,D_adj,...
                                l_now,1,0,type_kernel);
    % calculate CVE_imag, using q_real
    CVE_imag = norm(W_im*(imag(z_exp) - imag(K_adj*r_real)));
    
    %% INVERSION BASED ON IMAGINARY
    [r_imag] = DRT_Inversion(z_exp,t,K_adj,W_re,W_im,D_adj,...
                                l_now,0,1,type_kernel);
    % calculate CVE_real, using q_imag
    CVE_real = norm(W_re*(real(z_exp) - real(K_adj*r_imag)));
    
    %% Total CVE
    CVE  = CVE_imag + CVE_real;
        
end