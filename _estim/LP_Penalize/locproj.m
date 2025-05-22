function [IR, theta, gamma, IR_se] = locproj(y, x, w, H_min, H_max, r, lambda, varargin)
% Function for estimating IRF and coefficients in penalized LP
    % Reference: Smooth/penalized local projection (Barnichon & Brownlees, 2019)
    
    %%% Input %%%
    % y:       response variable
    % x:       impulse variable
    % w:       controls (contemperaneous and lagged)
    % H_min:   minimum horizon
    % H_max:   maximum horizon
    % r:       order of finite difference operator
    % lambda:  penalty strength
    % varargin: Newey-West lag length for standard errors (optional input)
    
    %%% Output %%%
    % IR:    impulse response
    % theta: penalized coef, for Xb. Warning: correspond to b_k in our paper
    % gamma: unpenalized coef, for W. Warning: correspond to \zeta_h and \phi_{h,l} in our paper
    % IR_se: standard errors for IR
    
    % Design data matrix
    w = [ ones(size(w,1),1) w ]; % Add intercept to controls
    [B, Xb, Y, Y_resw, Xb_resw] = locproj_design(y, x, w, H_min, H_max);
    
    % Compute impulse responses and coefficients using partitioned formula
    if isempty(varargin)
        L = H_max; % Default Newey-West lag length
    else
        L = varargin{1};
    end
    [IR, theta, gamma, IR_se] = locproj_partitioned(Y, B, Xb, w, Y_resw, Xb_resw, r, lambda, L);
    
end
