function [B, Xb, Y, Y_resw, Xb_resw] = locproj_design(y, x, w, H_min, H_max)
% Function for designing data matrix in penalized LP

    %%% Basic %%%
    % T:  sample size
    % K:  number of B-spline basis functions
    % HR: number of horizons
    % p:  number of controls
    
    %%% Input %%%
    % y:       response variable
    % x:       impulse variable
    % w:       controls (contemperaneous and lagged)
    % H_min:   minimum horizon
    % H_max:   maximum horizon
    
    %%% Output %%%
    % B:       B-spline values for each basis function at each horizon
    % Xb:      x * B
    % Y:       responses and leads thereof
    % Y_resw:  Y residualized by w
    % Xb_resw: Xb residualized by w
    
    % prepare
    T  = length(y);
    HR = H_max + 1 - H_min;

    % construct the B-spline basis functions
    B = bspline( (H_min:H_max)' , H_min , H_max+1 , H_max+1-H_min , 3 );
    K = size( B , 2 );

    % placeholder for output matrices
    Xb = nan(T,K,HR);
    Y = nan(T,HR);
    Y_resw = nan(T,HR);
    Xb_resw = nan(T,K,HR);

    % go thru each horizon
    for ih=1:HR

        h = H_min-1+ih;

        % Shift x and w relative to y, according to horizon h
        the_ylead = lagmatrix(y, -h);
        Y(:,ih) = the_ylead;
        Xb(:,:,ih) = x*B(ih,:); % x*B
        
        if nargout>3

            % Residualize y and x*B(h) on controls
            the_select = ~isnan(the_ylead);
            the_reg = w(the_select,:);
            the_outc = [the_ylead(the_select) Xb(the_select,:,ih)];
            the_beta = the_reg\the_outc;
            the_resw = the_outc - the_reg*the_beta;

            Y_resw(the_select,ih) = the_resw(:,1);
            Xb_resw(the_select,:,ih) = the_resw(:,2:end);
        
        end

    end

end



