function obj = locproj_cv(obj,K,lambda)

  L = length(lambda);  
  rss = zeros(L,1);
  aic = zeros(L,1);
  
  for l = 1:L       
      fprintf('.')
        S = obj.X * inv( obj.X'*obj.X + lambda(l) * obj.P ) * obj.X';
        rss(l) = sum( ( (obj.Y - S * obj.Y) ./ ( 1 - diag(S) ) ).^2 );
  end
  fprintf('.')
  
  obj.rss     = rss;
  obj.lambda  = lambda;

end