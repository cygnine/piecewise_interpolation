% Debugging poly_interpolate

clear
cd ..

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
xmin = 0;
xmax = 2*pi;
N = 100;
k = 6;

x = sort((xmax-xmin)*rand([N,1])+xmin);
x_fine = linspace(xmin,xmax,30*N);

%%%%%%%%%%%% Test smooth function
if true
  f = @(x) sin(x);
  df = @(x) cos(x);

  fx = f(x);
  dfx = df(x);

  d = poly_interpolate(x,fx,x_fine,'k',k);

  fprintf('Error is %3.3e\n', norm(f(x_fine)-d));

end

%%%%%%%%%%%%% Test discontinuous function
if true
  f_piecewise = @(x) double(x<1) + double(x>3);
  df_piecewise = @(x) 0*x;

  fx = f_piecewise(x);
  dfx = df_piecewise(x);
  
  d_piecewise = poly_interpolate(x,fx,x_fine,'k',k);

  fprintf('Error is %3.3f\n', norm(f_piecewise(x_fine)-d_piecewise));
end

cd debug
