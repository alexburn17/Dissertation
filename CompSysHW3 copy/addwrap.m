function index=addwrap(a,b,dim) %toroidal a+b with wrap, indexed 1..dim
index=rem(a-1+b+dim,dim)+1; %add or subtract around torus