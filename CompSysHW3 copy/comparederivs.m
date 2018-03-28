function comparederivs(allx,f,truedf)
% function comparederivs(allx,f,truedf)
% Compares first and second order approximations to derivatives as a
% function of step size h
%
% allx is a scalar or vector of what x-values to look at
%   This program produces one plot for each x value
% f is function handle of the function to differentiate
% truedf is a function handle for the analytic derivative of f
%
% sample call: comparederivs([0 5],@exp, @exp)

% Maggie Eppstein, 2/10/08

h=logspace(-20,-1,20); %step sizes to try

for xi=1:length(allx)
    figure
    x=allx(xi);
    
    % compute derivative and approximations
    df=truedf(x); %true derivative at x
    df1=backdiff(x,f,h); %approximate with 1st order backwards difference
    df2=centraldiff(x,f,h); %approximate with 2nd order central difference
    
    % compute the absolute errors for all h
    err1=abs(df-df1);
    err2=abs(df-df2);

    % plot the errors as a function of h
    loglog(h,err1,'ro',h,err2,'b*');
    
    % linear regression of log-log relationships (using last 5 points only)
    coef1=polyfit(log(h(end-4:end)),log(err1(end-4:end)),1);
    coef2=polyfit(log(h(end-4:end)),log(err2(end-4:end)),1);
    
    % plot the best-fit lines
    hold on
    loglog(h(end-4:end),exp(coef1(2))*h(end-4:end).^coef1(1),'r-');
    loglog(h(end-4:end),exp(coef2(2))*h(end-4:end).^coef2(1),'b-');
    
    % label the plots
    set(gca,'fontsize',14)
    xlabel('h')
    ylabel('error')
    legend('back diff','central diff',...
        ['slope = ',num2str(coef1(1),4)],['slope = ',num2str(coef2(1),4)],...
        'Location','BestOutside');
    title(['at x = ',num2str(x)])
    
    if xi<length(allx)
        pause
    end
end
figure(gcf)
    
% NOTE: the following functions are place here for convenience for this demo code, 
% but cannot be called from outside this file.

function df1=backdiff(x,f,h)
df1=(f(x)-f(x-h))./h;

function df2=centraldiff(x,f,h)
df2=(f(x+h)-f(x-h))./(2*h);