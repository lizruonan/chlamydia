%% Logistic Growth Model for  _C.Trachomatis_
% Authors: Elizabeth Zhao, Liam Vu
% Date: July, 2015
function [ERR,umax,t0] = logistic(x) 
%The input of the function logistic is x that stands for the number of entries 
%in the interval of approximate u_max and t0. This script document is
%designed to solve for ODEs of a logistic growth model we have and get
%the optimal value of u_max and t0 by the method of least squares. Then we
%plot the logistic model and the errors to help us comprehend the
%system. 
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 

%First of all, we set up all sorts of initial conditions for the inputs.
%init is the intial condition R_0=1 and E_0=0 at time t=12. "tspan" is a
%row vector of time t from 12 to 40. By analysing the datas, we got
%approximatae intervals for u_max and t0. "m" is a row vector of u_max 
%from 0.0 to 0.5 with x entries. "t0_vec" is also a row vector between 16
%and 32 with x entries. K is the carrying capacity where we got it from
%looking at datas of the largest average value of RB cells;
%alpha is the growth rate where we got it by solving the ODEs for R before the switch
%point t0 and u(t) is 0. 
init = [1,0]; % R_0
tspan = [12,16,20,24,28,32,36,40];
m = linspace(0.0,0.5,x); %u_max
t0_vec = linspace(16,32,x)';
K=260.1;
alpha=0.476;

%Thanks to the collaborators, we got the datas of time T, RB and EB to evaluate the
%least squares. Tdata, Rdata, and Edata are three column vectors with the
%same amount of entries.
Tdata = [12,16,20,24,28,32,36,40]';
Rdata = [0.58,3.90,16.82,43.33,195.46,260.1,151.8,129.5]';
Edata = [0, 0, 0, 0.33, 85, 384, 456.5, 547.75]';

% This nested for loop is to iterate every element in the vactors u_max and t0 to the
% function fun and myeror. We get the results of RB, EB and a matrix of
% errors of RB and EB. 
for i=1:x
    for j=1:x
        umax=m(i);
        t0=t0_vec(j);
        [tt,y] =ode45(@fun, tspan,init);
        Rf = y(:,1);
        Ef = y(:,2);
        ERR(i,j)=myerror(Rf,Ef,Rdata,Edata);
               
    end
end


% We initially set the maximum error value from the matrix ERR.
% Then in a nested for loop, we use an if statement to check
% if there are any smaller value of ERR and update it to teh
% maximum error. If there are no more values that smaller
% than the updated maximum, then it is the minimum error.
% This gives us the value and the position of the minimum
% error.
Min=max(ERR);
for i=1:x
    for j=1:j
        if ERR(i,j)<Min
            imin=i; jmin=j;
            Min=ERR(i,j);
        end
    end
end

imin
jmin
Min


% We want to keep track of the position of 
% the optimal t0 and u_max by the hint from
% the least squares. 
umax=m(imin);  t0=t0_vec(jmin); 


% We want to solve for the ODEs of dRdt and
% dEdt, and we already know the optimal
% t0 and u_max by least squares. "tt" is
% a column vector of all the points of
% tspan, and y is a matrix where the
% first column is RB and the second
% column is EB. 
[tt,y] =ode45(@fun, tspan,init); 
Rf = y(:,1);   
Ef = y(:,2);

% We plot the resulting logistic growth model, where the 
% optimal value of u_max and t0 is given by least squares. 
% We also plot the data to see how the values fit the data. 
%subplot(1,3,1);
figure;
plot(Tdata,Rdata,'go',Tdata, Edata, 'ro',tt, Rf,'g', tt, Ef,'r', 'LineWidth', 2);
grid on
axis square
xlabel('Time (hours)') 
ylabel('The Amount of RB cells and EB cells')
title('The Logistic Model of RB cells and EB cells')


% We want to see what the error value looks like
% on the surface in the domain. The meshgrid
% creates matrices of t0 and u_max that contains
% an nXn matrix. This helps us to create a 3-D
% diagram
% By the values of Error we got, it is not hard 
% to plot a 3-D diagram of error of u_max versus t0
% by "surf". Then we can get the approximate
% minimum error on the lowest point of the surface. 
%subplot(1,3,2);
figure;
[X,Y]=meshgrid(t0_vec,m);                
surf(X,Y,ERR)  ;                                   
ylabel('u_m_a_x');
xlabel('t_0');
zlabel('error');
title('Error versus t_0 and u_m_a_x')


%We can also get an indication of where the 
%minimum lies on the surface by creating a contour
%plot. So the minimum of t0 and u_max 
%approximatly lies in the blank circle on the map.  
%subplot(1,3,3);
figure;
contour(t0_vec,m,ERR,40)
hold on
plot(23.11, 0.22,'b.','MarkerSize',10)
xlabel('t_0')
ylabel('u_m_a_x')
title('A contour map of the errors versus u_m_a_x and t_0 ')

    function  res1=fun(t,v) 
    %The function has inputs t and v, where the variable v contains
    %the variables R and E. By calling the function
    %myerror, we get the value of u depending on t,
    %then we plug in the values alpha, K, and u in
    %to our logistic growth model. We will call this
    %function at the start to run the function
    %ode45 to solve for R and E.
    %This function returns a column vector of  
    %dRdt and dEdt.
    R = v(1); 
    E = v(2);
    u = ufun(t);
 
    dRdt = alpha*R*(1-R/K)-u*R; %Carrying Capacity K=260.1; a =0.15
    dEdt = u*R;
    res1 = [dRdt; dEdt];
    end

    function res=myerror(R,E,datax,datay)   
    %The inputs are R, E, datax, and datay, where R and E are
    %results of the ODEs R and E. Datax and Datay are the
    %datas of RB and EB cells. We want to the error at each point
    %of the tspan. We intent to call the function in a
    %nested for loop in the beginning to find the
    %minimal error value. The  evolution of  EB and RB
    %cells are a system, so the function returns
    %a result of the squre root of the sum of
    %square of error R and error E.                                      
    errR = R - datax; %x = Rf, datax = Rdata
    errE = E - datay; %y = Ef, datay = Edata       
    res=sqrt(sum(errR.^2)+ sum(errE.^2));
        
    end

    function res = ufun(t)
    %This function returns the value of u. 
    %We assume that the system is bang-bang
    %control.
    %If t is before the switch point t0,
    %then u is 0. Otherwise, u is u_max. 
    res = 0;
        if t>t0
            res= umax;  
        end
                                     
    end %of ufun

end %of main function





