% clear;
function geom_data_reader(file_string)
% close all;
% movietitle=input('Name the movie file' , 's');
% 2-D Case
% file_string = 'riemannh02p0nomuscl';
fid = fopen(strcat(file_string,'.dat'),'r');
InputText=textscan(fid,'%s %f',1);
nt = InputText{2};nt=8197;
HeaderMain=sprintf('nt = %s',num2str(nt));
% disp(HeaderMain);
t = zeros(nt,1); nx = zeros(nt,1); ny = zeros(nt,1);dt = zeros(nt,1);
% M = zeros(
n = 0;
tag = 0;
while(~feof(fid))
    n = n + 1
    % for n = 1:nt
    t(n)=n*dt(n);
    sprintf('t = %s', num2str(t(n))); % Display block number
    % Read in t, nx, ny from block header
    InputText=textscan(fid,'%s %f %s %f %s %f %s %f',1);
    t(n) = InputText{2}; nx(n) = InputText{4};
    
    ny(n) = InputText{6};dt(n) = InputText{8};
    
    if n ~= 1
        t(n)=t(n-1)+dt(n);
    else
        t(n)=dt(n);
    end
    NumRows=nx(n)*ny(n);
    x=zeros(ny(n),nx(n));
    y=zeros(ny(n),nx(n));
    u=zeros(ny(n),nx(n));
    v=zeros(ny(n),nx(n));
    rho=zeros(ny(n),nx(n));
    p=zeros(ny(n),nx(n));
    vmag =zeros(ny(n),nx(n));
    theta=zeros(ny(n),nx(n));
    a=zeros(ny(n),nx(n));
    b=zeros(ny(n),nx(n));
    l=zeros(ny(n),nx(n));
    m=zeros(ny(n),nx(n));
    h=zeros(ny(n),nx(n));
    delta=zeros(ny(n),nx(n));
    
    data = zeros(NumRows,11);
    
    InputText=textscan(fid,'%s',6);
    InputText=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
    data=cell2mat(InputText);
    for i = 1:nx(n)
        for j = 1:ny(n)
            x(j,i) = data(j+(i-1)*ny(n),1);
            y(j,i) = data(j+(i-1)*ny(n),2);
            u(j,i) = data(j+(i-1)*ny(n),3);
            v(j,i) = data(j+(i-1)*ny(n),4);
            rho(j,i) = data(j+(i-1)*ny(n),5);
            p(j,i) = data(j+(i-1)*ny(n),6);
            vmag(j,i) = sqrt(u(j,i)^2+v(j,i)^2);
            theta(j,i)= atan2(v(j,i),u(j,i));
            
            a(j,i) = data(j+(i-1)*ny(n),7 );
            b(j,i) = data(j+(i-1)*ny(n),8 );
            l(j,i) = data(j+(i-1)*ny(n),9 );
            m(j,i) = data(j+(i-1)*ny(n),10);
            delta(j,i) = a(j,i)*m(j,i) - b(j,i)*l(j,i);
            
            h(j,i) = data(j+(i-1)*ny(n),11);
            
        end
    end
    skip = 1;

        if size(x,2)==1
    %         plot(y,rho)
        elseif mod(n,skip) == 0
     if size(x,2)>1
           figure(1)
           colormap(jet(128))
           mach=vmag./sqrt(1.4*p./rho);
           e=0.5*rho.*vmag.^2+1.0/0.4*p;
           vec=.6:.2:1.8;
           surf(x,y,mach)
   %         contourf(x,y,mach,vec);
   %         caxis([.6,1.8])
           axis equal;
           xlim([0. 5]);ylim([-.5 1.0]);%zlim([-0.1 25]);
           colorbar%('YLim',[.8,3.4],'YTick',[.8,3.4]);
           view([0,0,1]);
           xlabel('x');ylabel('y')
           shading interp
           drawnow
           if n > 180
               pause
           end
%           set(gca,'CLim',[.8,3.4],'YTick',[-.5,0,.5,1.0])
%          pause
     end
    %         figure(2)
    %         [dtdx,dtdy]=gradient(theta);
    %         [dgdx,dgdy]=gradient(log(h.*vmag));
    %         eq=(l.^2+m.^2).*(a.*cos(theta)-b.*sin(theta)).*dgdx+(a.^2+b.^2).*(m.*cos(theta)-l.*sin(theta)).*dgdy...
    %             +(l.^2+m.^2).*(a.*cos(theta)+b.*sin(theta)).*dtdx-(a.^2+b.^2).*(l.*cos(theta)+m.*sin(theta)).*dtdy;
    %         contourf(x,y,eq);
    %         colorbar;
    
%             pause(.01)
    
    %         pause
    %         figure(2)
    %         surf(x,y,u);axis equal;xlim([0 0.55]);ylim([-0.5 0.5]);
    %         figure(2)
    %         contour(x,y,u(:,:));axis equal;%xlim([0 0.5]);ylim([0,1]);%zlim([-0.1 2.1]);colorbar;
    %         fprintf('Step')
    %         title('Mach Number')
    %         if(n>1)
    %             if tag == 0;
    %                 pause;tag = 1;
    %             end
    %             M(n-1) = getframe(gcf);
    %         end
    
    
        end
    end
    %                 figure(3)
    %                 plot(y(:,50)./x(:,50),rho(:,50),'.');ylim([.4 1.]);xlim([-0.6 0.4 ]);
    %         ylabel('\rho');xlabel('^y/_x')
    %         title('Density for riemann test problem')
    % movie2avi(M,movietitle);
    
    % y0=linspace(0,1,1000);
    % X=reshape(x,size(x,1)*size(x,2),1);
    
    % Y=reshape(y,size(x,1)*size(x,2),1);
    % RHO=reshape(rho,size(x,1)*size(x,2),1);
    % dens_interp=TriScatteredInterp(X,Y,RHO);
    % for n = 1:length(y0)
    %     dens(n)=dens_interp(.25,y0(n));
% end
% figure(4)
% plot(y0./.25,dens);
if file_string(1:5) == 'riema'
    statesin=[.25,.5,7*sqrt(1.4*.25/.5),0;1,1,2.4*sqrt(1.4*1/1),0]';
    [steady_riemann_ans,sampler]=exact_riemann_steady_main;
    [exact_states,exact_angles] = steady_riemann_ans(statesin);
    half_first_column_width = -(x(:,2)-x(:,1));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            exact_solution(i,j,:) = sampler((y(i,j)-.5)./(x(i,j)-half_first_column_width(i)),exact_states,exact_angles);
        end
    end
elseif file_string(1:5) == 'shock'
    beta = 45/180*pi;
    M1 = 1.8;
    d1 = 1.0;
    p1 = 1.0;
    gamma = 1.4;
    a1 = sqrt(gamma*p1/d1);
    Mn1=M1*sin(beta);
    rhs = 2*cot(beta)*((M1^2*sin(beta)^2-1)/(M1^2*(gamma+cos(2*beta))+2));
    d2 = d1*(gamma+1)*Mn1^2/((gamma-1)*Mn1^2+2);
    p2 = p1*(1+2*gamma/(gamma+1)*(Mn1^2-1));
    a2 = sqrt(gamma*p2/d2);
    Mn2= (Mn1^2+(2/(gamma-1)))/((2*gamma/(gamma-1))*Mn1^2-1);
    M2 = Mn2/sin(beta-atan(rhs));
    %     atan(rhs)
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            if y(i,j) < 1-(x(i,j)-.1)*tan(beta)
                exact_solution(i,j,:) = [p1,d1,a1*M1,0];
            else
                exact_solution(i,j,:) = [p2,d2,a2*M2*cos(atan(rhs)),-a2*M2*sin(atan(rhs))];
            end
        end
    end
elseif file_string(1:5) == 'expan'
    [steady_riemann_ans,sampler]=exact_riemann_steady_main;
    gamma = 1.4;
    M1 = 2/sqrt(3);
    P1 = 1.0;
    D1 = 1.0;
    Theta1 = 0.0;
    a1 = sqrt(gamma*P1/D1);
    M2 = 2;
    P2 = P1 / ( (1+(gamma-1)/2*M2^2) / (1+(gamma-1)/2*M1^2) )^(gamma/(gamma-1));
    D2 = D1 * (P2/P1)^(1/gamma);
    a2 = sqrt(gamma*P2/D2);
    Theta2 = PMfunction(M2)-PMfunction(M1);
    state1 = [1.0,1.0,M1*a1,0.0];
    state2 = [P2, D2, M2*a2*cos(Theta2), M2*a2*sin(Theta2)];
    states = zeros(4,4);
    states(:,1) = state2;
    states(:,2) = state2;
    states(:,3) = state2;
    states(:,4) = state1;
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            if x(i,j) > .1
                exact_solution(i,j,:) = sampler((y(i,j)-1.)/(x(i,j)-.1),states,[1,1,1,Theta2-asin(1/M2),Theta1-asin(1/M1)]);
            else
                exact_solution(i,j,:) = state1;
            end
        end
    end
end
% error=(cat(3,p,rho,u,v)-exact_solution)./exact_solution;
% for indx = 1:4
%     error(:,:,indx)=error(:,:,indx).*(a.*m-b.*l);
% end
% rmserror=reshape(sqrt(sum(sum((error).^2,1),2)/(size(error,1)*size(error,2))),4,1);
% save(strcat(file_string,'error.mat'),'error','rmserror')
save(strcat(file_string,'data.mat'))

% surf(x,y,sqrt((rho-steady_riemann_ans((y-.5)/x)).^2))
%surf(x,y,rho);view([0,0,1]);
end
function out = PMfunction(M)
gamma = 1.4;
out = sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M^2-1)))-atan(sqrt(M^2-1));
end
