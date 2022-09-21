clc; clear all; close all;
nm=1e-9; ev=1.60217e-19; m0=9.10938e-31; hc=1.05457e-34;
a=0.565325*nm; g1=6.8; g2=1.9; g3=2.73; del=0.34*ev; N=200; 
x1=0:N; x2=fliplr(-N:0); x3=x1;
E1=zeros(N+1,6); E2=E1; E3=E1;

%% Band along [100] direction
for i=0:N
   k=[1 0 0]*(i/N)*(2*pi/a);
   Q=((hc^2)*g2/(2*m0))*(k(1)^2+k(2)^2-2*k(3)^2);
   R=((hc^2)/(2*m0))*(-sqrt(3)*g2*(k(1)^2-k(2)^2)+sqrt(-1)*2*sqrt(3)*g3*k(1)*k(2));
   S=((hc^2)*g3/m0)*sqrt(3)*(k(1)-sqrt(-1)*k(2))*k(3);
   P=((hc^2)*g1/(2*m0))*(k(1)^2+k(2)^2+k(3)^2);
   H=-1.*[P+Q -S R 0 -S/sqrt(2) sqrt(2)*R;
           -S' P-Q 0 R -sqrt(2)*Q sqrt(1.5)*S;
           R' 0 P-Q S sqrt(1.5)*S' sqrt(2)*Q;
           0 R' S' P+Q -sqrt(2)*R' -S'/sqrt(2);
       -S'/sqrt(2) -sqrt(2)*Q' sqrt(1.5)*S -sqrt(2)*R P+del 0;
       sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q' -S/sqrt(2) 0 P+del];
   [wave,energy]=eig(H);
   E1(i+1,:)=sort(real(diag(energy)));
end

%% Band along [111] direction
for i=0:N
   k=[1/2 1/2 1/2]*(i/N)*(2*pi/a);
   Q=((hc^2)*g2/(2*m0))*(k(1)^2+k(2)^2-2*k(3)^2);
   R=((hc^2)/(2*m0))*(-sqrt(3)*g2*(k(1)^2-k(2)^2)+sqrt(-1)*2*sqrt(3)*g3*k(1)*k(2));
   S=((hc^2)*g3/m0)*sqrt(3)*(k(1)-sqrt(-1)*k(2))*k(3);
   P=((hc^2)*g1/(2*m0))*(k(1)^2+k(2)^2+k(3)^2);
   H=-1.*[P+Q -S R 0 -S/sqrt(2) sqrt(2)*R;
           -S' P-Q 0 R -sqrt(2)*Q sqrt(1.5)*S;
           R' 0 P-Q S sqrt(1.5)*S' sqrt(2)*Q;
           0 R' S' P+Q -sqrt(2)*R' -S'/sqrt(2);
       -S'/sqrt(2) -sqrt(2)*Q' sqrt(1.5)*S -sqrt(2)*R P+del 0;
       sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q' -S/sqrt(2) 0 P+del];
   [wave,energy]=eig(H);
   E2(i+1,:)=sort(real(diag(energy)));
end

%% Band along [110] direction
for i=0:N
   k=[1 1 0]*(i/N)*(2*pi/a);
   Q=((hc^2)*g2/(2*m0))*(k(1)^2+k(2)^2-2*k(3)^2);
   R=((hc^2)/(2*m0))*(-sqrt(3)*g2*(k(1)^2-k(2)^2)+sqrt(-1)*2*sqrt(3)*g3*k(1)*k(2));
   S=((hc^2)*g3/m0)*sqrt(3)*(k(1)-sqrt(-1)*k(2))*k(3);
   P=((hc^2)*g1/(2*m0))*(k(1)^2+k(2)^2+k(3)^2);
   H=-1.*[P+Q -S R 0 -S/sqrt(2) sqrt(2)*R;
           -S' P-Q 0 R -sqrt(2)*Q sqrt(1.5)*S;
           R' 0 P-Q S sqrt(1.5)*S' sqrt(2)*Q;
           0 R' S' P+Q -sqrt(2)*R' -S'/sqrt(2);
       -S'/sqrt(2) -sqrt(2)*Q' sqrt(1.5)*S -sqrt(2)*R P+del 0;
       sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q' -S/sqrt(2) 0 P+del];
   [wave,energy]=eig(H);
   E3(i+1,:)=sort(real(diag(energy)));
end

%% Figure band along [100] and [111] direction
figure(1)
plot(x1.*(2*pi/a/N)/(2*pi/a),E1./ev,'b','linewidth',1.5); hold on;
plot(x2.*(pi/a/N)/(2*pi/a),E2./ev,'b','linewidth',1.5);
axis([min(x2).*(pi/a/N)/(2*pi/a) max(x1).*(2*pi/a/N)/(2*pi/a) -14 0])
xlabel('\bfK');ylabel('\bfEnergy (eV)'); set(gcf,'color','white');
set(gca,'Xtick',[-.5 0 1]);set(gca,'Ytick',[-12,-10,-8,-6,-4,-2,0]);
set(gca,'Xticklabels',{'_{(^{1}/_{2},^{1}/_{2},^{1}/_{2})^{2\pi}/_{a}}L [111] ','{\Gamma} [000]','X [100]_{(1,0,0)^{2\pi}/_{a}}'});
set(gca,'fontsize',12);grid on;

%% curve fitting is done with these variables
% z=E1(:,6)./ev;
% x=x1.*(2*pi/a/N);

%% Figure band along [110] direction
figure(3)
plot(x3.*(2*pi/a/N)/(2*pi/a),E3./ev,'b','linewidth',1.5); 
axis([0 max(x3).*(3/4*2*pi/a/N)/(2*pi/a) -10 0])
xlabel('\bfK');ylabel('\bfEnergy (eV)'); set(gcf,'color','white');
set(gca,'Xtick',[0 .75]);set(gca,'Ytick',[-8,-6,-4,-2,0]);
set(gca,'Xticklabels',{'{\Gamma} [000]','K [110] _{(^{3}/_{4},^{3}/_{4}, 0)^{2\pi}/_{a}}'});
set(gca,'fontsize',12);grid on;
%% values at  k(0.2,0,0)
E=['heavy hole at k(0.2,0,0) is ',num2str(E1(N/5+1,6)/ev),' ev'];disp(E);
E=['light hole at k(0.2,0,0) is ',num2str(E1(N/5+1,4)/ev),' ev'];disp(E);
E=['split off band at k(0.2,0,0) is ',num2str(E1(N/5+1,2)/ev),' ev'];disp(E);

%% values at  k(0.4,0,0)
E=['heavy hole at k(0.4,0,0) is ',num2str(E1(N/2.5+1,6)/ev),' ev'];disp(E);
E=['light hole at k(0.4,0,0) is ',num2str(E1(N/2.5+1,4)/ev),' ev'];disp(E);
E=['split off band at k(0.4,0,0) is ',num2str(E1(N/2.5+1,2)/ev),' ev'];disp(E);

%% values at k(0.2,0.2,0)
E=['heavy hole at k(0.2,0.2,0) is ',num2str(E3(N/5+1,6)/ev),' ev'];disp(E);
E=['light hole at k(0.2,0.2,0) is ',num2str(E3(N/5+1,4)/ev),' ev'];disp(E);
E=['split off band at k(0.2,0.2,0) is ',num2str(E3(N/5+1,2)/ev),' ev'];disp(E); disp('Generating graphs... may take upto 5 minutes');

%% Iso surface generating and plotting
n = 200;   domain = .05*2*pi/a;
[x,y,z] = meshgrid(-domain/2:domain/n:domain/2, -domain/2:domain/n:domain/2, -domain/2:domain/n:domain/2);
% A=g1;B=2*g2;C=12*(g3^2-g2^2)*0.2;
% hh=-hc^2/(2*m0)*(A.*(x.^2+y.^2+z.^2)-sqrt((B.*(x.^2+y.^2+z.^2)).^2+C^2.*(x.^2.*y.^2+y.^2.*z.^2+z.^2.*x.^2)));
% lh=-hc^2/(2*m0)*(A.*(x.^2+y.^2+z.^2)+sqrt((B.*(x.^2+y.^2+z.^2)).^2+C^2.*(x.^2.*y.^2+y.^2.*z.^2+z.^2.*x.^2)));
hh=-hc^2/(2*m0*.33)*(x.^2+y.^2+z.^2);
lh=-hc^2/(2*m0*.0967)*(x.^2+y.^2+z.^2);

figure(4)
isosurface(x,y,z,hh,-.005*ev);
axis( [ -domain/2 domain/2 -domain/2 domain/2 -domain/2 domain/2 ]);
xlabel('\bfK_{x}');ylabel('\bfK_{y}');zlabel('\bfK_{z}');
set(gcf,'color','white'); box on; grid on; title('Heavy hole iso surface');

figure(5)
isosurface(x,y,z,lh,-.005*ev);
axis( [ -domain/2 domain/2 -domain/2 domain/2 -domain/2 domain/2 ]);
xlabel('\bfK_{x}');ylabel('\bfK_{y}');zlabel('\bfK_{z}');
set(gcf,'color','white'); box on; grid on;title('Light hole iso surface');

%% Density of states Stochastic approach
E=linspace(-1.5,0,76); del_E=E(2)-E(1);
D_Ehh=zeros(1,length(E));D_Elh=zeros(1,length(E));D_Eso=zeros(1,length(E));
for m=1:length(E)
    shh=0;slh=0;sso=0;N_total=10000;
    for N=1:N_total
        k=[2 0 0];
        while sum(abs(k))>1.5
            k=2*rand(1,3)-1;
        end
        k=k*(2*pi/a);
        Q=((hc^2)*g2/(2*m0))*(k(1)^2+k(2)^2-2*k(3)^2);
        R=((hc^2)/(2*m0))*(-sqrt(3)*g2*(k(1)^2-k(2)^2)+sqrt(-1)*2*sqrt(3)*g3*k(1)*k(2));
        S=((hc^2)*g3/m0)*sqrt(3)*(k(1)-sqrt(-1)*k(2))*k(3);
        P=((hc^2)*g1/(2*m0))*(k(1)^2+k(2)^2+k(3)^2);
        H=-1.*[P+Q -S R 0 -S/sqrt(2) sqrt(2)*R;
                -S' P-Q 0 R -sqrt(2)*Q sqrt(1.5)*S;
                R' 0 P-Q S sqrt(1.5)*S' sqrt(2)*Q;
                0 R' S' P+Q -sqrt(2)*R' -S'/sqrt(2);
            -S'/sqrt(2) -sqrt(2)*Q' sqrt(1.5)*S -sqrt(2)*R P+del 0;
            sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q' -S/sqrt(2) 0 P+del];
        [wave,energy]=eig(H);
        energy=sort(real(diag(energy)))./ev;

        if abs(E(m)-energy(4))<del_E
            slh=slh+1;
        end
        if abs(E(m)-energy(6))<del_E
            shh=shh+1;
        end  
        if abs(E(m)-energy(2))<del_E
            sso=sso+1;
        end  
    end
    D_Ehh(m)=(4/(a^3*del_E*N_total))*shh*1e-6;
    D_Elh(m)=(4/(a^3*del_E*N_total))*slh*1e-6;
    D_Eso(m)=(4/(a^3*del_E*N_total))*sso*1e-6;
end
D_E=D_Ehh+D_Elh+D_Eso;

%% Figure DOS stochastic approach
figure(6)

subplot(2,1,1)
plot(E,D_Ehh,'--','linewidth',1.5);hold on;
ghh=5.5e21*sqrt(0-E);
plot(E,ghh,'linewidth',1.5);hold on;title('Heavy hole band DOS');
xlabel('\bfEnergy (eV)');ylabel('\bfev^{-1}cm^{-3}');
set(gca,'Xtick',[-1.5 -1.25 -1 -.75 -.5 -.25 0]);
set(gcf,'color','white'); box on; grid off;
legend('Original DOS','Fitted Curve');

subplot(2,1,2)
plot(E,D_Elh,'--','linewidth',1.5);hold on;
glh=1.5e21*sqrt(0-E);
plot(E,glh,'linewidth',1.5);hold on; title('Light hole band DOS');
xlabel('\bfEnergy (eV)');ylabel('\bfev^{-1}cm^{-3}');
set(gca,'Xtick',[-1.5 -1.25 -1 -.75 -.5 -.25 0]);
legend('Original DOS','Fitted Curve');

%% DOS deterministic approach
kx=linspace(-1,1,100); ky=kx;kz=kx;t=length(kx)^3;
D_Ehhd=zeros(1,length(E));D_Elhd=zeros(1,length(E));D_Esod=zeros(1,length(E));
for i=1:length(kx)
    for j=1:length(ky)
        for m=1:length(kz)
            k=[kx(i) ky(j) kz(m)];
            if sum(abs(k))>1.5
                t=t-1;
                continue
            end 
            k=(2*pi/a).*k;
            Q=((hc^2)*g2/(2*m0))*(k(1)^2+k(2)^2-2*k(3)^2);
            R=((hc^2)/(2*m0))*(-sqrt(3)*g2*(k(1)^2-k(2)^2)+sqrt(-1)*2*sqrt(3)*g3*k(1)*k(2));
            S=((hc^2)*g3/m0)*sqrt(3)*(k(1)-sqrt(-1)*k(2))*k(3);
            P=((hc^2)*g1/(2*m0))*(k(1)^2+k(2)^2+k(3)^2);
            H=-1.*[P+Q -S R 0 -S/sqrt(2) sqrt(2)*R;
                    -S' P-Q 0 R -sqrt(2)*Q sqrt(1.5)*S;
                    R' 0 P-Q S sqrt(1.5)*S' sqrt(2)*Q;
                    0 R' S' P+Q -sqrt(2)*R' -S'/sqrt(2);
                -S'/sqrt(2) -sqrt(2)*Q' sqrt(1.5)*S -sqrt(2)*R P+del 0;
                sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q' -S/sqrt(2) 0 P+del];
            [wave,energy]=eig(H);
            energy=sort(real(diag(energy)))./ev;
            for n=1:length(E)
                if abs(E(n)-energy(4))<del_E
                    D_Elhd(n)=D_Elhd(n)+1;
                end
                if abs(E(n)-energy(6))<del_E
                    D_Ehhd(n)=D_Ehhd(n)+1;
                end  
                if abs(E(n)-energy(2))<del_E
                    D_Esod(n)=D_Esod(n)+1;
                end  
            end
        
        end
    end
end
D_Ehhd=(4/(a^3*del_E*t))*1e-6.*D_Ehhd;
D_Elhd=(4/(a^3*del_E*t))*1e-6.*D_Elhd;
D_Esod=(4/(a^3*del_E*t))*1e-6.*D_Esod;
D_Ed=D_Ehhd+D_Elhd+D_Esod;

%% Figure DOS deterministic vs stochastic approach
figure(7)
plot(E,D_Ehhd,'-.','linewidth',2);hold on;
plot(E,D_Ehh,'linewidth',1.5);
title('Heavy hole band DOS');
xlabel('\bfEnergy (eV)');ylabel('\bfev^{-1}cm^{-3}');
set(gca,'Xtick',[-1.5 -1.25 -1 -.75 -.5 -.25 0]);
set(gcf,'color','white'); box on; grid off; set(gca,'fontsize',12);
legend('Deterministic','Stochastic');