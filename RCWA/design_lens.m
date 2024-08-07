%script analyzes Kerolos library of circular pillars
clear all
folder=['C:\Users\Maryna\Dropbox (Harvard University)\Project_GasCamera\Design\Test_design\Circular Pillars Lib - Kerolos'];
ref=['REF_phase.mat'];
reff=fullfile(folder,ref);
load(reff);
refEx=REF_Ex;
for i=1:111
    filename=['phases_',num2str(i),'.mat'];
    f=fullfile(folder,filename);
    load(f);
    phase(:,i)=angle(Ex./refEx);
    ampl(:,i)=abs(Ex./refEx).^2;
end
%filter the values of amplitude: two conditions should be satisfied:
%1)amplitude is higher than 1.01
%2)amplitude is lower than 0.8
%for practice choose single wavelength lambda=540 nm element 15
ampl_f=ampl(15,:);
phase_f=phase(15,:);
[p,s]=size(ampl_f);
struc=1:s;%this variable connects name of the structure to the design
k=1;
for i=1:s
    if ampl_f(i)>1.01 || ampl_f(i)<0.8
        num(k)=i;
        k=k+1;
    end
end

[p,s]=size(num);
k=k-1;
for i=1:s
    ampl_f(num(k))=[];
    phase_f(num(k))=[];
    struc(num(k))=[];
    k=k-1;
end
  
X=linspace(-0.5*10^(-3),0.5*10^(-3),2500);
Y=linspace(-0.5*10^(-3),0.5*10^(-3),2500);

[p,s]=size(X);

phase=ones(s);
lambda_w=lambda(15);
f=5*10^(-3);
for i=1:s
    for j=1:s
        phase(i,j)=-((2*pi)/lambda_w)*(sqrt(X(i)^2+Y(j)^2+f^2)-f);
    end
end

phase_req = wrapTo2Pi(phase);%required phase
phase_lib = wrapTo2Pi(phase_f);%library phase

%optimize the design
for i=1:s
    for j=1:s
        [ttt,design(i,j)]=min(abs(exp(1i*phase_req(i,j))*ones(1,73)-exp(1i*phase_lib)));
        phase_design(i,j)=phase_lib(design(i,j));
    end
end

