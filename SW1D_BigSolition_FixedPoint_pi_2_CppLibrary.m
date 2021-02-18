%  function SW1D_BigSolition_FixedPoint_pi_2_CppLibrary(n)
%          % Try-catch expression that quits the Matlab session if your code crashes
%          try
%                  % Initialize the parallel pool
%                  c=parcluster();
%                  t=tempname();
%                  mkdir(t)
%                  c.JobStorageLocation=t;
%                  c.NumWorkers=8;
%                  parpool(c,n);
%                  % The actual program calls
%                  BFGS_para1D_SW_00000065;
%                  delete(gcp('nocreate'));
%          catch error
%                  getReport(error)
%                  disp('Error occured');
%                  exit(0)
%          end
%  end
% function SW1D_BigSolition_FixedPoint_pi_2_CppLibrary(n) should be
% uncommented if running as series script on computation cluster by shell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BFGS_para1D_SW_00000065 
format long
addpath('ClassCalculateA20InOneElement'); % add library to path 

xidt=1;xid=sqrt(2)*xidt;
Deltap=1;
Deltav1=0.1*Deltap;
gamma=Deltav1/Deltap;
gamma1=(1+gamma^2)/((1+gamma)^2);
gamma2=(gamma^2)/((1+gamma)^2);
G3=gamma/(1+gamma);

% Nelem=500; %500 points hard to arrive 10-6 tolerence
Nelem=100;
Np=Nelem+1;

xmin=-10*xidt;xmax=10*xidt;
x=xmin:((xmax-xmin)/Nelem):xmax;


A=sparse(Np,Np);

thetakn1t=sparse(ones(Np-2,1)*1.0);

thetak1=asin(-G3/2);thetakN1=pi-asin(-G3/2);
thetakn1=[thetak1; thetakn1t; thetakN1;];

Is=sparse(eye(Nelem-1));
Hkn1=Is;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=6.5e-6;
%option=optimset('TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',2000); This option
%will prevent iteration arrive tol 6e-6
gkn1=1;
D=[];
o=1;
O=[];
c=0;
C=[];
Lambda=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make the A matrix %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:Nelem
     
       clear functionA20;
     % functionA20 object
       functionA20=clib.ClassCalculateA20InOneElement.CalculateA20;

     % call the member function to calculate
       functionA20.CalculateA20InOneElement(x(i),x(i+1));

     % return the symmtrilized result
       A20=reshape(functionA20.returnResult(4),2,2);
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   A(i,i)=A20(1,1)+A(i,i);
   A(i,i+1)=A20(1,2)+A(i,i+1);
   A(i+1,i)=A20(2,1)+A(i+1,i);
   A(i+1,i+1)=A20(2,2)+A(i+1,i+1);

end

 
gkn1=gradient1DQN(A,thetakn1,Np,Nelem,x,gamma,gamma1,gamma2,xid);     
while norm(gkn1)>=tol 
      c=c+1;
    
      if norm(gkn1)<=5e-5 
         save('1DSW_00000065_xidt_MyEq_gamma0005_5e_5.mat','x','thetakn1','C','O','Lambda','gkn1','Hkn1','tol');    
      end
     
          C(c)=c;
          clc
          c
  
             O(c)=norm(gkn1);
             ABSgkn1=norm(gkn1)
             
             if norm(gkn1)<tol
                display('break:gkn1<tol'); 
                break;
                 
             else   
                    dkn1=-Hkn1*gkn1;   
                    
                    
                    [lambdakn1,fval,exitflag,output]=fminsearch(@(lambda) PithetakLambdaQN(lambda,A,thetakn1,dkn1,Nelem,x,gamma,gamma1,gamma2,xid),0);
                    %[lambdakn1,fval,exitflag,output]=fminbnd(@(lambda) PiukLambdaQNSinAL(lambda,A,ukn1,dkn1,Nelem,x),0,10);
       
         
                    thetakt=thetakn1(2:Nelem)+lambdakn1.*dkn1;
                    thetak=[thetak1; thetakt; thetakN1;];
             
                    sk=thetak(2:Nelem,:)-thetakn1(2:Nelem,:);
                    
                    gk=gradient1DQN(A,thetak,Np,Nelem,x,gamma,gamma1,gamma2,xid); 
                    yk=gk-gkn1;
                    rhok=1./((yk')*sk);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
                    Hk=InHassian(Is,rhok,sk,yk,Hkn1);
                    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
         
                    Hkn1=Hk;
                    thetakn1=thetak;
                    gkn1=gradient1DQN(A,thetakn1,Np,Nelem,x,gamma,gamma1,gamma2,xid); 
        
         
         
                    Lambda(c)=lambdakn1;
                    LAMBDA=lambdakn1
             
             end

    
     
     
     
end     
% save('1DSW_00000065_MyEq_100_gamma0005.mat','x','thetakn1','C','O','Lambda','gkn1','Hkn1','tol');
end

%%%%%%%%%%%%%%%%%%%%%%% PiukLambda %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PithetakL=PithetakLambdaQN(lambda,A,Thetak,d,Nelem,x,gamma,gamma1,gamma2,xid)


thetak1=Thetak(1);
thetakN1=Thetak(length(Thetak));

thetakt=Thetak(2:Nelem)+lambda.*d;

Thetak=[thetak1; thetakt; thetakN1;];

Pithetak1=((Thetak')*A)*Thetak;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pithetak3=zeros(1,Nelem);
parfor k=1:Nelem
    
       Pithetak3(k)=IntegralPithetakLambda(x,k,Thetak,gamma,gamma1,gamma2,xid);
end  
   

PithetakL=Pithetak1+sum(Pithetak3);
end

function Pithetak3k=IntegralPithetakLambda(x,k,Thetak,gamma,gamma1,gamma2,xid)
gamma3=gamma3Dv2(x(k+1),gamma); 
Pithetak3k=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) -gamma3*sin(1.*((x(k+1)-q)./(x(k+1)-x(k))).*Thetak(k)+1.*((q-x(k))./(x(k+1)-x(k))).*Thetak(k+1)),x(k),x(k+1),'RelTol',1e-9)+...
            (1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) -(1/2)*cos(2.*((x(k+1)-q)./(x(k+1)-x(k))).*Thetak(k)+2.*((q-x(k))./(x(k+1)-x(k))).*Thetak(k+1)),x(k),x(k+1),'RelTol',1e-9);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dPidk=gradient1DQN(A,thetak,Np,Nelem,x,gamma,gamma1,gamma2,xid)

     dPidkA=(2*A(2:Np-1,:))*thetak;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     dPidkB=zeros(Nelem-1,1);
     parfor k=2:Nelem
    
          dPidkB0=IntegralGrad(x,thetak,k,gamma,gamma1,gamma2,xid);
          dPidkB(k-1,1)=dPidkB0;
    
     end    
     dPidk=dPidkA+dPidkB; 
end

function dPidkB0=IntegralGrad(x,thetak,k,gamma,gamma1,gamma2,xid)
         
         if x(k+1)<=0
           gamma3=gamma3Dv2(x(k+1),gamma);
          
           dPidk1=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) (-gamma3*cos(1.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+1.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k))+...
                  1*sin(2.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+2.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k))).*((q-x(k-1))./(x(k)-x(k-1))),x(k-1),x(k),'RelTol',1e-9);
           
           dPidk2=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) (-gamma3*cos(1.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+1.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1))+...
                  1*sin(2.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+2.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1))).*((x(k+1)-q)./(x(k+1)-x(k))),x(k),x(k+1),'RelTol',1e-9);
              
        elseif x(k-1)<0 && x(k+1)>0 
         
               gamma3_1=gamma3Dv2(x(k-1),gamma);
               dPidk1=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) (-gamma3_1*cos(1.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+1.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k))+...
                      1*sin(2.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+2.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k))).*((q-x(k-1))./(x(k)-x(k-1))),x(k-1),x(k),'RelTol',1e-9);
        
               gamma3_2=gamma3Dv2(x(k+1),gamma);
               dPidk2=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) (-gamma3_2*cos(1.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+1.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1))+...
                      1*sin(2.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+2.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1))).*((x(k+1)-q)./(x(k+1)-x(k))),x(k),x(k+1),'RelTol',1e-9);
    
        elseif x(k-1)>=0 && x(k)>0
     
               gamma3=gamma3Dv2(x(k),gamma);
               dPidk1=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) (-gamma3*cos(1.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+1.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k))+...
                  1*sin(2.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+2.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k))).*((q-x(k-1))./(x(k)-x(k-1))),x(k-1),x(k),'RelTol',1e-9);
           
               dPidk2=(1/(xid*xid)*(gamma1+2*gamma2))*integral(@(q) (-gamma3*cos(1.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+1.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1))+...
                  1*sin(2.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+2.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1))).*((x(k+1)-q)./(x(k+1)-x(k))),x(k),x(k+1),'RelTol',1e-9);
        end  


%           dPidk1=(1/(2))*integral(@(q) sin(2.*((x(k)-q)./(x(k)-x(k-1))).*thetak(k-1)+2.*((q-x(k-1))./(x(k)-x(k-1))).*thetak(k)).*1.*((q-x(k-1))./(x(k)-x(k-1))),x(k-1),x(k),'RelTol',1e-9);
%           dPidk2=(1/(2))*integral(@(q) sin(2.*((x(k+1)-q)./(x(k+1)-x(k))).*thetak(k)+2.*((q-x(k))./(x(k+1)-x(k))).*thetak(k+1)).*1.*((x(k+1)-q)./(x(k+1)-x(k))),x(k),x(k+1),'RelTol',1e-9);

          dPidkB0=dPidk1+dPidk2;
end

%%%%%%%%%%%%%%%%%%%%%%%% Inverse Hassian %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Hk=InHassian(Is,rhok,sk,yk,Hkn1)

left=Is-rhok*(sk*(yk'));
right=Is-rhok*(yk*(sk'));

InHassian1=(left*Hkn1)*right;
Inhassian2=rhok*(sk*(sk'));

Hk=InHassian1+Inhassian2;
end

%%%%%%%%%%%%%%%%%%%%% gamma3 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma3=gamma3Dv2(X,gamma)
  

if X<=0
    omega=-gamma;
else
    omega=gamma;
end    
    gamma3=omega/(1-gamma);
end