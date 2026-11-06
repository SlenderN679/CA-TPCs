function out=aula1_demonstracao_observador(showPlots,stage)
% Demonstracao da aula 1: estimadores em malha aberta e fechada.
% A planta recebe u exogeno; a malha fechada refere-se ao OBSERVADOR.
% Requer Control System Toolbox. Executar: aula1_demonstracao_observador
if nargin<1,showPlots=true;end
if nargin<2,stage='comparar';end
assert(any(strcmp(stage,{'aberta','fechada','comparar'})),'Usar aberta, fechada ou comparar.');
% Ex.: aula1_demonstracao_observador(true,'aberta')
A=[-1 1.5;1 -2]; B=[1;0]; C=[1 0]; D=0;
sys=ss(A,B,C,D); %#ok<NASGU>
assert(rank(obsv(A,C))==2);
L=place(A',C',[-3 -4])';
x0=[-.5;-1]; xh0=[0;0];
t=(0:.1:10)';
u=[ones(15,1);-ones(15,1);.5*ones(15,1);-.5*ones(15,1);zeros(41,1)];
assert(numel(u)==numel(t));
Aol=blkdiag(A,A); Acl=[A zeros(2);L*C A-L*C];
Baug=[B;B]; Caug=eye(4); Daug=zeros(4,1);
% ZOH: entrada por patamares, igual nos dois estimadores.
% lsim: manter u constante entre amostras (ZOH).
[Zol,~,~]=lsim(ss(Aol,Baug,Caug,Daug),u,t,[x0;xh0],'zoh');
[Zcl,~,~]=lsim(ss(Acl,Baug,Caug,Daug),u,t,[x0;xh0],'zoh');
Eol=Zol(:,1:2)-Zol(:,3:4); Ecl=Zcl(:,1:2)-Zcl(:,3:4);
% Verificacao independente: e(t)=expm((A-LC)t)e(0).
exact=zeros(size(Ecl));
for j=1:numel(t),exact(j,:)=(expm((A-L*C)*t(j))*(x0-xh0))';end
assert(max(abs(Ecl-exact),[],'all')<1e-9);
fprintf('eig(A):\n');disp(eig(A));fprintf('L:\n');disp(L);
fprintf('eig(A-LC):\n');disp(eig(A-L*C));
fprintf('Erro final OL=%.4g; CL=%.4g\n',norm(Eol(end,:)),norm(Ecl(end,:)));
if showPlots
 modes=1:2;
 if strcmp(stage,'aberta'),modes=1;elseif strcmp(stage,'fechada'),modes=2;end
 for mode=modes
  if mode==1,Z=Zol;E=Eol;lab='Malha aberta';else,Z=Zcl;E=Ecl;lab='Malha fechada do observador';end
  figure('Name',lab,'Color','w');
  subplot(2,1,1);plot(t,Z(:,1:2),'-',t,Z(:,3:4),'--','LineWidth',1.4);
  grid on;title(lab);ylabel('estados');legend('x1','x2','xhat1','xhat2');
  subplot(2,1,2);plot(t,E,'LineWidth',1.4);grid on;xlabel('tempo [s]');ylabel('erro');legend('e1','e2');
 end
 if strcmp(stage,'comparar')
 figure('Name','Comparacao do erro','Color','w');
 semilogy(t,max(vecnorm(Eol,2,2),eps),t,max(vecnorm(Ecl,2,2),eps),'LineWidth',1.4);
 grid on;xlabel('tempo [s]');ylabel('norma do erro');legend('malha aberta','malha fechada');
 end
end
out=struct('t',t,'u',u,'Zol',Zol,'Zcl',Zcl,'Eol',Eol,'Ecl',Ecl,'L',L,'exact',exact);
end

