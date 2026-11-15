function out=aula2_ganho_L_students(showPlots,stage)
% Demonstracao da aula 2: como se calcula o ganho L do observador.
% Exemplo dos slides: A=[1 2 0;3 4 1;1 2 0], C=[0 1 1], polos {-10,-2+-2j}.
% Processo I = forma canonica; Processo II = dualidade (coeficientes ou Ackermann).
% Requer Control System Toolbox. MATLAB R2020b+.
%
%   aula2_ganho_L_students           -> tudo (calculo, simulacao, discreto)
%   aula2_ganho_L_students(true,'calculo')    -> L pelos 4 processos, lado a lado
%   aula2_ganho_L_students(true,'simular')    -> planta INSTAVEL: malha aberta vs fechada
%   aula2_ganho_L_students(true,'discreto')   -> mesmo projeto em z (liga ao hands-on)
%   r=aula2_ganho_L_students(false);          -> so calcula e verifica, sem janelas
%
% A planta NAO tem controlador: u e um sinal conhecido. A "malha fechada"
% refere-se ao OBSERVADOR (correcao pela inovacao y - C*xhat), como na aula 1.
if nargin<1, showPlots=true; end
if nargin<2, stage='tudo'; end
assert(any(strcmp(stage,{'calculo','simular','discreto','tudo'})), ...
    'Usar calculo, simular, discreto ou tudo.');

A=[1 2 0;3 4 1;1 2 0]; B=[1;-1;0]; C=[0 1 1]; D=0;
polos=[-10 -2+2i -2-2i];
n=size(A,1);
fprintf('eig(A) = %s\n',mat2str(eig(A).',4));
fprintf('-> a planta tem um polo em +5.70: e INSTAVEL.\n\n');

%% 0) Observabilidade
O=obsv(A,C);
assert(rank(O)==n,'Par (A,C) nao observavel: nao se podem colocar todos os polos.');

%% 1) PROCESSO I - forma canonica observavel
pd=real(poly(polos));        % [1 a_{n-1} ... a_0]  (desejado)
pc=poly(A);                  % [1 c_{n-1} ... c_0]  (malha aberta)
% TODO 1: ganho na forma canonica, Lo = pd - pc (vetor coluna, sem o 1 inicial).
Lo=nan(n,1);
Ao=[-pc(2:end).' [eye(n-1);zeros(1,n-1)]];   % forma canonica observavel
Co=[1 zeros(1,n-1)];
Oo=obsv(Ao,Co);
% TODO 2: matriz de transformacao T = O^-1*Oo e ganho L1 = T*Lo.
T=nan(n); L1=nan(n,1);

%% 2) PROCESSO II (a) - igualar coeficientes de det(sI-A+LC) a pd(s)
% det(sI-A+LC) e afim em L: coef(L) = coef(L=0) + M*L. A coluna i de M
% obtem-se com L = e_i. Resolver M*L = pd - pc e o sistema de equacoes do slide.
[L2,M]=ganhoPorCoeficientes(A,C,polos);

%% 3) PROCESSO II (b) - formula de Ackermann: L = pd(A) * O^-1 * [0;...;0;1]
% TODO 3: L3 = pd(A)*O^-1*[0;...;0;1]. Dica: polyvalm(pd,A).
L3=nan(n,1);

%% 4) MATLAB: dualidade (A',C') -> L'
L4=place(A',C',polos).';     % .' transposta simples (nao conjugada)
L5=acker(A',C',polos).';

Ls=[L1 L2 L3 L4 L5];
ok=~any(isnan(Ls),1);
assert(max(abs(Ls(:,ok)-L4),[],'all')<1e-6*max(1,norm(L4)),'Os processos nao coincidem: rever os TODOs.');
if ~all(ok), fprintf('Ainda por completar: TODO(s) com NaN na tabela.\n'); end
L=L4;
if any(strcmp(stage,{'calculo','tudo'}))
    fprintf('pd(s) coef.: %s   pc(s) coef.: %s\n',mat2str(pd,4),mat2str(pc,4));
    fprintf('Lo = pd - pc = %s\n',mat2str(Lo.',4));
    fprintf('M (Processo II) =\n'); disp(round(M,10));
    disp(array2table(Ls,'VariableNames',{'ProcI_canonica','ProcII_coef', ...
        'ProcII_Ackermann','place','acker'},'RowNames',{'l1','l2','l3'}));
    fprintf('eig(A-LC) = %s\n\n',mat2str(eig(A-L*C).',4));
end

%% 5) Simulacao: estimador em malha aberta (L=0) vs observador (L)
t=(0:1e-3:4)'; u=0.5*sin(2*pi*t);
x0=[0.1;0;-0.1]; xh0=[0;0;0];
Aol=blkdiag(A,A);  Acl=[A zeros(n);L*C A-L*C];
Baug=[B;B]; Caug=eye(2*n); Daug=zeros(2*n,1);
Zol=lsim(ss(Aol,Baug,Caug,Daug),u,t,[x0;xh0],'zoh');
Zcl=lsim(ss(Acl,Baug,Caug,Daug),u,t,[x0;xh0],'zoh');
Eol=Zol(:,1:n)-Zol(:,n+1:end); Ecl=Zcl(:,1:n)-Zcl(:,n+1:end);
% Verificacao independente: e(t)=expm((A-LC)t)e(0), a entrada cancela-se.
exact=zeros(size(Ecl));
for j=1:numel(t), exact(j,:)=(expm((A-L*C)*t(j))*(x0-xh0)).'; end
% Verificar ate 2 s: depois a planta (instavel) tem estados enormes e a
% subtracao x - xhat perde precisao numerica.
k=t<=2;
assert(max(abs(Ecl(k,:)-exact(k,:)),[],'all')<1e-9*max(1,max(abs(Zcl(k,:)),[],'all')));
fprintf('|e(%g s)|: malha aberta = %.4g ; observador = %.4g\n\n',t(end), ...
    norm(Eol(end,:)),norm(Ecl(end,:)));

%% 6) Discreto: o MESMO procedimento em z (como no massa-mola da aula 1)
Ts=0.01;
sysd=c2d(ss(A,B,C,D),Ts,'zoh'); Ad=sysd.A; Bd=sysd.B; Cd=sysd.C;
% TODO 4: mapear os polos continuos para z (ver hands-on da aula 1).
% Porque nao se podem usar -10 e -2+-2j diretamente como polos discretos?
pz=nan(size(polos));
discretoOK=~any(isnan(pz));
if discretoOK
    Ld=ganhoPorCoeficientes(Ad,Cd,pz);  % Processo II, agora com z
    Ldp=place(Ad',Cd',pz).';
    assert(max(abs(Ld-Ldp))<1e-6*max(1,norm(Ldp)));
    assert(all(abs(eig(Ad-Ld*Cd))<1),'Polos discretos fora do circulo unitario.');
else
    Ld=nan(n,1); Ldp=Ld; fprintf('Discreto: completar o TODO 4.\n');
end
if discretoOK && any(strcmp(stage,{'discreto','tudo'}))
    fprintf('Discreto (Ts=%g s): pz = %s\n',Ts,mat2str(pz,6));
    fprintf('Ld (coeficientes) = %s\n',mat2str(Ld.',5));
    fprintf('Ld (place)        = %s\n',mat2str(Ldp.',5));
    fprintf('Ts*L (continuo)   = %s   <- so uma aproximacao\n',mat2str(Ts*L.',5));
    fprintf('|eig(Ad-Ld*Cd)|   = %s\n\n',mat2str(abs(eig(Ad-Ld*Cd)).',5));
end

if showPlots && any(strcmp(stage,{'simular','tudo'}))
    figure('Name','Planta instavel: estados','Color','w');
    lab={'x1','x2','x3'}; tz=t<=1.5;
    for s=1:n
        subplot(n,1,s);
        plot(t(tz),Zcl(tz,s),'b',t(tz),Zcl(tz,n+s),'g--',t(tz),Zol(tz,n+s),'r:','LineWidth',1.3);
        grid on; ylabel(lab{s});
        if s==1, title('Planta instavel, sem controlador'); legend('real','observador','malha aberta','Location','best'); end
    end
    xlabel('tempo [s]');
    figure('Name','Planta instavel: norma do erro','Color','w');
    semilogy(t,max(vecnorm(Eol,2,2),eps),'r',t,max(vecnorm(Ecl,2,2),eps),'g','LineWidth',1.4);
    grid on; xlabel('tempo [s]'); ylabel('||x - x_{hat}||');
    legend('malha aberta (L=0): diverge','observador (polos -10, -2\pm2j)','Location','best');
    title('A malha aberta so convergia na aula 1 porque a planta era estavel');
end
if showPlots && discretoOK && any(strcmp(stage,{'discreto','tudo'}))
    figure('Name','Polos: continuo vs discreto','Color','w');
    subplot(1,2,1); plot(real(polos),imag(polos),'gx','MarkerSize',10,'LineWidth',2); hold on;
    plot(real(eig(A)),imag(eig(A)),'ro','MarkerSize',8); xline(0,'k');
    grid on; axis equal; title('plano s'); legend('A-LC','A','Location','best');
    subplot(1,2,2); th=linspace(0,2*pi,200); plot(cos(th),sin(th),'k'); hold on;
    plot(real(pz),imag(pz),'gx','MarkerSize',10,'LineWidth',2);
    plot(real(eig(Ad)),imag(eig(Ad)),'ro','MarkerSize',8);
    grid on; axis equal; title(sprintf('plano z (Ts=%g s)',Ts));
end
out=struct('A',A,'C',C,'polos',polos,'O',O,'Lo',Lo,'T',T,'M',M,'Ls',Ls,'L',L, ...
    't',t,'u',u,'Zol',Zol,'Zcl',Zcl,'Eol',Eol,'Ecl',Ecl,'Ad',Ad,'Cd',Cd,'pz',pz,'Ld',Ld);
end

function [L,M]=ganhoPorCoeficientes(A,C,polos)
% Igualar os coeficientes de det(xI-A+LC) aos do polinomio desejado.
% Serve para s (continuo) e para z (discreto). Apenas uma saida (C 1xn).
n=size(A,1);
pc=real(poly(A)); pd=real(poly(polos));
M=zeros(n);
for i=1:n
    e=zeros(n,1); e(i)=1;
    pe=real(poly(A-e*C));
    M(:,i)=(pe(2:end)-pc(2:end)).';
end
L=M\(pd(2:end)-pc(2:end)).';
end
