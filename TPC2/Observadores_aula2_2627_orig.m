function out = Observadores_aula2_2627(cfg)
% Bola-na-viga: realimentacao do estado estimado, u = -K*xhat.
% Exercicio 1 da aula 2. Completar apenas os TODOs entre as linhas EDITAR AQUI.
% Requer Control System Toolbox e Symbolic Math Toolbox. MATLAB R2020b+.
%   Executar:      Observadores_aula2_2627
%   Sem janelas:   r=Observadores_aula2_2627(struct('plots',false,'animate',false));
%   Outros polos:  Observadores_aula2_2627(struct('polosL',[-1 -1 -1 -1]));
% Estados: x1 = posicao da bola [m], x2 = velocidade [m/s],
%          x3 = angulo da viga [rad], x4 = velocidade angular [rad/s].
% Mede-se apenas a posicao da bola: y = x1.
if nargin<1, cfg=struct; end
cfg=defaults(cfg,struct('polosK',[-2+1j -2-1j -2 -2],'polosL',[-5 -5 -5 -5], ...
    'x0',[0.1;0.5;deg2rad(-5);deg2rad(2)],'xh0',zeros(4,1), ...
    'T',8,'dt',0.01,'plots',true,'animate',true));
p=struct('m',0.2,'r',0.05,'g',9.81,'I',2e-4,'J',2,'mu',0);  % m, r: massa e raio da bola; J: inercia da viga; mu: atrito
[A,B,C,D]=findStateSpace(p); %#ok<ASGLU>
n=size(A,1);
fprintf('Modelo linearizado em x=0, u=0:\nA =\n');disp(A);fprintf('B =\n');disp(B);
fprintf('eig(A) = %s\n',mat2str(eig(A).',4));

%% EDITAR AQUI ---------------------------------------------------------
% CONTROLADOR (revisao da parte anterior da UC): u = -K*x
% TODO 1: verificar a controlabilidade de (A,B) e calcular K com os polos
%         cfg.polosK. Atencao: ha polos repetidos -> usar acker, nao place.
K=zeros(1,n);

% OBSERVADOR
% TODO 2: verificar a observabilidade de (A,C) e calcular L com os polos
%         cfg.polosL, usando a dualidade (A',C'). Quatro polos iguais -> acker.
L=zeros(n,1);
%% ----------------------------------------------------------------------

if ~any(K) || ~any(L)
    warning('K e/ou L ainda a zero: completar os TODOs 1 e 2. Sem K a bola cai (planta instavel).');
end
% Polos repetidos: comparar polinomios (os valores proprios de uma raiz
% multipla sao numericamente sensiveis).
if any(K), assert(norm(poly(A-B*K)-poly(cfg.polosK))<1e-6*norm(poly(cfg.polosK)),'Rever K.'); end
if any(L), assert(norm(poly(A-L*C)-poly(cfg.polosL))<1e-6*norm(poly(cfg.polosL)),'Rever L.'); end
fprintf('K = %s\nL = %s\n',mat2str(K,5),mat2str(L.',5));

%% Simulacao (modelo linear): planta + observador, u = -K*xhat
t=(0:cfg.dt:cfg.T)';
opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,XX]=ode45(@(t,z) sistema(z,A,B,C,K,L),t,[cfg.x0;cfg.xh0],opts);
X=XX(:,1:n); Xhat=XX(:,n+1:end); E=X-Xhat; u=-(K*Xhat.').';
fprintf('|e(%g s)| = %.3g ; |x(%g s)| = %.3g ; max|u| = %.3g N.m\n', ...
    t(end),norm(E(end,:)),t(end),norm(X(end,:)),max(abs(u)));

if cfg.plots
    nomes={'x  [m]','x'' [m/s]','\theta [graus]','\theta'' [graus/s]'};
    esc=[1 1 180/pi 180/pi];          % angulos em graus apenas nos graficos
    figure('Name','Estados reais e estimados','Color','w');
    for s=1:n
        subplot(2,2,s);
        plot(t,esc(s)*X(:,s),'b',t,esc(s)*Xhat(:,s),'r--','LineWidth',1.3);
        grid on; ylabel(nomes{s}); xlabel('tempo [s]');
        if s==1, legend('real','observador','Location','best'); end
    end
    figure('Name','Erro de estimacao e entrada','Color','w');
    subplot(2,1,1); semilogy(t,max(vecnorm(E,2,2),eps),'LineWidth',1.3); grid on;
    ylabel('||x - x_{hat}||'); title(sprintf('polos do observador: %s',mat2str(cfg.polosL)));
    subplot(2,1,2); plot(t,u,'LineWidth',1.3); grid on;
    ylabel('u = -K x_{hat} [N.m]'); xlabel('tempo [s]');
end
if cfg.animate
    figure('Name','Bola-na-viga: real (verde) e estimativa (vermelho)','Color','w');
    Animate_BaB(X,t,Xhat);
end
out=struct('A',A,'B',B,'C',C,'K',K,'L',L,'t',t,'X',X,'Xhat',Xhat,'E',E,'u',u,'cfg',cfg);
end

function dz = sistema(z,A,B,C,K,L)
n=size(A,1);
x=z(1:n); xhat=z(n+1:end);
u=-K*xhat;                       % o controlador so conhece xhat
y=C*x;                           % unica medida: posicao da bola
dx=A*x+B*u;
dxhat=A*xhat+B*u+L*(y-C*xhat);   % observador de Luenberger
dz=[dx;dxhat];
end

function [A,B,C,D] = findStateSpace(p)
% Linearizacao do modelo nao linear em x=0, u=0.
syms x1 x2 x3 x4 u real
m=p.m; g=p.g; r=p.r; I=p.I; J=p.J; mu=p.mu;
% Viga: binario do peso da bola. O centro da bola esta em (x1, r) no
% referencial da viga (r = raio da bola), logo o binario e
% m*g*(x1*cos(x3) + r*sin(x3)). O termo em x1 e o dominante.
f1=(m*g*(x1*cos(x3)+r*sin(x3))-u)/(J+I+m*(x1^2+r^2));   % aceleracao angular da viga
f2=m*g*cos(x3)-m*(x1*f1+2*x2*x4);           % forca normal
dx=[x2;
    x1*x4^2+g*sin(x3)-mu/m*f2;              % mu: coeficiente de atrito (0 aqui)
    x4;
    f1];
x=[x1;x2;x3;x4];
A=double(subs(jacobian(dx,x),[x1 x2 x3 x4 u],[0 0 0 0 0]));
B=double(subs(jacobian(dx,u),[x1 x2 x3 x4 u],[0 0 0 0 0]));
C=[1 0 0 0];
D=0;
end

function c=defaults(c,d)
f=fieldnames(d);for j=1:numel(f),if ~isfield(c,f{j}),c.(f{j})=d.(f{j});end,end
end
