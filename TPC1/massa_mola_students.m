function out = massa_mola_students(cfg)
% Observadores CA | massa-mola-amortecedor | MATLAB R2020b+.
% Requer Control System Toolbox. Executar: massa_mola_students
% Automatizar: r=massa_mola_students(struct('animate',false,'plots',false));
if nargin<1, cfg=struct; end
cfg=defaults(cfg,struct('Ts',0.01,'T',10,'poles',[-3 -4], ...
    'sigma',0,'mReal',1,'cReal',1.5,'kReal',20, ...
    'x0',[0.15;0],'xh0',[0;0],'animate',true,'plots',true, ...
    'saveVideo',true,'videoFile','anim_observador.mp4','seed',7));
assert(cfg.Ts>0 && cfg.T>2,'Usar Ts>0 e T>2.');
rng(cfg.seed);

m=1; c=1.5; k=20;
A=[0    1;
   -k/m -c/m]; 
B=[0;
   1/m]; 
C=[1 0]; 
D=0;
n=2;

sysd=c2d(ss(A,B,C,D),cfg.Ts,'zoh');
Ad=sysd.A; Bd=sysd.B; Cd=sysd.C;
realD=c2d(ss([0                    1; ...
              -cfg.kReal/cfg.mReal -cfg.cReal/cfg.mReal], ...
             [0;1/cfg.mReal], ...
             C, ...
             0),cfg.Ts,'zoh');
t=0:cfg.Ts:cfg.T; 
N=numel(t); 
u=2*sin(2*pi*0.3*t);

%% BLOCO 1: observabilidade (modelo continuo).
% TODO 1: construir O e verificar se rank(O)==2.

O = obsv(A, C) 
r = rank(O)
if r == n;
    fprintf('O sistema é observável\n')
else
    fprintf('O sistema não é observável\n');
end

%% BLOCO 2: polos do observador preditor discreto.
% TODO 2: mapear cfg.poles para pd=exp(...), depois usar place.

pd=exp(cfg.poles*cfg.Ts)
L=place(Ad',Cd',pd)'
ei = eig(Ad-L*Cd)
if abs(sort(ei') - sort(pd)) < 1e-10
    fprintf('Os pólos coincidem\n')
else
    fprintf('Os pólos não coincidem\n')
end

X=zeros(2,N); H=X; OL=X; Y=zeros(1,N); innovation=zeros(1,N);
X(:,1)=cfg.x0; H(:,1)=cfg.xh0; OL(:,1)=cfg.xh0;
noise=cfg.sigma*randn(1,N);
for j=1:N
    % Primeiro guardar/medir os estados na MESMA amostra j.
    Y(j)=Cd*X(:,j)+noise(j);
    %% BLOCO 3: unica medida disponivel ao observador: Y(j).
    % TODO 3: substituir 0 pela inovacao y_j - Cd*xhat_j.
    
    innovation(j)=Y(j)-Cd*H(:,j);
    if j<N
        X(:,j+1)=realD.A*X(:,j)+realD.B*u(j);
        OL(:,j+1)=Ad*OL(:,j)+Bd*u(j);
        %% BLOCO 4: observador preditor. Nao usar X na atualizacao.
        % TODO 4: acrescentar o termo de correcao L*innovation(j).
        
        H(:,j+1)=Ad*H(:,j)+Bd*u(j)+L*innovation(j);
    end
end
E=X-H; Eol=X-OL; mask=t>=2;
rmse=sqrt(mean(E(:,mask).^2,2));
fprintf('Posto O: %d/2; RMSE apos 2 s: pos=%.5g m, vel=%.5g m/s\n', ...
    rank(obsv(A,C)),rmse(1),rmse(2));

if cfg.plots
    figure('Name','Estados e erro','Color','k');
    labels={'posicao [m]','velocidade [m/s]'};
    for s=1:2
        subplot(3,1,s); plot(t,X(s,:),'b',t,H(s,:),'g--',t,OL(s,:),'r:','LineWidth',1.2);
        grid on; ylabel(labels{s}); legend('real','observador','malha aberta');
    end
    subplot(3,1,3); plot(t,E,'LineWidth',1.2);grid on;
    xlabel('tempo [s]');ylabel('erro');legend('e posicao [m]','e velocidade [m/s]');
    figure('Name','Medicao e inovacao','Color','k');
    subplot(2,1,1);plot(t,Y,'.','Color',[.6 .6 .6]);hold on;plot(t,H(1,:),'g');
    ylabel('posicao [m]');legend('medicao','estimativa');grid on;
    subplot(2,1,2);plot(t,innovation);ylabel('y - C xhat [m]');xlabel('tempo [s]');grid on;
end
if cfg.animate
    animate(t,X,H,OL,cfg);
end
out=struct('t',t,'u',u,'X',X,'H',H,'OL',OL,'Y',Y,'E',E, ...
    'Eol',Eol,'innovation',innovation,'L',L,'Ad',Ad,'Cd',Cd, ...
    'poles',eig(Ad-L*Cd),'rmse',rmse,'cfg',cfg);
end
function c=defaults(c,d)
f=fieldnames(d);for j=1:numel(f),if ~isfield(c,f{j}),c.(f{j})=d.(f{j});end,end
end
function animate(t,X,H,OL,cfg)
f=figure('Name','Massa-mola: real e estimativas','Color','k','Position',[80 80 1000 420]);
ax=axes(f);hold(ax,'on');grid(ax,'on');
base=0.65; bw=.10; bh=.12;
allp=base+[X(1,:) H(1,:) OL(1,:)];
xlim([min(-.08,min(allp)-.18),max(1.1,max(allp)+.18)]);ylim([-.28 .28]);
plot([0 0],[-.20 .20],'k','LineWidth',3);
block=rectangle('Position',[base+X(1,1)-bw/2,-bh/2,bw,bh], ...
    'FaceColor',[.5 .7 1],'EdgeColor','b');
spring=plot(nan,nan,'b','LineWidth',1.5);
damper=plot(nan,nan,'k','LineWidth',2);
h=plot([base base],[-.15 .15],'g--','LineWidth',2);
o=plot([base base],[-.20 .20],'r:','LineWidth',2);
legend([h o],{'estimativa','malha aberta'},'Location','southoutside');
xlabel('posicao visual = 0.65 + deslocamento [m]');set(ax,'YTick',[]);
writer=[];
if cfg.saveVideo
    writer=VideoWriter(cfg.videoFile,'MPEG-4');writer.FrameRate=30;open(writer);
    cleanup=onCleanup(@()close(writer)); 
end
% 30 frames/s; cada frame representa o tempo indicado no titulo.
frames=unique([round(linspace(1,numel(t),round(t(end)*30)+1)) numel(t)]);
clock0=tic;
for j=frames
    if ~isgraphics(f),break,end
    pos=base+X(1,j);block.Position=[pos-bw/2,-bh/2,bw,bh];
    sx=linspace(.03,pos-bw/2,25);sy=.025*(-1).^(0:24);sy([1 end])=0;
    set(spring,'XData',sx,'YData',sy);
    set(damper,'XData',[0 .12 .12 pos-bw/2],'YData',[-.13 -.13 -.10 -.10]);
    set(h,'XData',base+H(1,j)*[1 1]);set(o,'XData',base+OL(1,j)*[1 1]);
    title(sprintf('t = %.2f s | azul: massa real | verde: observador | vermelho: malha aberta',t(j)));
    drawnow;
    if cfg.saveVideo,writeVideo(writer,getframe(f));else,pause(max(0,t(j)-toc(clock0)));end
end
end
