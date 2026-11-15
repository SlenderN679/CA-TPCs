function Animate_BaB(X,t,Xhat)
% Animacao da bola-na-viga.
%   Animate_BaB(X,t)        -> so o sistema real (como antes)
%   Animate_BaB(X,t,Xhat)   -> real (azul/verde) + estimativa do observador (vermelho)
% X e Xhat: N x 4, colunas [x x' theta theta'] (posicao em m, angulo em rad).
% Desenho a escala: viga de 1.2 m; a posicao da bola e visivel.
if nargin<3, Xhat=[]; end
Lb=1.2;  R=0.04;                 % comprimento da viga e raio da bola [m]
f=gcf; ax=gca; cla(ax); hold(ax,'on'); grid(ax,'on');
axis(ax,'equal'); axis(ax,[-0.8 0.8 -0.55 0.55]);
plot(ax,[-Lb/2 Lb/2],[0 0],'k-.');
plot(ax,0,0,'ko','MarkerFaceColor','k');
th=linspace(0,2*pi,60);
hViga=plot(ax,nan,nan,'b-','LineWidth',5);
hBola=plot(ax,nan,nan,'g-','LineWidth',3);
if ~isempty(Xhat)
    hVigaH=plot(ax,nan,nan,'r--','LineWidth',1.5);
    hBolaH=plot(ax,nan,nan,'r--','LineWidth',2);
    legend(ax,[hBola hBolaH],{'real','estimativa do observador'},'Location','southoutside');
end
xlabel(ax,'[m]');
% ~30 imagens por segundo de simulacao, ao ritmo do tempo real
idx=unique([round(linspace(1,numel(t),max(2,round((t(end)-t(1))*30)+1))) numel(t)]);
clock0=tic;
for i=idx
    if ~isgraphics(f), break, end
    desenha(hViga,hBola,X(i,1),X(i,3),Lb,R,th);
    if ~isempty(Xhat), desenha(hVigaH,hBolaH,Xhat(i,1),Xhat(i,3),Lb,R,th); end
    txt=sprintf('t = %.2f s | x = %.3f m | \\theta = %.1f graus',t(i),X(i,1),rad2deg(X(i,3)));
    if ~isempty(Xhat)
        txt=[txt sprintf(' | erro x = %.3f m',X(i,1)-Xhat(i,1))]; %#ok<AGROW>
    end
    title(ax,txt);
    drawnow;
    pause(max(0,(t(i)-t(1))-toc(clock0)));
end
end

function desenha(hViga,hBola,x,teta,Lb,R,th)
set(hViga,'XData',Lb/2*cos(teta)*[-1 1],'YData',Lb/2*sin(teta)*[-1 1]);
cx=x*cos(teta)-R*sin(teta); cy=x*sin(teta)+R*cos(teta);   % bola em cima da viga
set(hBola,'XData',cx+R*cos(th),'YData',cy+R*sin(th));
end
