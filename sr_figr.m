
function h = sr_figr(W,H,xs,ys)
% W : width (inches)
% H : hight (inches)
% xs : x-loc
% ys : y-loc

scCorr_stz = get(0,'ScreenSize');
defaultSCR=[scCorr_stz(3) scCorr_stz(4)];
SCR = defaultSCR;
nxs = xs;
nys = ys;
WH = [W H];
asprat = WH(1)/WH(2);
h=figure;
US = get(h,'papersize');
xp = (US(1)-WH(1))/2;
yp = (US(2)-WH(2))/2;
xs = SCR(2)*xp/US(2);
ys = SCR(2)*yp/US(2);
Hs = SCR(2)*WH(2)/US(2);
Ws = Hs*asprat;

set(h,'position',[nxs nys Ws Hs],'paperposition',[xp yp WH]);

end
