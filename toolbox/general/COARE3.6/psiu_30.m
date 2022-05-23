
function psi=psiuo(zet)

	x=(1-15*zet).^.25;
	psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
	x=(1-10.15*zet).^.3333;
	psic=1.5*log((1+x+x.*x)/3)-sqrt(3)*atan((1+2*x)/sqrt(3))+4*atan(1)/sqrt(3);
	f=zet.*zet./(1+zet.*zet);
	psi=(1-f).*psik+f.*psic;                                               
   ii=find(zet>0);
   if ~isempty(ii);

	%psi(ii)=-4.7*zet(ii);
  	%c(ii)=min(50,.35*zet(ii));
   c=min(50,.35*zet);
	psi(ii)=-((1+1.0*zet(ii)).^1.0+.667*(zet(ii)-14.28)./exp(c(ii))+8.525);
	end;
