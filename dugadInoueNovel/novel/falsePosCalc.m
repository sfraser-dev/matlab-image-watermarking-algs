function total=falsePosCalc(T,Nw);
%T=0.35;
%Nw=140;

format long
m=floor (Nw*(T+1)/2);
total=0.0;

for i=m:Nw
	% The maximum factorial value in matlab is: factorial(170)
	FACTVAL = factorial(Nw) / (factorial(i) * factorial(Nw-i)) ;
	total = total + (FACTVAL * (0.5 ^ Nw));
	%fprintf('%.10f\n',total);
end; clear i;

while 0
	fprintf('T=%.2f, Nw=%d\n',T,Nw);
	fprintf('Pfp = %.30f\n',total);
	total
end
