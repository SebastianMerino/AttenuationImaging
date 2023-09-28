
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sintaxis: [out]=lectura
% Salidas: out.RF = Datos RF
%          out.Bmode = Datos de brillo
%          out.fs = Frecuencia de muestreo
%          out.fc = Frecuencia central
%          out.x,out.z = Dimensiones laterales y axiales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out]=lectura_OK(nombre)

%nombre = input('Escriba el nombre de archivo: ','s');

[RF, Settings] = RPread(nombre, '6.0.5');


z=1/Settings.sf*(1:Settings.h)/2*1540;
x=linspace(0,0.038,size(RF,2));

h = waitbar(0, 'Generando las imagenes de brillo');
for index=1:size(RF,3)
US = 20*log10(abs(hilbert(RF(:,:,index)))); 
iUS(:,:,index) = US - max(US(:));
waitbar(index/size(RF,3),h)
end
close(h);

out.RF=RF;
out.Bmode=iUS;
out.fs=Settings.sf;
out.fc=Settings.txf;
out.x=x;
out.z=z;
end