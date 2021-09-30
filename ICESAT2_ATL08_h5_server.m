% Copyright (C) 2019 by The HDF Group.
% All rights reserved.
%
%  This example code illustrates how to access and visualize
%  ICESat-2 ATL06 HDF5 file in MATLAB. 
%
%  If you have any questions, suggestions, comments on this example, please 
% use the HDF-EOS Forum (http://hdfeos.org/forums). 
%
%  If you would like to see an  example of any other NASA HDF/HDF-EOS data 
% product that is not listed in the HDF-EOS Comprehensive Examples page 
% (http://hdfeos.org/zoo), feel free to contact us at eoshelp@hdfgroup.org or 
% post it at the HDF-EOS Forum (http://hdfeos.org/forums).
% 
% Usage:save this script and run (without .m at the end)
%
%  $matlab -nosplash -nodesktop -r ATL06_20190223232535_08780212_001_01_h5
% 
%

%addpath('/home/nico/Documents/matlab_script')

% Tested under: MATLAB R2018a
% Last updated: 2020-09-10
%cd /home/nico/Documents/cecs/ICEsat_2/dat_cecs_lake
%cd /home/nico/Documents/dat_cecs_lake/datos_icesat2_julio_2019/
%cd /media/nico/data_nicod_/icesat2

close all; clear all
format long;

path_icesat2_files = '/home/ndonso/ice_sat2/python_icesat2/olivares/data';   % where are icesat2 .h5 files
math_file          = 'Olivares_septiembre2021_icesat2_ATL08.mat';    % name of matlab file that contain all data 
math_filter        = 'Olivares_septiembre2021_icesat2_ATL08_filtrado.mat'; % name of mathlab file that contain filtered data

cd(path_icesat2_files);

nom=dir('*.h5');
nombre_proyecto = 'Olivares_ATL08_sept2021'; %give a name for your project 
print_figures = 0;% 1 (yes) or 0 (not)
%############################################ IMPORTANTE
%pasos: opciones
%1 Read and save data from h5 files
%2 Graficar trazas del satélite sin filtrar. Se realiza por partes para no saturar la memoria
%3 Filtrado de datos para hacerlos mas manejables
%4 Gráficos
%############################################
%##### "paso" DEFINE EL SWITCH CASE
paso = 3;
disp(['Se ejecutará el paso Nº',num2str(paso)])
%############################################
%Vetices del entorno general del área de estudio, se tomaron desde google earth. Se utilizan para recortar los datos de ICESAT-2 y hacerlos más manejables.

v_lat = [-33.174692,-33.289772,-33.237035,-33.045642];       %vertices lat
v_lon = [-70.298994,-70.266591,-70.040984,-70.099149];   %vertices lon

%poligono del área de estudio definido, obtenido de un GIS

area_lat = [-33.174692,-33.289772,-33.237035,-33.045642];
area_lon = [-70.298994,-70.266591,-70.040984,-70.099149];

%% Mensajes a imprimir
if print_figures == 1
    disp('Se guaradarán las figuras')
elseif print_figures == 0
    disp('No se guardarán las figuras')
end
%%
switch paso 
    case 1 %1 Read and save data from h5 files
	for ii = 1:length(nom)
		disp(['File number: ',num2str(ii)]); 
    
		% Open the HDF5 File.
		FILE_NAME = nom(ii).name;
		file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
		% Open the datasets.
		%LATFIELD_NAME='gt1l/land_ice_segments/latitude';
		LATFIELD_NAME='gt1l/land_segments/latitude';
		lat_id=H5D.open(file_id, LATFIELD_NAME);

		%LONFIELD_NAME='gt1l/land_ice_segments/longitude';
		LONFIELD_NAME='gt1l/land_segments/longitude';
		lon_id=H5D.open(file_id, LONFIELD_NAME);

		%DATAFIELD_NAME='gt1r/land_ice_segments/h_li';
		DATAFIELD_NAME='gt1l/land_segments/dem_h';
		data_id=H5D.open(file_id, DATAFIELD_NAME);

		% Read the datasets.
		latitude=H5D.read(lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL',...
                'H5P_DEFAULT');
		lon=H5D.read(lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
             	'H5P_DEFAULT');
		data=H5D.read(data_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
             	'H5P_DEFAULT');

		% Read the attributes.
		ATTRIBUTE = 'units';
		attr_id = H5A.open_name (data_id, ATTRIBUTE);
		units_data = H5A.read(attr_id, 'H5ML_DEFAULT');

		ATTRIBUTE = 'long_name';
		attr_id = H5A.open_name (data_id, ATTRIBUTE);
		long_name_data = H5A.read(attr_id, 'H5ML_DEFAULT');

		% Read the fill value attribute.
		ATTRIBUTE = '_FillValue';
		attr_id = H5A.open_name(data_id, ATTRIBUTE);
		fillvalue=H5A.read(attr_id, 'H5T_NATIVE_DOUBLE');

		% Close and release resources.
		H5A.close (attr_id)
		H5D.close (data_id);
		H5D.close (lon_id);
		H5D.close (lat_id);
		H5F.close (file_id);

		% Replace the filled value with NaN.
		data(data==fillvalue) = NaN;

		%Store data in A matrix

		A(1:length(data),ii) = data;
		LAT(1:length(data),ii) = latitude;
		LON(1:length(data),ii) = lon;
		names_file(ii)={nom(ii).name};
		end

	save(math_file,'A','LAT','LON','names_file')
	disp('Done!, go to step 2 or 3')

	case 2 %2 Graficar trazas del satélite sin filtrar. Se realiza por partes para no saturar la memoria
	load('icesat2_lago_data_V2_14_02_2020.mat')
	%% Figura de la serie completa por partes, debido a que se llena la memoria del computador si se hace completa
	% Figura creadas y guardadas

	f = figure(1)
	axesm('MapProjection','eqdcylin','Frame','on','Grid','on', ...
     	'MeridianLabel','on','ParallelLabel','on','MLabelParallel','south')
	for ss = 100:150;%length(nom) %No puedo correr todas las trazas, se pega
		disp(['S: ',num2str(ss)])
 		hold on
 		scatterm(LAT(:,ss),LON(:,ss),1,A(:,ss))
	end 
	coast = load('coast.mat');
	plotm(coast.lat, coast.long, 'k');
	hold on
	plotm(-79.25,-87.56,'or')
	tightmap;
	toc
	saveas(f, 'Traza_100_150.png');

	% Create the graphics figure.
	f = figure('Name', FILE_NAME, ...
           'Renderer', 'zbuffer', ...
           'Position', [0,0,800,600], ...
           'visible','off');

	% Put title.
	var_name = sprintf('%s', long_name_data);
	tstring = {FILE_NAME;var_name; ii};
	title(tstring,...
      	'Interpreter', 'none', 'FontSize', 16, ...
      	'FontWeight','bold');
	axesm('MapProjection','eqdcylin','Frame','on','Grid','on', ...
      	'MeridianLabel','on','ParallelLabel','on','MLabelParallel','south')

	% Plot data.
	%scatterm(latitude, lon, 1, data);
	plotm(latitude, lon, data);
	hold all
	h = colorbar();
	units_str = sprintf('%s', char(units_data));
	set (get(h, 'title'), 'string', units_str, 'FontSize', 8, ...
                   'Interpreter', 'None', ...
                   'FontWeight','bold');

	% Plot world map coast line.
	coast = load('coast.mat');
	plotm(coast.lat, coast.long, 'k');
	hold on
	plotm(-79.25,-87.56,'or')
	tightmap;

	%saveas(f, 'FINAL.m.png');
 
%% Filtrado e datos
    case 3 %3 Filtrado de datos para hacerlos mas manejables
	disp('############################################################')
	disp('############################################################')
	disp('######################### PASO 3 ###########################')
	cd(path_icesat2_files)
	load(math_file);
        [a,b]=size(A);% obtengo tamaño de matriz A
	
	for ii=1:b
    		for jj=1:a
        		disp(['ii,jj = ',num2str(ii),'/',num2str(b),' , ',num2str(jj),'/',num2str(a)]);
        		in(jj,ii) = inpolygon(LAT(jj,ii),LON(jj,ii),v_lat,v_lon);%es la matriz similar a A que tiene 1 y 0 logicos indicando si están dentro o fuera del poligono
    		end
	end


	%% Filtramos los datos que están dentro del lago datos por los puntos lake_lon y lake_lan
	[a,b]=size(A);% obtengo tamaño de matriz A
	for ii=1:b
    		for jj=1:a
        		%disp(['ii,jj dentro del lago = ',num2str(ii),',',num2str(jj)]);
			disp(['ii,jj = ',num2str(ii),'/',num2str(b),' , ',num2str(jj),'/',num2str(a)]);
        		in_lake(jj,ii) = inpolygon(LAT(jj,ii),LON(jj,ii),area_lat,area_lon);%es la matriz similar a A que tiene 1 y 0 logicos indicando si están dentro o fuera del poligono
    		end
	end


	disp(['Se encontraron ',num2str(sum(sum(in))),' puntos dentro del poligono.']);
	disp(['Se encontraron ',num2str(sum(sum(in_lake))),' dentro de ',nombre_proyecto]);

	save(math_filter,'in','A','LAT','LON','in_lake','names_file');

%#####################################################################################################
%% Gráficos




    case 4
	disp('Paso 4')
	%load(math_file)
	load('icesat2_lago_data_fitrado_V2.mat')
	% %graficamos todas las tramas en el mapa mundial. NO SE VE NADA BIEN
	% ff1 = figure(1)
	% axesm('MapProjection','eqdcylin','Frame','on','Grid','on', ...
	%      'MeridianLabel','on','ParallelLabel','on','MLabelParallel','south')
	% for ss = 1:length(nom) 
	% disp(['S: ',num2str(ss)])
	% hold on
	% scatterm(LAT(in(:,ss),ss),LON(in(:,ss),ss),1,A(in(:,ss),ss))
	% end 
	% coast = load('coast.mat');
	% plotm(coast.lat, coast.long, 'k');
	% hold on
	% plotm(-79.25,-87.56,'or');% centro del lago
	% plotm(lake_lat,lake_lon,'-b');% perimetro
	% h = colorbar;

%#############################################################################################
	%% Grafico con datos solo dentro del lago
	
	%Buscamos los nombres de los archivos que tienen datos dentro del lago
	nombres = names_file(any(in_lake));%guardo los nombres en "nombres"
	disp(['Los archivos con datos dentro del LAGO son: ',nombres])

	disp(['Los indices que tienen datos dentro de la zona de estudio son: ',num2str(find(any(in_lake)))])

	%% Eliminacion de datos
	% Si se detectó una medicion dentro del área de estudio se puede eliminar
	% indicando el indice 
	
	%indx_erase = [126];
	%in_lake(:,indx_erase) = false;
	%indx_erase = [126];
	%in_lake(:,indx_erase) = false;
	
	%disp(['Se eliminó las mediciones del día: ',names_file(indx_erase)])



	%%

	for i = 1:length(nombres)
		name(i).nombre = nombres{i};
	end
%name(1).nombre(7:14)% me da las fechas del archivo
%disp([name(1).nombre(7:10),'/',name(1).nombre(11:12),'/',name(1).nombre(13:14)]);% me da: 2018/11/17

%cd('/home/nico/Documents/dat_cecs_lake/Datos_gps')
%M = csvread('datos_gps.csv');
%cd('/media/nico/data_nicod_/icesat2')% volvemos a donde estabamos

	ff2 = figure(2);
	clf
	worldmap([-79.3 -79.19],[-87.8 -87.3])
	%scatterm(M(:,1),M(:,2),15,M(:,3),'fill')%Datos GPS de balizas puntos
	hold on

	for ss = 1:length(nom) 
		%disp(['S: ',num2str(ss)])
		hold on
		scatterm(LAT(in_lake(:,ss),ss),LON(in_lake(:,ss),ss),5,A(in_lake(:,ss),ss))
	end 

	plotm(lake_lat,lake_lon);

	%plotm(-79.25,-87.56,'or');% centro del lago
	plotm(lake_lat,lake_lon,'-b');% perimetro
	h = colorbar;
	caxis([2020 2030])
	scaleruler
	ylabel(h,'m.a.s.l.');


if print_figures == 1
print(ff2,'-djpeg','ICESAT_2_DENTRO_LAGO_CECs_2','-r1000')
end
%###############################################################################################



%Gráfico de los valosres de altura
	aa1 = 1;
	ff3=figure(3);
	clf
	parfor ss = 1:length(nom) 
		disp(['plot in: ',num2str(ss)])
		hold on
		plot(A(in_lake(:,ss),ss),'.','LineWidth',aa1)%,'Color',[0.4 0.4 0.4])
		ylim([2018 2032])% elimina del gráfico valores que están fuera de rango
		xlabel('Number of measurements by satellite pass','FontSize',15)
		ylabel('m.a.s.l.','FontSize',15)
	end 
	grid on
	legend([name(1).nombre(7:10),'/',name(1).nombre(11:12),'/',name(1).nombre(13:14)],[name(2).nombre(7:10),'/',name(2).nombre(11:12),'/',name(2).nombre(13:14)],[name(3).nombre(7:10),'/',name(3).nombre(11:12),'/',name(3).nombre(13:14)],[name(4).nombre(7:10),'/',name(4).nombre(11:12),'/',name(4).nombre(13:14)],'Location','southeast')

if print_figures == 1
print(ff3,'-djpeg','Mediciones_ICESAT_2_DENTRO_LAGO_CECs_2','-r1000')
end

%###############################################################################################



%Gráfico de los valosres de altura, contra latitud
aa1 = 1;
ff4=figure(4);
clf
for ss = 1:length(nom) 
disp(['plot in: ',num2str(ss)])
hold on
plot(LAT(in_lake(:,ss),ss),A(in_lake(:,ss),ss),'.','LineWidth',aa1)%
ylim([2018 2032])
xlabel('Latitude','FontSize',15)
ylabel('m.a.s.l.','FontSize',15)
end 
grid on
legend([name(1).nombre(7:10),'/',name(1).nombre(11:12),'/',name(1).nombre(13:14)],[name(2).nombre(7:10),'/',name(2).nombre(11:12),'/',name(2).nombre(13:14)],[name(3).nombre(7:10),'/',name(3).nombre(11:12),'/',name(3).nombre(13:14)],[name(4).nombre(7:10),'/',name(4).nombre(11:12),'/',name(4).nombre(13:14)],'Location','southeast')

if print_figures == 1
print(ff4,'-djpeg','Mediciones_ICESAT_2_DENTRO_LAGO_CECs_latitude_2','-r1000')
end


% %% Trabajando con los nombre de los archivos
% 
% nombres_archivos_usados =  {names_file(any(in_lake))};
% 
% nombres_archivos_usados =  names_file(any(in_lake)); 
% 
% for xx = 1:length(names_file(any(in_lake)))
% end


%% Realizo un análisis de las balizas
%convierto a UTM
%[X_sat,Y_sat]=ll2utm(LAT,LON);%ICE_SAT2
%[X_gps,Y_gps]=ll2utm(M(:,1),M(:,2));%

%M(:,4) = ones(length(M),1);
% for jj = 1:length(M)
%     M(jj,4) = sqrt((X)^2+()^2); 
% end



end %switch
