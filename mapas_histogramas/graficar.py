'''
author: mlizarazo
La rutina graficar llama al modulo selec_SCPySIS que construye un reporte de sismicidad bien sea de la base de datos bd (sesisan y seicomp) o ap (archivo propio)
luego al modulo funcion_mapa que elabora el mapa a partir de basemap, arcgis, el wms del SGC, y un shape de fallas relevantes, asi mismo el modulo histogramas 
que construye graficos estadisticos. Esta rutina es un archivo de control de usuario en el que se editan las condiciones de busqueda y parametros a tener en cuenta
en la busqueda y arquitectura del mapa. Como resultado genera un reporte de sismicidad y uno alternativo en formato simple de ZMAP asi como los graficos basicos de la 
busqueda

'''

import os, glob
import matplotlib.pyplot as plt
import numpy as np
from funcion_mapa import *
from selec_SCPySIS import *
from histogramas import *
from mapa_kml import r2kml
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import csv

import pickle

#------borra imagnes y csv's de anteriorres busquedas para no confundir el resultado de la actual consulta
for f in glob.glob("*.png"):
	os.remove(f)

for f in glob.glob("*.csv"):
	os.remove(f)

for f in glob.glob("*.kml"):
	os.remove(f)
#----------------------

#-----------CONDICIONES DE BUSQUEDA (editar coordenadas a las necesarias y elegir bd (base de datos) o ap (archivo propio)
#---Eleccion de base de datos, si el archivo de sismos es propio op="ap", si requiere consultar bd de seisan y seiscomp  op="bd". Si op="ap" ir a linea ~125 para dar la ruta del archivo propio y leer las instrucciones 
op ="bd"#	
#op ="ap"	

incluir_No_loc=False #Incluye los no localizables de Seiscmop. Sirve solamente para consultas de 2018-03-01 en adelante y en opcion base de datos op="bd"

#-----------------Elegir tipo de busqueda
tipo_busqueda="c" #Elegir c: circular, r: rectangular, p: poligono o "sin_filtro_loc"
	#si se elige opcion sin_filtro_loc, esta solo servira cuando la base de datos es propia op = "ap". En este caso, obtiene lat, lon y prof minimas y maximas del propio catalogo, con el fin de que se analice todo el contenido 
	#del catalogo sin imponer restricciones espaciales, sin embargo, si tiene habilitado el filtro por minimos y maximos de tiempo, mag y errores (funciona unicamente para base de datos propia)
	#si se elige la opcion "p" poligono, se busca la sismicidad restringida a dicho poligono. Dar path del poligno (archivo_poligono e linea 64 aprox), ver campos del archivo en la linea 180 aprox

hora_local=False #True o False. Sirve para sacar el report y los histogramas en hora local. Funciona solo para op=bd

#-------------------------------CONDICIONES GEOMETRICAS DE BUSQUEDA		
#------Editar para busquedas rectangulares
lat_min=9.660
lat_max=9.896
lon_min=-73.809
lon_max=-73.575
#------Editar para busquedas circulares
lat_centro=3.85 
lon_centro=-73.68
radio=20 #en kilometros 
#------Editar para busquedas poligono (archivo sin header cada linea es una cordenada (lon lat) separado por espacio. ver campos del archivo en la linea 180 aprox)
#archivo_poligono="/home/mlizarazo/Documentos/codigos_python/consulta_modelo_y_tipo_mag/poligonos_zonas/el_paso.txt"
archivo_poligono="/opt/rutinas/mapas_histogramas/poligonos/Poligono_tecpetro_V3.txt"
#------------------------------------------------------

#------------------------------PARAMETROS de DE BUSQUEDA 

time_i="2011-06-01 00:00:00" #fecha UT
time_f="2023-01-12 23:59:59" #fecha UT

prof_min=-50 #para que tome sismos aereos
prof_max=1200
mag_min=-10
mag_max=10
elat_max=100
elon_max=100
eprof_max=100
#----------------------------

#--------------------------OPCIONES DE RESULTADOS: True o False para graficar fallas y/o estaciones en el mapa generar kml e histogramas, indicar centro de busqueda...
generar_mapa=True #mapa de sismicidad
generar_map_gap=False
generar_map_rms=False
generar_map_time=False
generar_map_calor=False
generar_kml=False

incluir_fallas_sgc=True
incluir_fallas_movmasa=False
graficar_estaciones=True
incluir_rectangulo_leyenda=True #es un rectangulo blanco que se pone abajo del mapa, sirve para resaltar la leyenda
incluir_mapa_ref=True #incluye mapa de referncia en una de las esquinas superiores del mapa, este mapa tiene una escala mayor al mapa original (configurarlo abajo)
indica_centro_busqueda=False #valido solo para busquedas circulares

generar_histrogramas=True

prof_cinco_superficial=False #True para escribir en el reporte de salida en la columna de profundid "superficilal" para los eventos con profundidad menor a 5 km (sirve para las solicitudes), habilitado solo para busquedas en bd
#-------------------------

#--------------------------PARAMETROS DE MAPA

alpha=0.7 #transparencia sismos 0 transparentes 1 solidos

#-----ELECCION DE MAPA BASE----elegir alguna de las opciones de abajo, debe ir entre comillas (ir a http://server.arcgisonline.com/ArcGIS/rest/services a ver actualizaciones y mapas bases )
tipo_mapa="Elevation/World_Hillshade" 
	#favoritos "World_Shaded_Relief" , "NatGeo_World_Map"   "Elevation/World_Hillshade"   "ESRI_Imagery_World_2D"    "ESRI_StreetMap_World_2D"    "NatGeo_World_Map"    "Ocean/World_Ocean_Base"   "World_Imagery"  
	#"World_Physical_Map"   "World_Shaded_Relief" # "World_Street_Map"   "World_Terrain_Base"   "World_Topo_Map"  "Canvas/World_Light_Gray_Base", "Elevation/World_Hillshade"  "Canvas/World_Dark_Gray_Base"

#----CONDICIONES DE LEYENDAS
#---------LOGO SGC
pos_ley_x=-0.02 #se usa para localizar EN X el cuadro del logo del SGC, comienzan en el vertice inferior izquierdo
pos_ley_y=-0.03 #se usa para localizar EN Y el cuadro del logo del SGC, comienzan en el vertice inferior izquierdo
zoom_leyenda=0.71#normalmente tiene valores entre 0.5 y 0.8 (0.72)

#-----mapa de referencia (si esta activo con true)
zoom_map_ref=0.45 #el zoom con el que se pega el mapa de referencia sobre el mapa deseado (buenos valores 0.4 a 0.7), que tan grande es este mapa en la imagen de salida
tam_map_ref=14 #n veces la escala del mapa original (buenos valores 6 a 12)

#----------LOCALIZACION barra de ESCALA 
loc_escala_x=0.78  #localizacion en x, de 0 a 1, 1 es el fin del mapa en la horizontal
loc_escala_y=0.06  #para subir o bajar la barra de escala

#--------- PARAMETROS MAGNITUD, elejir h o v para horizontal o vertical 
FAM=1.6 #factor que amplifica el tamaño de los circulos que representan la magnitud (tipicamente 1.4)
exp_mag=5.2 #exponente que aumenta el tamaño de las pepas (mags_map.append((((mag/8)*(8-FAM))+FAM)**exp_mag --- tipicamente 4.8) se puede disminuir y subir un poco FAM para mapas mas pequeños que logran mejor visibilidad

d_mag=1.05 #acondicionarlo para separar o unir la pepas de la leyenda de magnitud (bajarlo para mas unidadas, subirlo para mas separadas)
y_mag=-0.001 #jugar con el para subir o bajar las pepas de magnitud en la leyenda (pos hacia arriba, neg hacia abajo en grados)
pos_mag_x=0.34 # 0.38 en discreto. localizacion en x, de 0 a 1, 1 es el fin del mapa en la horizontal. Afecta tanto la palabra Magnitud como la cordenda x de incio de las pepas, ajustar la palabra adicionalmente con loc_x_txt)
loc_x_txt=-0.02 #mover la palabra MAGNITD en x respecto a pos_mag_x (en grados), para ajustar
loc_y_txt=0.07 #localizacion en y, de 0 a 0.1 (e grados), de la palabra MAGNITUD
tam_text=30 #palabra MAGNITUD
dir_ley_mag="h"


#---------PARAMETROS PROFUNDIDAD
color_prof="discreto" # "continuo" o "discreto". Grafica los sismos por profundidad asignando un color discreto: rojo, amarillo, azull... de acuerdo a rangos de profundidad, o continuo dando una degrade en cmap a las profundidades
lat_pos_prof=-0.005 #sube o baja las pepas de profundidad (en la leyenda) respecto a la palabra PROFUNDIDAD (aplica en modo discreto)
cmap="cool" #activo para cuando se elige color continuo -- opcines: "cool", "spring", "hot", "summer", "seismic", "coolwarm, "jet" , "turbo", "plasma", "viridis", "cividis"

#---------PARAMETROS DE ESTACIONES
TX=0 #mover en km en X y Y los nombres de las estaciones para su visualizacion (en km) en el MAPA
TY=1
ver_temporales=True
ver_permanentes=True
imprimir_cod_estacion=True #escribe en el mapa el nombre de la estacion
tam_text_est=25 #tamaño del txt en el MAPA del codigo de las estaciones (BAR2)
tam_triangulo=1.2 #tamaño del triangulo en el MAPA (en km) que representa a la estacion
pos_triag_x=0.5 #posicion del triangulo de LEYENDA (x e y estan noralizados entre 0 y 1)
pos_triag_y=0.04 #
pos_tx_est_x=0.5 # #posicion del texto de la LEYENDA que dice ESTACION (x e y estan noralizados entre 0 y 1)
pos_tx_est_y=0.07 #
tam_triang_leyen=40 #Tamaño del triangulo en la LEYENDA

#--------PARAMETROS FALLAS DISPONIBLES DE REMOSION EN MASA
#---no son las del WMS del mapa geologico colombiano sino un shape interno que esta en la ruta: "/home/mlizarazo/Documentos/fallas_activas_col/FALLAS"
#---estas fallas se usan pues hay un convenio sobre las mas peligrosas, tener cuidado pues el label (nombre de la falla)fue diseñado con este codigo, alta probabilidad de fallar
ancho_fallas=0.5 #ancho del trazo de la falla
tam_tex_falla=11 #tamaño del texto con el nombre de la falla
en_negrilla=False #pone en negrilla el texto de la falla

#--------HACE PERFILES ...editar cuando se requira. Activar con True incluir_perfiles. 
# perfiles (lon_W,lon_E,lat_W,lat_E,"letra para denominar el perfil")
incluir_perfil=False
p1=(-74.35,-73.45,7.365,7.365,"A")
p2=(-74.35,-73.45,7.257,7.257,"B")
perfiles=[p1,p2] #agrgar a la lista los perfiles que se necesiten [p1,p2,p3,....pn], previamente definirlos arriba
ancho=0.08 #ancho del perfil en grados


#------------------------------------------------------------------


#--------***** FUNCIONAMIENTO DEL CODIGO ****------------
#--- Editar solo para modificar funcionamiento del codigo
#------------

if tipo_busqueda=="p":
	
	lat_min,lat_max, lon_min, lon_max,poligono=busqueda_poligono(archivo_poligono)
	print("Leyendo en Poligono: "+archivo_poligono)
	'''
	El archivo no debe tener header, cada linea es una cordenada lon lat (separados por espacio o tabulacion)
	-75.15471017	3.160705975
	-74.08543988	4.020286263
	-73.47272867	3.508702147
'''
#------------
if tipo_busqueda!="p":
	archivo_poligono=""

#-------------------Obtiene reporte bn sea de la base de datos de la rsnc o de un reporte propio_ (si es propio ver las siguientes lineas y editar arriba la opcion op = "ap")

if op == "ap": #sirve para graficar sismos que estan en una tabla, por ejemplo los sismos destacados en el turno
	#el archivo propio debe tener la siguiente informacion minima por columnas 1. fecha, 2. hora 3. lat 4. lon 5. prof 6. mag
	#si se quiere en formato estandar agregar columanas: 7. EPICENTRO, 8. RMS, 9. GAP,  10. ERROR-LAT, 11. ERROR-LON, 12. ERROR-Z
	#debe tener encabezado (la primera fila) en formato csv (separado por comas, decimales con punto)
	
	#formato: 	#fecha: aaaa-mm-dd  #hora:hh:mm:ss 	#lat y lon (separador decimal punto)
	
	
	archivo="/opt/rutinas/mapas_histogramas/archivos_propios/linea_base_local_platero.csv"
	

	
	r=selec_catalogo_propio(archivo,time_i,time_f,lat_min,lat_max,lon_min,lon_max,lon_centro,lat_centro,radio,
							prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono)


if op == "bd": #consulta las bases de datos de seissan y sesicomp

	r=report_catalogo(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,
					mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,prof_cinco_superficial,hora_local,archivo_poligono,incluir_No_loc)




#------------------------------si se desea adicional algo al mapa editar

#------ si se quieren agregar cosas al mapa, como estaciones, pozos o demas, agregar lineas de este tipo


def adiciones_mapa(mapa,ax):
	'''

	#------------ejemplo1------ adiciona dos estaciones
	lats_a=[4.4798,4.9053]
	lons_a=[-75.6317, -75.2867]
	nom_a=["PBLC","PIRM"]

	x_a,y_a=mapa(lons_a, lats_a)
	ax.text(x_a[0],y_a[0]+1000, nom_a[0], fontsize=35, horizontalalignment='center', color="k", style='italic')
	ax.text(x_a[1],y_a[1]-3000, nom_a[1], fontsize=35, horizontalalignment='center', color="k", style='italic')
	est = mapa.scatter(x_a, y_a, marker='^',color='k', edgecolor='k', s=1500, zorder=10) 
	#------------fin ejemplo1------ adiciona dos estaciones
	'''

	'''
	#------------ejemplo2------ adiciona un poligono 
	archivo_poligono="/home/mlizarazo/Documentos/codigos_python/consulta_modelo_y_tipo_mag/poligonos_zonas/el_paso.txt"

	arc_pol=open(archivo_poligono,"r")
	lineas_pol=arc_pol.readlines()
	lons_pol=[]
	lats_pol=[]
	for lc in lineas_pol:
		vec_lc=lc.split()
		lons_pol.append(float(vec_lc[0]))
		lats_pol.append(float(vec_lc[1]))
	x_c,y_c=mapa(lons_pol, lats_pol)
	cir = mapa.plot(x_c, y_c, 'r', linewidth=3)
	#--------fin ejemplo2------ adiciona un poligono  
	'''
	return mapa,ax	


#***+++++----- codigo necesario para la genreacion de mapas E HISTOGRAMAS ****----

#---- genera histogramas Y KML

if generar_kml==True:
    mapa_kml=r2kml(r)
    
if generar_histrogramas==True:
    histo_sismos(r,time_i,time_f,op,lat_centro,lon_centro,tipo_busqueda) #cambiar false por True, si se quiere el histograma de magnitudes en los meses analizados

#---- genera MAPAS


bool_mapas_deseados=[generar_mapa,generar_map_gap,generar_map_rms,generar_map_time,generar_map_calor] #es una lista de boleanos con true o false, dependiendo de que mapa se quiere hacer
mapas_posibles=["SISMICIDAD","GAP","RMS","TIEMPO","CALOR"]  

if True in bool_mapas_deseados:

	#guarda las coordenadas iniciales de busqueda, pues al generar el mapa base se modifican en la funcion mapa crudo, las ira modificando para adaptar el tamaño del mapa a su contenido
	#con lo cual, la segunda vez que se llame a la funcion mapa_base o mapa_calor (que a su vez llama la fucnion mapa_crudo) es bueno meterle las coordendas originales
	#no hay problema con coordendas circulares, pues dentro de la funcion se convierte ese cirtculo en un rectangulo haciendo modificaciones uncamente en esas variables    
	lat_min_o=lat_min
	lat_max_o=lat_max
	lon_min_o=lon_min
	lon_max_o=lon_max
	#---------------------

	#-----------genera mapa crudo con adiciones, guarda los objetos mapa y ax en memoria local y genera los mapas

	#lat_min,lat_max,lon_min,lon_max=restituir_geografia(lat_min_o,lat_max_o,lon_min_o,lon_max_o)
	mapa_c, ax_c,lat_min_c,lat_max_c,lon_min_c,lon_max_c,long_barra,x_map_c,y_map_c=mapa_crudo(r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,time_i,
																								time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,
																								zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,
																								incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,
																								FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,
																								indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,
																								tam_map_ref,color_prof)
	    
	if graficar_estaciones==True:
	    grafica_estaciones(lat_min_c,lat_max_c,lon_min_c,lon_max_c,lat_centro,lon_centro, radio,tipo_busqueda, mapa_c,ax_c,TX,TY,
	    					tam_text_est,tam_text,tam_triangulo,pos_triag_x,pos_triag_y,pos_tx_est_x,pos_tx_est_y,tam_triang_leyen,ver_temporales,
	    					ver_permanentes,imprimir_cod_estacion)

	mapa_c,ax_c=adiciones_mapa(mapa_c,ax_c)

	if tipo_busqueda=="p":
	    mapa_c,ax_c=adicion_contorno_poligono(mapa_c,ax_c,archivo_poligono)

	#------------guarda ax y map en binarios para ser usados en la recurrencia
	nombre_archivo_salida="mapa_c"
	with open(nombre_archivo_salida, "wb") as fp:
	    pickle.dump(mapa_c, fp)     

	nombre_archivo_salida="ax_c"
	with open(nombre_archivo_salida, "wb") as fp:
	    pickle.dump(ax_c, fp) 
	#---------fin-guarda ax y map en binarios para ser usados en la recurrencia

	if incluir_mapa_ref==True:
	    mapa_referencia(tipo_mapa,tam_map_ref,tipo_busqueda,radio,lon_centro,lat_centro,lat_max,lat_min,lon_max,lon_min,op,r) 


	generacion_mapas(bool_mapas_deseados,mapas_posibles,perfiles,incluir_perfil,long_barra,x_map_c,y_map_c,r,lon_centro,lat_centro,
		lat_min_o,lat_max_o,lon_min_o,lon_max_o,prof_max,prof_min,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,
		zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,
		FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,
		color_prof,ancho,lat_min_c,lat_max_c,lon_min_c,lon_max_c,alpha,exp_mag,d_mag,y_mag,lat_pos_prof,cmap)

print(f"\n¡LISTOnes!\nLos restultados del analisis se escribieron en: {os.getcwd()}\n\n")