'''
Se trajo el archivo seisan_1993_2018_02_28_web_edit.csv desde windows (separador ;). 
Este archivo corresponde a la base de datos de seisan que esta en la pagina web, se desacargo desde 1993-06-01 a 2018-02-28 (el 28/03/2022) 
EN total se tienen 174303 eventos (coincidentes con los de la pag web). Se suprimio la columna revisado y se organizo en las ultimas columnas minicipio y departamento
el archivo tiene los siguientes campos

FECHA	HORA_UTC	LATITUD (grados)	LONGITUD (grados)	PROFUNDIDAD (Km)	MAGNITUD Ml	MAGNITUD Mw	RMS (Seg)	GAP (grados)	ERROR LATITUD (Km)	ERROR LONGITUD (Km)	ERROR PROFUNDIDAD (Km)	MUNICIPIO	DEPARTAMENTO	# FASES

esta version tiene editado el municipio y el departamento, manualmente se pusieron eñes y tildes a las palabras mas comunes (seguramente faltan)

Este codigo se creo porque hay numerosos sismos sin los campos de rms, gap, o errores (las celdas estan vacias), asi que el codigo reescribe
en esas celdas un 999, que despues le indicara a quien interese que esos campos estan dañados. Los eventos con celdas vacias mayaritariamente son 
de 1994 a 1996, que tienen dañado el sfile y no cuentan con linea tipo E (error). Despues de asiganar 999, genera un nuevo archivo
que se llama 1993-2018_con_erores_lat_y_lon.csv delimitado por comas. Este archivo, tiene este nombre ya que asi esta en todos los codigos que lo leen. 
Si se cambia este nombre alterara numerosos codigos generando error
'''

import csv
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime #importante para el condicional de tiempo en la fuincion selec_catalogo_seisan
import MySQLdb
from operator import itemgetter #tiene la funcion sort que sevira para clasificar los sismos de acurdo a sus profundidades o magnitudes
from io import open
import copy

import shapely.affinity
from shapely.geometry import Point, Polygon
from descartes import PolygonPatch

from functools import partial
import pyproj
from shapely.ops import transform


arc_e=open("/home/mlizarazo/Documentos/codigos_python/mapas_histogramas/leyends_base_seisan/seisan_1993_2018_02_28_web_edit.csv","r")
arc_s=open("/home/mlizarazo/Documentos/codigos_python/mapas_histogramas/leyends_base_seisan/1993-2018_con_erores_lat_y_lon.csv","w")
arc_d=open("/home/mlizarazo/Documentos/codigos_python/mapas_histogramas/leyends_base_seisan/seisan_averiados.csv","w")


filas_e=arc_e.readlines()[1:]

arc_s.writelines("FECHA,HORA,LAT,LON,PROF,MAG,MW,RMS,GAP,ELAT,ELON,EPROF,MUNICIPIO,DEPARTAMENTO\n")
cont=0
for fe in filas_e:
	vec_sis =fe.split(';') #convierte cada sismo en un vector 
	#print(cont,vec_sis)

	fecha=vec_sis[0]
	hora=vec_sis[1]
	fecha_hora=fecha+" "+hora
	utc_fecha=UTCDateTime(fecha_hora)				
	
	lat=float(vec_sis[2])
	lon=float(vec_sis[3])
	prof=float(vec_sis[4])
	mag=float(vec_sis[5]) #ojo esta es ml, y vec_sis[6] es mw, 
	magW=vec_sis[6]
	municipio=vec_sis[12]
	departamento=vec_sis[13]
	epicentro=municipio+" - "+departamento
	
	if vec_sis[7] == "" or vec_sis[8] == "" or vec_sis[9] == "" or vec_sis[10] == "" or vec_sis[11] == "":
		arc_d.writelines("\n"+str(fecha)+","+str(hora)+","+str(lat)+","+str(lon)+","+str(prof)+","+str(mag)+","+str(magW))

	if vec_sis[7] == "" or float(vec_sis[7]) == 0 or float(vec_sis[7]) == 99.9 or float(vec_sis[7]) == 99:
		print("sin rms",fe)
		vec_sis[7]=999
		arc_d.writelines(",sin_rms")

		
	rms=float(vec_sis[7])

	if vec_sis[8] == "" or float(vec_sis[8]) == 0 or float(vec_sis[8]) == 99.9: #hay mucho sismo con 99.0 que parece ser real, asi que no se tomo 
		print("sin gap",fe)
		vec_sis[8]=999
		arc_d.writelines(",sin_gap")
		
	gap=float(vec_sis[8])

	if vec_sis[9] == "" or float(vec_sis[9]) == 0 or float(vec_sis[9]) == 99.9 or float(vec_sis[9]) == 99:
		print("sin elat",fe)
		vec_sis[9]=999
		arc_d.writelines(",sin_elat")
		
	elat=float(vec_sis[9])
	
	if vec_sis[10] == "" or float(vec_sis[10]) == 0 or float(vec_sis[10]) == 99.9 or float(vec_sis[10]) == 99:
		print("sin elon",fe)
		vec_sis[10]=999
		arc_d.writelines(",sin_elon")
		
	elon=float(vec_sis[10])

	if vec_sis[11] == "" or float(vec_sis[11]) == 0 or float(vec_sis[11]) == 99.9 or float(vec_sis[11]) == 99:
		print("sin eprof",fe)
		vec_sis[11]=999
		arc_d.writelines(",sin_eprof")
		
	eprof=float(vec_sis[11])
	
	
	arc_s.writelines(str(fecha)+","+str(hora)+","+str(lat)+","+str(lon)+","+str(prof)+","+str(mag)+","+str(magW)+","+str(rms)+","+str(gap)+","+str(elat)+","+str(elon)+","+str(eprof)+","+str(municipio)+","+str(departamento)+"\n")
	cont=cont+1

arc_e.close()
arc_s.close()
arc_d.close()

'''


def selec_catalogo_seisan(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda):
		
	print("Leyendo base de datos Seisan...")
	
	if tipo_busqueda=="c":
		
		circulo_p=genera_circulo(radio,lon_centro,lat_centro) #circulo de interes para intersecar eventos

	archivo="/home/mlizarazo/Documentos/codigos_python/mapas_histogramas/leyends_base_seisan/base_seisan_1993_06_01__2018_02_28_epic_editado_2.csv" #archivo obtenido de la base de datos de seisaan DEBE ESTAR EN LA CARPETA DESDE LA QUE SE CORRE EL CODIGO
	#el archivo contiene en sus columnas: 0.an,mes,dia 1.hora utc 2.lat 3.lon 4.prof 5.ml 6.mw 7.rms 8.gap 9.error lat 10.error lon 11.error prof
	utc_time_i = UTCDateTime(time_i)
	utc_time_f = UTCDateTime(time_f)

	with open(archivo, 'r') as catalogo:
		catalogo = catalogo.readlines() #lee cada fila del catalogo, cada fila contiene un sismo

		sismos=catalogo[1:] 				
		
		r_seisan=[] #contiene la matriz con los sismos	
		cont=0
		for i in sismos:

			
			vec_sis =i.split(',') #convierte cada sismo en un vector (ya que viene separado por comas) #sismos[0] es el encabezado, desde sismos[1] inicia la lectura de los sismos
			#print(cont,vec_sis)

			fecha=vec_sis[0]
			hora=vec_sis[1]
			fecha_hora=fecha+" "+hora
			utc_fecha=UTCDateTime(fecha_hora)				
			
			lat=float(vec_sis[2])
			lon=float(vec_sis[3])
			prof=float(vec_sis[4])
			mag=float(vec_sis[5]) #ojo esta es ml, y vec_sis[6] es mw, 
			magW=vec_sis[6]

			if len (magW)>0:
				if magW[0].isalnum() == True: 
					mag=float(magW)
			

			rms=float(vec_sis[7])
			gap=float(vec_sis[8])
			elat=float(vec_sis[9])
			elon=float(vec_sis[10])
			eprof=float(vec_sis[11])
			
			municipio=vec_sis[12]
			departamento=vec_sis[13][:-1]
			epicentro=municipio+" - "+departamento
				
			#Fecha  hora  lat  lon  Z  Ml  Mw  RMS  GAP  e_lat  e_lon  e_prof epicentro
			if utc_time_i<=utc_fecha<=utc_time_f:
				
				if tipo_busqueda=="r":

					if lat_min<=lat<=lat_max and lon_min<=lon<=lon_max and prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:
						vec_sis_seisan=[fecha,hora,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]

						r_seisan.append(vec_sis_seisan) #escribe una matriz con los sismos seleccionados
				
				if tipo_busqueda=="c":
					
					evento=Point(lon, lat)
					
					if evento.within(circulo_p)==True:				

						if  prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:
							vec_sis_seisan=[fecha,hora,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]
							r_seisan.append(vec_sis_seisan) #escribe una matriz con los sismos seleccionados

			cont=cont+1
	return r_seisan #r contiene la matriz con los sismos filtrados


#-----------CONDICIONES DE BUSQUEDA (editar coordenadas a las necesarias y elegir bd (base de datos) o ap (archivo propio)
#---Eleccion de base de datos, si el archivo de sismos es propio op="ap", si requiere consultar bd de seisan y seiscomp  op="bd". Si op="ap" ir a linea ~125 para dar la ruta del archivo propio y leer las instrucciones 
op ="bd"#	
#op ="ap"	

#-----------------Elegir tipo de busqueda
tipo_busqueda="r" #Elegir c: circular, r: rectangular, "sin_filtro_loc"
	#si se elige opcion sin_filtro_loc, esta solo servira cuando la base de datos es propia op = "ap". En este caso, obtiene lat, lon y prof minimas y maximas del propio catalogo, con el fin de que se analice todo el contenido 
	#del catalogo sin imponer restricciones espaciales, sin embargo, si tiene habilitado el filtro por minimos y maximos de tiempo, mag y errores (funciona unicamente para base de datos propia)
		
#------Editar para busquedas rectangulares
lat_min=-7
lat_max=17
lon_min=-89
lon_max=-65


#------Editar para busquedas circulares

lat_centro=7.84
lon_centro=-74.3
radio=30 #en kilometros
#kale lat_centro=7.3655, lon_centro=-73.8566
#platero lat_centro=7.2571, lon_centro=-73.8938
#-----------------------------

#------------PARAMETROS DE BUSQUEDA 

time_i="1993-06-01 05:00:00" #fecha UT
time_f="2018-02-28 23:59:59" #fecha UT 

prof_min=-100 #para que tome sismos aereos
prof_max=1000
mag_min=0
mag_max=10
elat_max=1000
elon_max=1000
eprof_max=1000
#----------------------------
r_seisan = selec_catalogo_seisan(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda)

for rs in r_seisan:
	print(rs)

'''