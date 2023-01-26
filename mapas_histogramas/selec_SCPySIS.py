# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:39:46 2018

@author: mlizarazo
La funcion selec_catalogo_seisan(time_i,time_f,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max) 
exporta el archivo "1993-2018_con_erores_lat_y_lon.csv que contiene todos los sismos en la base de datos de seisan desde 1993 a 2018 
La funcion selec_catalogo_seiscomp(time_i,time_f,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max)
hace una consulta a la base de datos MySQLdb host="10.100.100.232",    # your host, usually localhost
	                     					user="consulta",         # your username
	                     					passwd="consulta",  # your password
	                     					db="seiscomp3")
de la cual se extraen los sismos de SEISCOMP del 01-03-2018 hasta la fecha final que decida el usuario
El codigo genera un unico report.cvs que se guarda en la carpeta donde se corre el codigo. 
Esta configurado para hacer una seleccion por tiempo, lat, lon, prof y sus errores
se deben tener los paquetes de python: Basemap, obspy y MySQLdb intalatados y el codigo funcion_mapa del cual se exportara la funcion mapa_base

La funcion report_catalogos debe ser usada en otros codigos como: from selec_SCPySIS import report_catalogos
debe estar grabado el archivo selec_SCPySIS.py en la misma carpeta
requiede los siguientes parametros de entrada report_catalogos(time_i,time_f,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max)
devuelve r, que es una matriz que contiene los sismos escritos asi: Fecha  hora  lat  lon  Z  Ml  Mw  RMS  GAP  e_lat  e_lon  e_prof
para usar la funcion:

r_SCPySIS=report_catalogos(time_i,time_f,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max)
"""
import csv
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime 
from datetime import date, datetime, timedelta
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

#dos tipos de busqueda r o c:

#if tipo_busqueda=="r":

def genera_circulo(radio,lon_centro,lat_centro): 
    '''
    Genera una superficie circular geografica con un radio (km) centrado en unas cordenadas dadas.
    Para lograr una mayor precision se usa una proyección AEQD centrada en el pozo, que proyecta 
    distancias iguales en todas las direcciones (la proyección de Mercator genera mayor distorsion) 
    La proyección equidistante azimutal es una proyección cartográfica azimutal. Tiene las útiles propiedades de que todos los puntos en el mapa están 
    a distancias proporcionalmente correctas desde el punto central, y que todos los puntos en el mapa están en el acimut (dirección) correcto desde 
    el punto central.
    '''      
    lon_p=lon_centro
    lat_p=lat_centro
    r=radio
    
    
    proj_wgs84 = pyproj.Proj('+proj=longlat +datum=WGS84')
    proj_aeqd = pyproj.Proj('+proj=aeqd +lat_0='+str(lat_p)+' +lon_0='+str(lon_p)+' +x_0=0 +y_0=0')
    
    project = partial(pyproj.transform,proj_aeqd,proj_wgs84)
        
    buf_cir = Point(0, 0).buffer(r* 1000)
       
    circulo_p=transform(project, buf_cir)    
    
    return circulo_p


#if hora_local==True:
def convierte_hora_local(r):
				
	r_HL=[]
	for i in r:
		fecha=str(i[0])
		horas=str(i[1])
		HL_fecha=str(UTCDateTime(fecha+horas).datetime -timedelta(hours=5)) #para hora local restar 5 horas
		vec_fecha=HL_fecha.split(" ")
		an_me_di=vec_fecha[0]
		hora_HL=vec_fecha[1]

		r_HL.append([an_me_di,hora_HL]+i[2:])

	r=r_HL

	return r
	
def busqueda_poligono(archivo_poligono):
	#archivo_poligono="/home/mlizarazo/Documentos/codigos_python/seiscomp2phasedat/poligonos/poli_alg_norte.txt"
	arc=open(archivo_poligono,"r")
	lineas=arc.readlines()    
	coords=[]
	lons_poli=[]
	lats_poli=[]
	for l in lineas:
		
		vec_l=l.split()
		coords.append((float(vec_l[0]), float(vec_l[1])))
		lons_poli.append(float(vec_l[0]))
		lats_poli.append(float(vec_l[1]))
	poligono = Polygon(coords)

	#define coordendas de busqueda rectangulares traves de los maximos y minimos del poligno (a esa busqueda posteriormente se le filtrara en el polgno como tal)
	lat_min=min(lats_poli)
	lat_max=max(lats_poli)
	lon_min=min(lons_poli)
	lon_max=max(lons_poli)

	return lat_min,	lat_max, lon_min, lon_max,poligono

def selec_catalogo_seisan(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono):
		
	print("Leyendo base de datos Seisan...")
	
	if tipo_busqueda=="c":
		
		circulo_p=genera_circulo(radio,lon_centro,lat_centro) #circulo de interes para intersecar eventos

	if tipo_busqueda=="p":

		lat_min,lat_max, lon_min, lon_max,poligono=busqueda_poligono(archivo_poligono)

	archivo="/opt/rutinas/mapas_histogramas/catSeis_ests_leyens/1993-2018_con_erores_lat_y_lon.csv" #archivo obtenido de la base de datos de seisaan DEBE ESTAR EN LA CARPETA DESDE LA QUE SE CORRE EL CODIGO
	#el archivo contiene en sus columnas: 0.an,mes,dia 1.hora utc 2.lat 3.lon 4.prof 5.ml 6.mw 7.rms 8.gap 9.error lat 10.error lon 11.error prof
	utc_time_i = UTCDateTime(time_i)
	utc_time_f = UTCDateTime(time_f)

	with open(archivo, 'r') as catalogo:
		catalogo = catalogo.readlines() #lee cada fila del catalogo, cada fila contiene un sismo

		sismos=catalogo[1:] 				
		
		r_seisan=[] #contiene la matriz con los sismos	
		for i in sismos:

			vec_sis =i.split(',') #convierte cada sismo en un vector (ya que viene separado por comas) #sismos[0] es el encabezado, desde sismos[1] inicia la lectura de los sismos
			
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
			departamento=vec_sis[13].replace("\n","")
			epicentro=municipio+" - "+departamento
				
			#Fecha  hora  lat  lon  Z  Ml  Mw  RMS  GAP  e_lat  e_lon  e_prof epicentro
			if utc_time_i<=utc_fecha<=utc_time_f:
				
				if tipo_busqueda=="r":

					if lat_min<=lat<=lat_max and lon_min<=lon<=lon_max and prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:
						vec_sis_seisan=[fecha,hora,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]

						r_seisan.append(vec_sis_seisan) #escribe una matriz con los sismos seleccionados

				if tipo_busqueda=="p":

					p = Point(lon, lat)

					if p.within(poligono)==True:

						if  prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:
							vec_sis_seisan=[fecha,hora,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]
							r_seisan.append(vec_sis_seisan) #escribe una matriz con los sismos seleccionados
				
				if tipo_busqueda=="c":
					
					evento=Point(lon, lat)
					
					if evento.within(circulo_p)==True:				

						if  prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:
							vec_sis_seisan=[fecha,hora,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]
							r_seisan.append(vec_sis_seisan) #escribe una matriz con los sismos seleccionados

	
	return r_seisan #r contiene la matriz con los sismos filtrados

def selec_catalogo_seiscomp(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono,incluir_No_loc):

	print("Leyendo base de datos seiscomp3...")
	
	if tipo_busqueda=="c":
		
		circulo_p=genera_circulo(radio,lon_centro,lat_centro)

		cords=circulo_p.exterior.coords[:]
		lats_c=[]
		lons_c=[]

		for c in cords:
		    lons_c.append(c[0])
		    lats_c.append(c[1])
		
		#modifica los minimos que ingresan a la consulta que por defecto es rectangular   
		lat_min=min(lats_c)-0.1
		lat_max=max(lats_c)+0.1
		lon_min=min(lons_c)-0.1
		lon_max=max(lons_c)+0.1
	
	if tipo_busqueda=="p":

		lat_min,lat_max, lon_min, lon_max,poligono=busqueda_poligono(archivo_poligono)		

	#db = MySQLdb.connect(host="10.100.100.232",user="consulta", passwd="consulta", db="seiscomp3")
	db = MySQLdb.connect(host="172.25.3.135", user="consulta", passwd="consulta", db="seiscomp3")       
	cur = db.cursor()

	if incluir_No_loc==False:
		
		cur.execute("Select \
		POEv.publicID, \
		Origin.time_value as 'Fecha', \
		ROUND(Origin.latitude_value,6) as 'Lat', \
		ROUND(Origin.longitude_value,6) as 'Lon', \
		ROUND(Origin.depth_value,4) as 'Z', \
		ROUND(Magnitude.magnitude_value,1) as 'Mag', \
		Magnitude.type ,\
		ROUND(Origin.quality_standardError,1) as 'RMS', \
		CONCAT(ROUND(Origin.quality_azimuthalGap,1)) as 'GAP', \
		quality_usedStationCount as 'NF', \
		convert(cast(convert(EventDescription.text using latin1) as binary) using utf8) as 'UBICACION', \
		Origin.latitude_uncertainty, \
		Origin.longitude_uncertainty, \
		Origin.depth_uncertainty, \
		Origin.creationInfo_author \
		from \
		Event AS EvMF \
		left join PublicObject AS POEv ON EvMF._oid = POEv._oid \
		left join PublicObject as POOri ON EvMF.preferredOriginID=POOri.publicID \
		left join Origin ON POOri._oid=Origin._oid \
		left join PublicObject as POMag on EvMF.preferredMagnitudeID=POMag.publicID \
		left join Magnitude ON Magnitude._oid = POMag._oid \
		left join EventDescription \
		ON EvMF._oid=EventDescription._parent_oid \
		where \
		(Magnitude.magnitude_value between "+str(mag_min)+" and "+str(mag_max)+" and \
		Origin.latitude_value between "+str(lat_min)+" and "+str(lat_max)+" and \
		Origin.longitude_value between "+str(lon_min)+" and "+str(lon_max)+" and \
		Origin.depth_value between "+str(prof_min)+" and "+str(prof_max)+" and \
		Origin.latitude_uncertainty<"+str(elat_max)+" and \
		Origin.longitude_uncertainty<"+str(elon_max)+" and \
		Origin.depth_uncertainty<"+str(eprof_max)+" and \
		Origin.evaluationMode like 'manual' and \
		Origin.time_value between '"+time_i+"' and '"+time_f+"')\
		and (EvMF.type like 'earthquake' or EvMF.type like 'volcanic eruption' or EvMF.type is null )") #ultima linea para filtrar los no localizables (para casos como el de providencia quitar esta linea, pues interesan los no localizables)

	if incluir_No_loc==True:
		
		print("Incluyendo NO localizables")

		cur.execute("Select \
		POEv.publicID, \
		Origin.time_value as 'Fecha', \
		ROUND(Origin.latitude_value,6) as 'Lat', \
		ROUND(Origin.longitude_value,6) as 'Lon', \
		ROUND(Origin.depth_value,4) as 'Z', \
		ROUND(Magnitude.magnitude_value,1) as 'Mag', \
		Magnitude.type ,\
		ROUND(Origin.quality_standardError,1) as 'RMS', \
		CONCAT(ROUND(Origin.quality_azimuthalGap,1)) as 'GAP', \
		quality_usedStationCount as 'NF', \
		convert(cast(convert(EventDescription.text using latin1) as binary) using utf8) as 'UBICACION', \
		Origin.latitude_uncertainty, \
		Origin.longitude_uncertainty, \
		Origin.depth_uncertainty, \
		Origin.creationInfo_author \
		from \
		Event AS EvMF \
		left join PublicObject AS POEv ON EvMF._oid = POEv._oid \
		left join PublicObject as POOri ON EvMF.preferredOriginID=POOri.publicID \
		left join Origin ON POOri._oid=Origin._oid \
		left join PublicObject as POMag on EvMF.preferredMagnitudeID=POMag.publicID \
		left join Magnitude ON Magnitude._oid = POMag._oid \
		left join EventDescription \
		ON EvMF._oid=EventDescription._parent_oid \
		where \
		(Magnitude.magnitude_value between "+str(mag_min)+" and "+str(mag_max)+" and \
		Origin.latitude_value between "+str(lat_min)+" and "+str(lat_max)+" and \
		Origin.longitude_value between "+str(lon_min)+" and "+str(lon_max)+" and \
		Origin.depth_value between "+str(prof_min)+" and "+str(prof_max)+" and \
		Origin.latitude_uncertainty<"+str(elat_max)+" and \
		Origin.longitude_uncertainty<"+str(elon_max)+" and \
		Origin.depth_uncertainty<"+str(eprof_max)+" and \
		Origin.evaluationMode like 'manual' and \
		Origin.time_value between '"+time_i+"' and '"+time_f+"')") 

	#0.id sismo 1.datatime.datatime (año ...seg) 2.lat 3.lon. 4.prof 5.mag 6.tipo magnitud 7.rms 8.gap 9.N fses 10.ubicacion 11.elat, 12.elon, 13.eprof
	
	r_seiscomp=[] #contiene la matriz con los sismos	seleccionados
	for i in cur.fetchall():
		
		Ids=i[0]
		fecha=str(i[1])
		vec_fecha=fecha.split(" ")

		an_me_di=vec_fecha[0]
		hora_utc=vec_fecha[1]
		
		lat=float(i[2])
		lon=float(i[3])
		prof=float(i[4])
		mag=float(i[5]) #ojo esta es m indistinta si es ml mlr mb mw (es con la que se fijo el evento)
		
		rms=float(i[7])
		gap=float(i[8])
		elat=float(i[11])
		elon=float(i[12])
		eprof=float(i[13])


		autor=i[14]
		'''
		if autor=="ptogaitan@proc1" or autor=="estella@proc1":
			print(fecha,autor)
		'''

		epicentro=i[10].split(",")[0].upper()

		if tipo_busqueda=="r":

			vec_sis_sesicomp=[an_me_di,hora_utc,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]									
			r_seiscomp.append(vec_sis_sesicomp) #escribe una matriz con los sismos seleccionado

		if tipo_busqueda=="p":

			p = Point(lon, lat)

			if p.within(poligono)==True:

				vec_sis_sesicomp=[an_me_di,hora_utc,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]
				r_seiscomp.append(vec_sis_sesicomp) #escribe una matriz con los sismos seleccionado

		if tipo_busqueda=="c":
			
			evento=Point(lon, lat)

			if evento.within(circulo_p)==True:			
														
				vec_sis_sesicomp=[an_me_di,hora_utc,lat,lon,prof,mag,rms,gap,elat,elon,eprof,epicentro]									
				r_seiscomp.append(vec_sis_sesicomp) #escribe una matriz con los sismos seleccionados
				

	cur.close()
	db.close()
	
	
	return r_seiscomp


def report_catalogo(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,prof_cinco_superficial,hora_local,archivo_poligono,incluir_No_loc):
	
	if incluir_No_loc==True:
		if UTCDateTime(time_i) < UTCDateTime("2018-03-01"):
			incluir_No_loc=False
			print("No se incluiran los No localizables... Intente esta opcion solo despues de 2018-03-01")

	print("Generando reporte...")
	#busca en los catalogos de seisan y seiscomp y unifica el report

	#-------------------------
	
	if UTCDateTime(time_i) >= UTCDateTime("2018-03-01"):
		r_seisan=[]
		r_seiscomp = selec_catalogo_seiscomp(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono,incluir_No_loc)
			
	if UTCDateTime(time_i) < UTCDateTime("2018-03-01") and UTCDateTime(time_f) < UTCDateTime("2018-03-01"):
		r_seisan = selec_catalogo_seisan(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono)
		r_seiscomp = []

	if UTCDateTime(time_i) < UTCDateTime("2018-03-01") and UTCDateTime(time_f) >= UTCDateTime("2018-03-01"):
		r_seisan = selec_catalogo_seisan(time_i,time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono)
		r_seiscomp = selec_catalogo_seiscomp("2018-03-01",time_f,lat_centro,lon_centro,radio,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono,incluir_No_loc)
	

	#Se obliga a Seiscomp a iniciar la busqueda desde 2018-03-01 en caso que el usauario este buscando desde un tiempo inicial anterior a esta fecha ya que por la Migracion SEiscomp ya tiene en su base de datos 
	#sismos anteriores a esa fecha por tanto repite la busqueda que ya se hizo en Seisan y se repiten eventos. Sin embargo, si el usuario, inicia la busqueda porsterior a esta fecha, no se puede fijar "2018-03-01"
	#como tiempo inicial, sino que ahora si la funcion recibe el time_i dado por el usuario
	r_SCPySIS=r_seisan+r_seiscomp
	
	#-------------------------------
	#escribe un archivo reporte.csv 

	if hora_local==True:
		r_SCPySIS=convierte_hora_local(r_SCPySIS)
	
	reporte_salida=r_SCPySIS
	archivo_salida = open('reporte.csv', 'w') #escribe archivo 
	archivo_salida_ZMAP = open('reporte_f_ZMAP.csv', 'w') #escribe archivo 
	writer = csv.writer(archivo_salida)
	writer_ZMAP = csv.writer(archivo_salida_ZMAP)

	#---header
	if hora_local==True:
		writer.writerow(["FECHA ", "HORA LOCAL", "LATITUD º", "LONGITUD º", "PROFUNDIDAD (km)", "MAGNITUD", "EPICENTRO","RMS", "GAP", "ERROR-LAT", "ERROR-LON", "ERROR-Z"])

	else:
		writer.writerow(["FECHA ", "HORA UTC", "LATITUD º", "LONGITUD º", "PROFUNDIDAD (km)", "MAGNITUD", "EPICENTRO","RMS", "GAP", "ERROR-LAT", "ERROR-LON", "ERROR-Z"])
	#---

	for s in reporte_salida:

		prof_rep=round(s[4], 3)
		
		if prof_cinco_superficial==True:
			#para solicitudes se requiere que para Z<5 escriba superficial, por ello se activa el if, para que escriba el numero comentar el if (con la linea que engloba el if)
			if prof_rep <= 5: 
				prof_rep="superficial"			 	
		
		writer.writerow([s[0],s[1],round(s[2], 3),round(s[3], 3),prof_rep,s[5],s[11],s[6],s[7],s[8],s[9],s[10]])

		#escribe reporte en formato ZMAP
		#Lon lat year Month day magnitude Depth [km] Hour Minute Second
		
		fecha=s[0]
		horas=s[1]
		an=fecha.split("-")[0]
		mes=fecha.split("-")[1]
		dia=fecha.split("-")[2]
		hora=horas.split(":")[0]
		minu=horas.split(":")[1]
		seg=horas.split(":")[2]
		writer_ZMAP.writerow([s[3],s[2],  an, mes, dia, s[5], s[4], hora, minu, seg])

	archivo_salida.close()
	archivo_salida_ZMAP.close()
	#-----------fin escribe formato zmap
		

	if tipo_busqueda=="r":

		print(f"\n'Se escribio reporte.csv' (en la ruta {os.getcwd()}) contiene {len(r_SCPySIS)} sismos\n \
			{len(r_seisan)} Sismos en SEISAN y {len(r_seiscomp)} en SeisComp. \n\n\
			Se filtraron las bases de datos con los siguientes parametros: \n\n\
			tiempo: {time_i} a {time_f} \n\
			latitud: {lat_min}º a {lat_max}º, \n\
			longitud: {lon_min}º a {lon_max}º \n\
			profundidad: {prof_min}km a {prof_max}km \n\
			magnitud: {mag_min} a {mag_max},\n \
			maximo error en latitud: {elat_max}km \n\
			maximo error en longitud {elon_max}km \n\
			maximo error en profundidad {eprof_max}km\n\n")

	if tipo_busqueda=="p":

		print(f"\n'Se escribio reporte.csv' (en la ruta {os.getcwd()}) contiene {len(r_SCPySIS)} sismos\n \
			{len(r_seisan)} Sismos en SEISAN y {len(r_seiscomp)} en SeisComp. \n\n\
			Se filtraron las bases de datos con los siguientes parametros: \n\n\
			Poligono: {archivo_poligono} \n\
			tiempo: {time_i} a {time_f} \n\
			profundidad: {prof_min}km a {prof_max}km \n\
			magnitud: {mag_min} a {mag_max},\n \
			maximo error en latitud: {elat_max}km \n\
			maximo error en longitud {elon_max}km \n\
			maximo error en profundidad {eprof_max}km\n\n")
	
	if tipo_busqueda=="c":

		print(f"\n'Se escribio reporte.csv' (en la ruta {os.getcwd()}) contiene {len(r_SCPySIS)} sismos\n \
			{len(r_seisan)} Sismos en SEISAN y {len(r_seiscomp)} en SeisComp. \n\n\
			Se filtraron las bases de datos con los siguientes parametros: \n\n\
			tiempo: {time_i} a {time_f} \n\
			latitud central: {lat_centro}º, \n\
			longitud central: {lon_centro}º \n\
			radio: {radio} km \n\
			profundidad: {prof_min}km a {prof_max} km \n\
			magnitud: {mag_min} a {mag_max},\n \
			maximo error en latitud: {elat_max} km \n\
			maximo error en longitud {elon_max} km \n\
			maximo error en profundidad {eprof_max} km\n\n")

	return r_SCPySIS



def selec_catalogo_propio(archivo,time_i,time_f,lat_min,lat_max,lon_min,lon_max,lon_centro,lat_centro,radio,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,tipo_busqueda,archivo_poligono):
	
	print("Generando reporte...")
	print("Leyendo "+archivo+"...\n")
	#si se elige base propia "ap"
	#el archivo propio debe tener la siguiente informacion minima por columnas 1. fecha, 2. hora 3. lat 4. lon 5. prof 6. mag
	#si se quiere en formato estandar agregar columanas: 7. EPICENTRO, 8. RMS, 9. GAP,  10. ERROR-LAT, 11. ERROR-LON, 12. ERROR-Z
	#No debe tener encabezado (desde la primera fila se presenta la informacion) en formato csv (separado por comas, si los decimales estan en comas tranformarlo a puntos)
	
	catalogo = open(archivo, 'r')
	catalogo = catalogo.readlines() #lee cada fila del catalogo, cada fila contiene un sismo		
	
	
	if ";" in catalogo[1] and "," not in catalogo[1]: 
		separador=";"

	if "," in catalogo[1] and ";" not in catalogo[1]: 
		separador=","

	if "," not in catalogo[1] and ";" not in catalogo[1]: 
		print("no se encontro separador por comas en el archivo (reviselo)")

	if "," in catalogo[1] and ";" in catalogo[1]:
		separador="," 
		print("se encontraron , y ; en el archivo, se tomara separador por comas")

	if tipo_busqueda=="sin_filtro_loc":
		#obtiene lat, lon y prof 
		#minimas y maximas del propio catalogo, con el fin de que se analice todo el contenido del catalogo sin imponer restricciones espaciales, 
		#sin embargo, filtra por minimos y maximos de tiempo, mag y errores (funciona unicamente para base de datos propia)
		
		sismos=catalogo[1:] 
		tam_catalogo = len(sismos[1].split(separador) )
		
		lats=[]
		lons=[]
		profs=[]
		for i in sismos:
				
			vec_sis =i.split(separador) #convierte cada sismo en un vector (ya que viene separado por comas) #sismos[0] es el encabezado, desde sismos[1] inicia la lectura de los sismos
				
			lats.append(float(vec_sis[2]))
			lons.append(float(vec_sis[3]))
			profs.append(float(vec_sis[4]))
		
		lat_min=min(lats)
		lat_max=max(lats)
		lon_min=min(lons)
		lon_max=max(lons)
		prof_min=min(profs)
		prof_max=max(profs)

	if tipo_busqueda=="c":
		
		circulo_p=genera_circulo(radio,lon_centro,lat_centro) #circulo de interes para intersecar eventos

	if tipo_busqueda=="p":

		lat_min,lat_max, lon_min, lon_max,poligono=busqueda_poligono(archivo_poligono)	

	utc_time_i = UTCDateTime(time_i)
	utc_time_f = UTCDateTime(time_f)


	with open(archivo, 'r') as catalogo:		
		
		catalogo = catalogo.readlines() #lee cada fila del catalogo, cada fila contiene un sismo		
		
		sismos=catalogo[1:] 
		tam_catalogo = len(sismos[1].split(separador) )		

		#-------------------------------
		#escribe un archivo reporte_propio.csv 
		archivo_salida = open('reporte_propio.csv', 'w') #escribe un archivo llamado simos_filtrados.csv en la carpeta actual
		writer = csv.writer(archivo_salida)
		archivo_salida_ZMAP = open('reporte_propio_ZMAP.csv', 'w') #escribe un archivo llamado simos_filtrados.csv en la carpeta actual
		writer_ZMAP = csv.writer(archivo_salida_ZMAP)
		header=catalogo[0].split(separador)
		if "\n" in header[-1]:
			header[-1]=header[-1].replace("\n","")
		writer.writerow(header)

		r_propio=[] #contiene la matriz con los sismos	
		
		if tam_catalogo >= 12: #comprueba sei el catalogo leido tiene 12 campos: #el archivo propio debe tener la siguiente informacion minima por columnas 1. fecha, 2. hora 3. lat 4. lon 5. prof 6. mag
			#si se quiere en formato estandar agregar columanas: 7. EPICENTRO, 8. RMS, 9. GAP,  10. ERROR-LAT, 11. ERROR-LON, 12. ERROR-Z
			for i in sismos:
				
				vec_sis =i.split(separador) #convierte cada sismo en un vector (ya que viene separado por comas) #sismos[0] es el encabezado, desde sismos[1] inicia la lectura de los sismos
				
				if "\n" in vec_sis[-1]:
					vec_sis[-1]=vec_sis[-1].replace("\n","")

				try:

					fecha=vec_sis[0]
					hora=vec_sis[1]
					fecha_hora=fecha+" "+hora			
					utc_fecha=UTCDateTime(fecha_hora)
					
					an=fecha.split("-")[0]
					mes=fecha.split("-")[1]
					dia=fecha.split("-")[2]
					ho=hora.split(":")[0]
					minu=hora.split(":")[1]
					seg=hora.split(":")[2]

					lat=float(vec_sis[2])
					lon=float(vec_sis[3])
					prof=float(vec_sis[4])
					mag=float(vec_sis[5])				
					
					elat=float(vec_sis[9])
					elon=float(vec_sis[10])
					eprof=float(vec_sis[11])
					
					
					if utc_time_i<=utc_fecha<=utc_time_f:
						
						if tipo_busqueda=="r" or tipo_busqueda=="sin_filtro_loc": 
						#obtiene lat, lon y prof 
						#minimas y maximas del propio catalogo, con el fin de que se analice todo el contenido del catalogo sin imponer restricciones espaciales, 
						#sin embargo, filtra por minimos y maximos de tiempo, mag y errores (funciona unicamente para base de datos propia)

							if lat_min<=lat<=lat_max and lon_min<=lon<=lon_max and prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:
								
								r_propio.append(vec_sis) #escribe una matriz con los sismos seleccionados
								writer.writerow(vec_sis)
								writer_ZMAP.writerow([lon,lat,  an, mes, dia, mag, prof, ho, minu, seg])

						if tipo_busqueda=="p":

							p = Point(lon, lat)

							if p.within(poligono)==True:

								if prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max:

									r_propio.append(vec_sis) #escribe una matriz con los sismos seleccionados
									writer.writerow(vec_sis)
									writer_ZMAP.writerow([lon,lat,  an, mes, dia, mag, prof, ho, minu, seg])


						if tipo_busqueda=="c":
								
							evento=Point(lon, lat)

							if evento.within(circulo_p)==True:

								if prof_min<=prof<=prof_max and mag_min<=mag<=mag_max and elat<elat_max and elon<elon_max and eprof<eprof_max: 
									
									r_propio.append(vec_sis) #escribe una matriz con los sismos seleccionados
									writer.writerow(vec_sis)
									writer_ZMAP.writerow([lon,lat,  an, mes, dia, mag, prof, ho, minu, seg])

				except:
						
						print("No se pudo leer la fila:\n",i)	
											
		
		elif 6 <= tam_catalogo < 12: #comprueba si el catalogo leido tiene 6 campos:  1. fecha, 2. hora 3. lat 4. lon 5. prof 6. mag
			
			for i in sismos:
				

				vec_sis =i.split(separador) #convierte cada sismo en un vector (ya que viene separado por comas) #sismos[0] es el encabezado, desde sismos[1] inicia la lectura de los sismos
				
				if "\n" in vec_sis[-1]:
					vec_sis[-1]=vec_sis[-1].replace("\n","")

				try:
					fecha=vec_sis[0]
					hora=vec_sis[1]
					fecha_hora=fecha+" "+hora
					utc_fecha=UTCDateTime(fecha_hora)

					an=fecha.split("-")[0]
					mes=fecha.split("-")[1]
					dia=fecha.split("-")[2]
					ho=hora.split(":")[0]
					minu=hora.split(":")[1]
					seg=hora.split(":")[2]

					lat=float(vec_sis[2])
					lon=float(vec_sis[3])
					prof=float(vec_sis[4])				
					mag=float(vec_sis[5])
					
										
					if utc_time_i<=utc_fecha<=utc_time_f:
					
						if tipo_busqueda=="r" or tipo_busqueda=="sin_filtro_loc":

							if lat_min<=lat<=lat_max and lon_min<=lon<=lon_max and prof_min<=prof<=prof_max and mag_min<=mag<=mag_max:
								
								r_propio.append(vec_sis) #escribe una matriz con los sismos seleccionados
								writer.writerow(vec_sis)
								writer_ZMAP.writerow([lon,lat,  an, mes, dia, mag, prof, ho, minu, seg])
													
						if tipo_busqueda=="p":

							p = Point(lon, lat)

							if p.within(poligono)==True:

								if prof_min<=prof<=prof_max and mag_min<=mag<=mag_max: 
									
									r_propio.append(vec_sis) #escribe una matriz con los sismos seleccionados
									writer.writerow(vec_sis)
									writer_ZMAP.writerow([lon,lat,  an, mes, dia, mag, prof, ho, minu, seg])

						if tipo_busqueda=="c":
								
							evento=Point(lon, lat)

							if evento.within(circulo_p)==True:
								
								if prof_min<=prof<=prof_max and mag_min<=mag<=mag_max: 
									
									r_propio.append(vec_sis) #escribe una matriz con los sismos seleccionados
									writer.writerow(vec_sis)
									writer_ZMAP.writerow([lon,lat,  an, mes, dia, mag, prof, ho, minu, seg])

				except:
						print("\nNo se pudo leer la fila:\n",i)	
									
								
		else:
			print("\nRevise el catalogo ingresado \n")
			print("\nEl catalogo ingresado debe tener encabezado (la primera fila presenta la informacion fecha, hora, lat...) \nLos campos minimos son:\n\n1. fecha, 2. hora 3. lat 4. lon 5. prof 6. mag \n\nLos campos ideales son:\n\n\
					     1. fecha, 2. hora 3. lat 4. lon 5. prof 6. mag 7. EPICENTRO, 8. RMS, 9. GAP,  10. ERROR-LAT, 11. ERROR-LON, 12. ERROR-Z\n\nSe require en formato csv (separado por comas). Seguir el formato: \n\n\
						fecha: aaaa-mm-dd \n hora:hh:mm:ss \n lat y lon (separador decimal punto)\n\n")

	print("Se escribio el reporte_propio que contiene "+str(len(r_propio))+" sismos\n")
	
		
	archivo_salida.close()
	archivo_salida_ZMAP.close()

	if tipo_busqueda=="r" or tipo_busqueda=="sin_filtro_loc":

		print(f"\n'Se escribio reporte_propio.csv' (en la ruta {os.getcwd()}) contiene {len(r_propio)} sismos\n\n \
			Se filtraron las bases de datos con los siguientes parametros: \n\n\
			tiempo: de {time_i} a {time_f} \n\
			latitud: {lat_min}º a {lat_max}º, \n\
			longitud: {lon_min}º a {lon_max}º \n\
			profundidad: {prof_min}km a {prof_max}km \n\
			magnitud: {mag_min} a {mag_max},\n \
			maximo error en latitud: {elat_max}km \n\
			maximo error en longitud {elon_max}km \n\
			maximo error en profundidad {eprof_max}km\n\n")
	
	if tipo_busqueda=="c":

		print(f"\n'Se escribio reporte_propio.csv' (en la ruta {os.getcwd()}) contiene {len(r_propio)} sismos\n\n \
			Se filtraron las bases de datos con los siguientes parametros: \n\n\
			tiempo: de {time_i} a {time_f} \n\
			latitud central: {lat_centro}º, \n\
			longitud central: {lon_centro}º \n\
			radio: {radio} km \n\
			profundidad: {prof_min}km a {prof_max} km \n\
			magnitud: {mag_min} a {mag_max},\n \
			maximo error en latitud: {elat_max} km \n\
			maximo error en longitud {elon_max} km \n\
			maximo error en profundidad {eprof_max} km\n\n")

	return r_propio #r contiene la matriz con los sismos filtrados



