# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:39:46 2018

@author: mlizarazo
#funcion que crea mapa base, pide cordenadas de recangulo en geograficas y como resultado retorna el mapa con coordenadas en los ejes y escala lateral izquierda
V3. Se estructuro por funciones, se incorporo mapa gap y rms, afecto graficar.py
"""



from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from operator import itemgetter
import os
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,AnnotationBbox)
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
import matplotlib.patches as patches
import matplotlib.dates as mdates

from obspy import UTCDateTime #importante para el condicional de tiempo en la fuincion selec_catalogo_seisan

import math
import MySQLdb
import matplotlib.colors
from matplotlib.colors import LogNorm
from copy import copy
import matplotlib.cm as cm
from matplotlib.dates import date2num

from scipy.ndimage.filters import gaussian_filter

from PIL import Image
from PIL import Image as Ima

from shapely.geometry import Point, Polygon
import utm

import requests
import urllib

import math
import geopy.distance

import pickle


def distancia_utm(lon1,lat1,lon2,lat2):
    
    #-------- halla distancia epicentral entre 2 eventos en utm
    coords_1 = (lat1,lon1)
    coords_2 = (lat2,lon2)
    dist_epi=round ( geopy.distance.geodesic(coords_1, coords_2).km , 3)
    return dist_epi


def graficar_fallas_rem_masa(mapa,ax,ancho_fallas,tam_tex_falla,en_negrilla):


    dir_fallas="/opt/rutinas/mapas_histogramas/catSeis_ests_leyens/FALLAS"
    myshape =(dir_fallas+'/Fallas')

    #sf = shp.Reader(myshape) #otra manera de leer shapes 
    #mapa.readshapefile(myshape,"faults", linewidth=ancho_fallas)
    mapa.readshapefile(myshape,"faults", linewidth=ancho_fallas)
        
    shapefile_vector = getattr(mapa, "faults") #trae las cordenadas
    shapefile_info = getattr(mapa, "faults_info") #trae la informacion del shape

    nombres=[]
    i=0
    #print("Fallas encontradas en shape de riesgo potencial:")
    for vector in shapefile_vector:
        id_falla=shapefile_info[i]["RINGNUM"]
        #print(shapefile_info[i])
        name_falla=shapefile_info[i]["STR_NAM"]
        #print(name_falla)
        x=[]
        y=[]

        for s in vector:
            x.append(s[0])
            y.append(s[1])
        
        eps_map_x=0.08*(ax.get_xlim()[1]-ax.get_xlim()[0])
        eps_map_y=0.08*(ax.get_ylim()[1]-ax.get_ylim()[0])
        
        lim_izq=ax.get_xlim()[0]+eps_map_x
        lim_der=ax.get_xlim()[1]-eps_map_x
        lim_inf=ax.get_ylim()[0]+(3*eps_map_y)
        lim_sup=ax.get_ylim()[1]-eps_map_y

        x_validos=[]
        y_validos=[]
        for j in range(0,len(x)-1):
            if lim_izq < x[j] < lim_der and lim_inf < y[j] < lim_sup:
                x_validos.append(x[j])
                y_validos.append(y[j])
        
        if len(x_validos)>=5:

            #ax.plot(x_validos,y_validos,"r",linewidth=5)

            if (x_validos[-1]-x_validos[0])!=0:
                angulo=np.arctan((y_validos[-1]-y_validos[0])/(x_validos[-1]-x_validos[0]))*(180/3.16)
            if (x_validos[-1]-x_validos[0])==0:
                angulo=90
            #print(angulo)
            pos=int(len(x_validos)/2)

            if lim_izq < x_validos[pos] < lim_der and lim_inf < y_validos[pos] < lim_sup:
                
                if name_falla not in nombres:

                    #print(name_falla)
                    if en_negrilla==True:
                        tx = ax.text(x_validos[pos], y_validos[pos], "F. "+name_falla, horizontalalignment='center', verticalalignment="center",fontsize=tam_tex_falla,rotation=angulo,fontweight="bold")
                    else:
                        tx = ax.text(x_validos[pos], y_validos[pos], "F. "+name_falla, horizontalalignment='center', verticalalignment="center",fontsize=tam_tex_falla,rotation=angulo)

                nombres.append(name_falla)
            
        i=i+1

def mapa_referencia(tipo_mapa,tam_map_ref,tipo_busqueda,radio,lon_centro,lat_centro,lat_max,lat_min,lon_max,lon_min,op,r):

    lon_min,lon_max,lat_min,lat_max,lat_min_ini=configura_geografia(tipo_busqueda,radio,lon_centro,lat_centro,lat_max,lat_min,lon_max,lon_min,op,r)
    
    print("Generando mapa_referencia...")
    
    tmg=tam_map_ref #numero de veces mas grande que el original
    del_lat=lat_max-lat_min
    del_lon=lon_max-lon_min

    lon_min_ref=round( lon_min-((tmg/2)*(del_lon)) , 3)
    lon_max_ref=round( lon_max+((tmg/2)*(del_lon)) , 3)
    lat_min_ref=round( lat_min-((tmg/2)*(del_lat)) , 3)
    lat_max_ref=round( lat_max+((tmg/2)*(del_lat)) , 3)

    #valores limite para hacer un mapa de referencia que conserve el epsg 3117
    if lon_min_ref < -105:
        lon_min_ref=-105

    if lon_max_ref > -28:
        lon_max_ref=-28

    if lat_min_ref < -34:
        lat_min_ref=-34

    if lat_max_ref > 48:
        lat_max_ref=48


    

    #---------------------
    fig, ax = plt.subplots(figsize=(15, 18)) #tamaño de la ventana del mapa    
    mapa_ref = Basemap(lon_min_ref,lat_min_ref,lon_max_ref,lat_max_ref, epsg=3117, resolution=None,projection='merc') #grafica el rectangulo de lat y lon max y min
        
    
    try:
    
        mapa_ref.arcgisimage(service="World_Shaded_Relief", xpixels = 1200, verbose= False, alpha=0.5) # # en http://server.arcgisonline.com/ArcGIS/rest/services hay varios provedores de imagenes
            
    except Exception as e:
        #continue
        print(e)
        print("Hubo problemas al generar el mapa de referencia---")
        mapa_ref.arcgisimage(service="NatGeo_World_Map", xpixels = 900, verbose= True, alpha=0.4)
    
    
    wms_server_sgc = "http://srvags.sgc.gov.co/arcgis/services/Mapa_Geologico_Colombia/Mapa_Geologico_Colombia/MapServer/WMSServer?request=GetCapabilities&service=WMS"
    
    try:    
        mapa_ref.wmsimage(wms_server_sgc, xpixels=3500, ypixels=3500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18"],verbose=False,alpha=1, transparent=True)
    
    except Exception as e:
        #continue
        print(e)
        print("Hubo problemas al generar el mapa de referencia---")
        mapa_ref.arcgisimage(service="NatGeo_World_Map", xpixels = 900, verbose= True, alpha=0.3)  
    
    

        
    x_c, y_c = mapa_ref( [lon_min, lon_max], [lat_min, lat_max] )
    
    def rec(mapa_ref,x_c,y_c,width,color):
        r1 = mapa_ref.plot([x_c[0],x_c[1]], [y_c[0],y_c[0]] , color,linewidth=width)
        r2 = mapa_ref.plot([x_c[0],x_c[1]], [y_c[1],y_c[1]] , color,linewidth=width)
        r3 = mapa_ref.plot([x_c[0],x_c[0]], [y_c[0],y_c[1]] , color,linewidth=width)
        r4 = mapa_ref.plot([x_c[1],x_c[1]], [y_c[0],y_c[1]] , color,linewidth=width)

    rec(mapa_ref,x_c,y_c,width=5,color="-k")

    plt.savefig('MAPA_REF.png', bbox_inches='tight')
    plt.close()

def agrega_mapa_ref(path_mapa,path_mapa_ref,zoom_map_ref,color_prof):

    pos_map_ref_x=0.82
    pos_map_ref_y=0.12

    if color_prof=="continuo":
        pos_map_ref_x=0.0
    

    path_mapa_sis=path_mapa
    path_mapa_ref=path_mapa_ref
    img1 = Ima.open(path_mapa_sis)

    img2 = Ima.open(path_mapa_ref)
    anch_original=float(img2.size[0])
    basewidth=int(zoom_map_ref*anch_original)

    if anch_original > basewidth:
        wpercent = (basewidth/float(img2.size[0]))
        hsize = int((float(img2.size[1])*float(wpercent)))
        img2 = img2.resize((basewidth,hsize), Ima.ANTIALIAS)

    x_offset = int(pos_map_ref_x*img1.size[0] + img2.size[0] ) #relacion valida para cuando el mapa ref excede el limite derecho del mapa de sismicidad
    if x_offset < img1.size[0]: #si el mapa ref esta muy a la izquierda la relacion de arriba no es correcta, porque el mapa ref esta dentro de los limites del mapa de sismicidad
        x_offset = img1.size[0]
    y_offset = int(img1.size[1]+ (pos_map_ref_y*img2.size[1]))

    new_im = Image.new('RGB', (x_offset, y_offset),color="white")    

    new_im.paste(img1, (0,int (pos_map_ref_y*img2.size[1]) ) )
    new_im.paste(img2, (int(pos_map_ref_x*img1.size[0]),0) ) 

    new_im.save(path_mapa)

def configura_geografia(tipo_busqueda,radio,lon_centro,lat_centro,lat_max,lat_min,lon_max,lon_min,op,r):

    if tipo_busqueda=="c":
        
        if radio < 10: #hace correcccion para generar un mapa mas grande que el pretendido, pues esri no puede traer mapas tan pequeños, 8 es el minimo encontrado
            radio = 10

        epsil=(2*radio/111) #grafica el rectangulo de lat y lon muy grande para alli buscar los datos, para hacer la tranformacion a planas, genera un mapa aproximadamente del doble del tamaño de la busqueda
        mapa_g = Basemap(lon_centro-epsil,lat_centro-epsil,lon_centro+epsil,lat_centro+epsil, epsg=3117, resolution=None) #grafica el rectangulo de lat y lon max y min
        x_centro,y_centro=mapa_g(lon_centro,lat_centro)

        lon_min_cuad, lat_min_cuad = mapa_g(x_centro - (radio*1000),y_centro - (radio*1000),inverse=True)
        lon_max_cuad, lat_max_cuad = mapa_g(x_centro + (radio*1000),y_centro + (radio*1000),inverse=True)
        
        epsilon=0.15 #15% aumenta hacia arriba, abajo, derecha e izquierda despecto al tamaño del circulo, si se quiere un mapa mas grande aumentar
        epsilon_lat=(lat_max_cuad-lat_min_cuad)*epsilon
        lat_min=lat_min_cuad-epsilon_lat
        lat_max=lat_max_cuad+epsilon_lat

        epsilon_lon=(lon_max_cuad-lon_min_cuad)*epsilon
        lon_min=lon_min_cuad-epsilon_lon
        lon_max=lon_max_cuad+epsilon_lon

    if op=="ap" and tipo_busqueda=="sin_filtro_loc":
        #sin filtro loc, toma las lat, lon y prof min y max como las del propio catalogo con el fin de que se analice todo el contenido
        #del catalogo sin imponer restricciones espaciales, deja restrictivo el resto de variables menos las de posicion, funciona unicamente para base de datos propia
        
        lats=[]
        lons=[]
        profs=[]
        for i in r:

            lats.append(float(i[2]))
            lons.append(float(i[3]))
            profs.append(float(i[4]))
        
        lat_min=min(lats)
        lat_max=max(lats)
        lon_min=min(lons)
        lon_max=max(lons)

        del_lat=lat_max-lat_min
        del_lon=lon_max-lon_min

        if del_lon < 0.28 or del_lat < 0.28: #corrige dimensiones muy pequeñas del mapa, para que tenga como minimo 0.2 grados
            agrg_lat=(0.28 - del_lat)/2
            agrg_lon=(0.28 - del_lon)/2

            lat_min=lat_min - agrg_lat
            lat_max=lat_max + agrg_lat
            lon_min=lon_min - agrg_lon
            lon_max=lon_max + agrg_lon


        prof_min=min(profs)
        prof_max=max(profs)

    #para el caso rectangular ha buscado los eventos hasta unos limites minimos y maximos, luego para graficar necesita exceder esos limites para que queden todos los sismos dentro del mapa, agregara por defecto el 5% de las mediad
    #en cada extremo
    if tipo_busqueda=="r" or tipo_busqueda=="p" or tipo_busqueda=="sin_filtro_loc": 
        eps_lat_rec=(lat_max-lat_min)*0.05
        eps_lon_rec=(lon_max-lon_min)*0.05

        lat_min=lat_min-eps_lat_rec
        lat_max=lat_max+eps_lat_rec
        lon_min=lon_min-eps_lon_rec
        lon_max=lon_max+eps_lon_rec

        del_lat=lat_max-lat_min
        del_lon=lon_max-lon_min

        if del_lon < 0.28 or del_lat < 0.28: #corrige dimensiones muy pequeñas del mapa, para que tenga como minimo 0.2 grados
            agrg_lat=(0.28 - del_lat)/2
            agrg_lon=(0.28 - del_lon)/2

            lat_min=lat_min - agrg_lat
            lat_max=lat_max + agrg_lat
            lon_min=lon_min - agrg_lon
            lon_max=lon_max + agrg_lon
    #---------------------    
    #re calcula maximos y minimos de las coordenadas lat y lon cuando estan muy desproporcionadas generando mapas muy 
    #rectangulares, se toma como desproporcionado una relaciones mayores a 1.7, se arreglan los nuevos maximos y minimos generando relaciones de 1.2
    if lat_max-lat_min > 1.7*(lon_max-lon_min):
        
        del_lon_new=(lat_max-lat_min)/1.2
        eps_lon_new=(del_lon_new-(lon_max-lon_min))/2
        lon_min=lon_min-eps_lon_new
        lon_max=lon_max+eps_lon_new

    if 1.7*(lat_max-lat_min) < (lon_max-lon_min):
        del_lat_new=(lon_max-lon_min)/1.2
        eps_lat_new=(del_lat_new-(lat_max-lat_min))/2
        lat_min=lat_min-eps_lat_new
        lat_max=lat_max+eps_lat_new
    #----------------

    #se agrega una porcion (1/8 parte) de mapa hacia abajo, para que alli grafique la leyenda
    lat_min_ini=lat_min      
    lat_min=lat_min-((lat_max-lat_min)/8) 

    return lon_min,lon_max,lat_min,lat_max,lat_min_ini

def asigna_color_pepa(r_sort): #colors y existe_prof son listas que se entregan vacias, para que esta funcion las llene
    
    colors=[]
    existe_prof=[] #sera util para graficar leyenda de profundidad
        
    for i in r_sort:

        prof=float(i[4])                
                
        if prof <= 30:
            color="#ff0000" #rojo
            colors.append(color)
            existe_prof.append(30)
        if 30 < prof <= 70:
            color="#ffff00" #amarillo
            colors.append(color)
            existe_prof.append(70)
        if 70 < prof <= 120:
            color="#32cd32" #verde
            colors.append(color)
            existe_prof.append(120)
        if 120 < prof <= 180:
            color="#0000ff" #azul
            colors.append(color)
            existe_prof.append(180)
        if 180 < prof <= 1000:
            color="#8000ff" #morado
            colors.append(color)
            existe_prof.append(1000)

    return colors, existe_prof

def colores_y_texto_ley_prof(vec_colores_unic,long_barra):
    colores_ley_prof=[]
    txt_ley_prof=[]            
    
    for e in vec_colores_unic:
        if e==30:
            colores_ley_prof.append("#ff0000")
            txt_ley_prof.append("0-30")
        if e==70:
            colores_ley_prof.append("#ffff00")
            txt_ley_prof.append("30-70")
        if e==120:
            colores_ley_prof.append("#32cd32")
            txt_ley_prof.append("70-120")
        if e==180:
            colores_ley_prof.append("#0000ff")
            txt_ley_prof.append("120-180")
        if e==1000:
            colores_ley_prof.append("#8000ff")
            txt_ley_prof.append(">180")

    if len(vec_colores_unic)==1:

        del_prof=0.001*long_barra ##
        s_scatt=570
        text_font=26

    if len(vec_colores_unic)==2:

        del_prof=0.00045*long_barra ##
        s_scatt=560
        text_font=26

    if len(vec_colores_unic)==3:

        del_prof=0.00037*long_barra ##
        s_scatt=550
        text_font=26

    if len(vec_colores_unic)==4:

        del_prof=0.00033*long_barra ##
        s_scatt=int(590-(0.45*long_barra))
        text_font=25
        #text_font=int(29-((3/100)*long_barra))

    if len(vec_colores_unic)==5:

        del_prof=0.00031*long_barra ##
        s_scatt=int(585-(0.5*long_barra))
        text_font=24
        #text_font=int(29-((3/100)*long_barra))
    
    return colores_ley_prof,txt_ley_prof,del_prof,s_scatt,text_font

def mapa_profundidad(r_sort,long_barra,lon_min,lon_max,loc_x_txt,pos_mag_x,loc_y_txt,mapa,x,y,lat_pos_prof,lat_min,alpha,profs,mags_map,ax,color_prof,tam_text,y_mag,cmap):
    
    print ("Generando mapa de sismicidad...")
    pos_mag_y=y_mag+(0.00104*long_barra)  ##

    if color_prof=="discreto":

        colors, existe_prof=asigna_color_pepa(r_sort)
        aleja_prof_mag=1.7 #entre mas grande sea corre la leyenda hacia la izquierda, alejandose proporcionalmente de la leyenda de magnitud
        
        lon_pos_lp=lon_min+((pos_mag_x/aleja_prof_mag)*(lon_max-lon_min))
        evs = mapa.scatter(x, y, color=colors, edgecolor='k', s=mags_map, zorder=10,alpha=alpha)
        
        vec_colores_unic=sorted(set(existe_prof))            
        #-------colores y texto leyenda prof
        colores_ley_prof,txt_ley_prof,del_prof,s_scatt,text_font=colores_y_texto_ley_prof(vec_colores_unic,long_barra)

        for i in range(0,len(colores_ley_prof)):
            lat_pos_lp=lat_pos_prof+lat_min+(1.6*pos_mag_y)-(del_prof*i)
            x_lp,y_lp=mapa(lon_pos_lp,lat_pos_lp)
            mapa.scatter(x_lp, y_lp,color=colores_ley_prof[i], s=s_scatt, edgecolors='k',zorder=20, linewidth=1,alpha=0.8)
            plt.text(x_lp+((30*long_barra)-100), y_lp,txt_ley_prof[i], fontsize=text_font,horizontalalignment='left',verticalalignment='center') #5 tamaño texto de MAGNITUD
        #ley
        ax.text(0.68*(pos_mag_x/aleja_prof_mag),loc_y_txt, "PROFUNDIDAD (Km)", transform=ax.transAxes, fontsize=tam_text, horizontalalignment='left')
        ax.text(loc_x_txt+pos_mag_x,loc_y_txt, "MAGNITUD", transform=ax.transAxes, fontsize=tam_text, horizontalalignment='left')
    
    if color_prof=="continuo":        
        
        evs = mapa.scatter(x, y, c=profs, edgecolor='k', s=mags_map, zorder=10,alpha=alpha,cmap=cmap)            
        cbar = plt.colorbar(orientation='vertical', shrink=0.6, aspect=30, fraction=0.2,pad=0.01)
        cbar.set_label("\nProfundidad (km)",size=40)
        cbar.ax.tick_params(labelsize=30)
        #ley
        ax.text(loc_x_txt+pos_mag_x,loc_y_txt, "MAGNITUD", transform=ax.transAxes, fontsize=tam_text, horizontalalignment='left')        

    return mapa, ax

def mapa_gap(y_mag,long_barra,mapa,ax,x,y,gaps,mags_map,alpha,loc_x_txt,pos_mag_x,loc_y_txt,tam_text,cmap):
    
    print ("Generando mapa de GAP...")
                
    
    evs = mapa.scatter(x, y, c=gaps, edgecolor='k', s=mags_map, zorder=10,alpha=alpha,cmap=cmap)            
    cbar = plt.colorbar(orientation='vertical', shrink=0.6, aspect=30, fraction=0.2,pad=0.01)
    cbar.set_label("\nGap (°)",size=40)
    cbar.ax.tick_params(labelsize=30)
    #ley
    ax.text(loc_x_txt+pos_mag_x,loc_y_txt, "MAGNITUD", transform=ax.transAxes, fontsize=tam_text, horizontalalignment='left')        

    return mapa, ax,

def mapa_rms(y_mag,long_barra,mapa,ax,x,y,rmss,mags_map,alpha,loc_x_txt,pos_mag_x,loc_y_txt,tam_text,cmap):
    
    print ("Generando mapa de RMS...")
                
    
    evs = mapa.scatter(x, y, c=rmss, edgecolor='k', s=mags_map, zorder=10,alpha=alpha,cmap=cmap)            
    cbar = plt.colorbar(orientation='vertical', shrink=0.6, aspect=30, fraction=0.2,pad=0.01)
    cbar.set_label("\nRMS (s)",size=40)
    cbar.ax.tick_params(labelsize=30)
    #ley
    ax.text(loc_x_txt+pos_mag_x,loc_y_txt, "MAGNITUD", transform=ax.transAxes, fontsize=tam_text, horizontalalignment='left')        

    return mapa, ax

def mapa_time(y_mag,long_barra,mapa,ax,x,y,time,mags_map,alpha,loc_x_txt,pos_mag_x,loc_y_txt,tam_text,cmap):
    
    print ("Generando mapa de tiempo...")

    del_t=(max(time)- min(time)).total_seconds()/60/60/24 #en dias              
    
    datenums=[]    
    for t in time:
        datenums.append(date2num(t))        

    dic={"tiempo": time}
    df_tiempo = pd.DataFrame(dic)     

    evs = mapa.scatter(x, y, c=datenums, edgecolor='k', s=mags_map, zorder=10,alpha=alpha,cmap=cmap)

    
    ticks=[min(datenums), max(datenums)]
    cbar=plt.colorbar(evs,ticks=ticks, shrink=0.6, aspect=30, fraction=0.2,pad=0.01)
        
    if del_t < 5:
        escala_t="%Y-%m-%d %H:%M:%S"
    if 5 <= del_t < 1000:
        escala_t="%Y-%m-%d"
    if 1000 <=del_t < 3000:
        escala_t="%Y-%m"
    if del_t >= 3000:
        escala_t="%Y"

    fech_tiks=[min(time).strftime(escala_t),max(time).strftime(escala_t)]
    cbar.ax.set_yticklabels(fech_tiks,fontsize=30)  
    cbar.set_label("\nFecha",size=40)   
    

    #ley
    ax.text(loc_x_txt+pos_mag_x,loc_y_txt, "MAGNITUD", transform=ax.transAxes, fontsize=tam_text, horizontalalignment='left')        

    return mapa, ax

def leyenda_magnitudes(lon_min,lon_max,lat_min,pos_mag_x,d_mag,long_barra,mags,FAM,exp_mag,mapa,dir_ley_mag,y_mag):
    pos_mag_y=y_mag+(0.00104*long_barra)  ##
    lon_pos=lon_min+(pos_mag_x*(lon_max-lon_min))
    del_mag=d_mag*0.000084*long_barra ##

    mag_ley=[]
    lat_ley=[]
    lon_ley=[]
    a=0
    b=0
    
    max_mag=math.ceil(max(mags)) #busca el maximo valor de magnitud y lo aproxima al entero proximo mayor

    if min(mags) < 1: #decide cual es la menor magnitud para ponerla como leyenda, para no graficar el cero si la menor magnitud es menor a 1, en el mapa la menor magnitud sera un 1
        min_mag=1
    else:
        min_mag=int(min(mags)) 
    
    if dir_ley_mag=="h":
    #leyenda H        
        for i in range(min_mag,max_mag+1):
            
            lat_ley.append(lat_min+pos_mag_y) 
            lon_ley.append(lon_pos+(i*a))
            
            #-------------------------------
            mag_ley.append((((i/8)*(8-FAM))+FAM)**exp_mag) 
            #-------------------------------        

            x_tx,y_tx=mapa(lon_ley[b],lat_ley[b])
            if i<4:
                plt.text(x_tx, y_tx-((80*long_barra)-100),str(i),fontsize=5*(i+3),horizontalalignment='center',verticalalignment='center') #80*long_barra)-100 es una relacion lineal que baja los numeros 1,2 y 3 en la escala de magnitud dependiendo del tamaño del mapa
            else:
                plt.text(x_tx, y_tx,str(i),fontsize=5*(i+2),horizontalalignment='center',verticalalignment='center') #5 tamaño texto de MAGNITUD
            mapa.scatter(x_tx, y_tx,facecolors='none', s=mag_ley[b], edgecolors='k',zorder=20, linewidth=2)

            a=a+del_mag 
            b=b+1
    
    if dir_ley_mag=="v":
        #leyenda V
        for i in range(min_mag,max_mag+1):
            lat_ley.append(lat_min+((0.005+a)*(lat_max-lat_min)))
            lon_ley.append(lon_pos)
            
            #-------------------------------
            mag_ley.append((((i/8)*(8-FAM))+FAM)**exp_mag) 
            #---------------------------------
            
            x_tx,y_tx=mapa(lon_ley[b],lat_ley[b])
            plt.text(x_tx, y_tx,str(i),fontsize=4*i,horizontalalignment='center',verticalalignment='right')
            mapa.scatter(x_tx, y_tx,facecolors='none', s=mag_ley[b], edgecolors='k',zorder=1)

            a=a+0.01
            b=b+1
            
    x_ley, y_ley = mapa(lon_ley, lat_ley)
    puntos_mag = mapa.scatter(x_ley, y_ley, c='b', alpha=0.2, s=mag_ley, zorder=20, linewidth=2)

    return mapa

def proyeccion_sismo_a_perfil(lon_1,lat_1,lon_2,lat_2,lat,lon):
    dis_sis=distancia_utm(lon_1,lat_1,lon,lat)
    longitud_perfil=distancia_utm(lon_1,lat_1,lon_2,lat_2)
    longitud_perfil_x=distancia_utm(lon_1,lat_1,lon_2,lat_1)
    longitud_perfil_y=distancia_utm(lon_1,lat_1,lon_1,lat_2)
    dis_sis_x=distancia_utm(lon_1,lat_1,lon,lat_1)
    dis_sis_y=distancia_utm(lon_1,lat_1,lon_1,lat)
    
    if lon_1 == lon_2:        
        theta=3.1416/2
    else:
        theta=math.atan(longitud_perfil_y/longitud_perfil_x)
    
    beta=math.atan(dis_sis_y/dis_sis_x)

    if lat_1 != lat_2 and lon_1 != lon_2:
    
        if lat_1 > lat_2: #perfil con pendiente negativa
            
            if lat_1<=lat:
                angulo=theta+beta

            if lat_1>lat and theta>=beta:
                angulo=theta-beta

            if lat_1>lat and theta<beta:
                angulo=-(theta-beta)

        if lat_1 < lat_2: #perfil con pendiente positiva
            
            if lat_1>=lat:
                angulo=theta+beta

            if lat_1<lat and theta>=beta:
                angulo=theta-beta

            if lat_1<lat and theta<beta:
                angulo=-(theta-beta)

        proyec=round( dis_sis*(math.cos(angulo)) , 2)

    if lat_1 == lat_2: #perfil horizontal
        proyec=dis_sis_x

    if lon_1 == lon_2: #perfil vetical
        proyec=dis_sis_y


    
    #print(theta,beta,dis_sis,proyec)       

    return proyec

def construye_gepgrafia_perfil(ancho,lon_1,lat_1,lon_2,lat_2):
    
    #construye cordenadas del perfil

    longitud_perfil_x=distancia_utm(lon_1,lat_1,lon_2,lat_1)
    longitud_perfil_y=distancia_utm(lon_1,lat_1,lon_1,lat_2)

    if lon_1 == lon_2:        
        theta=3.1416/2
    else:
        theta=math.atan(longitud_perfil_y/longitud_perfil_x)

    del_lon_e=round( (ancho/2)*math.cos( (3.1416/2) - theta) , 3)
    del_lat_e=round( (ancho/2)*math.sin( (3.1416/2) - theta)  ,3)

    #print(longitud_perfil_x,longitud_perfil_y,theta,del_lat_e)

    if lat_1 > lat_2: #perfil pendiente negativa 
        lon_pol_1=lon_1+del_lon_e
        lat_pol_1=lat_1+del_lat_e

        lon_pol_2=lon_2+del_lon_e
        lat_pol_2=lat_2+del_lat_e

        lon_pol_3=lon_2-del_lon_e
        lat_pol_3=lat_2-del_lat_e

        lon_pol_4=lon_1-del_lon_e
        lat_pol_4=lat_1-del_lat_e

    if lat_1 < lat_2: #perfil pendiente positivo
        lon_pol_1=lon_1-del_lon_e
        lat_pol_1=lat_1+del_lat_e

        lon_pol_2=lon_2-del_lon_e
        lat_pol_2=lat_2+del_lat_e

        lon_pol_3=lon_2+del_lon_e
        lat_pol_3=lat_2-del_lat_e

        lon_pol_4=lon_1+del_lon_e
        lat_pol_4=lat_1-del_lat_e

    if lat_1 == lat_2: #perfil pendiente horizontal
        lon_pol_1=lon_1
        lat_pol_1=lat_1+del_lat_e

        lon_pol_2=lon_2
        lat_pol_2=lat_2+del_lat_e

        lon_pol_3=lon_2
        lat_pol_3=lat_2-del_lat_e

        lon_pol_4=lon_1
        lat_pol_4=lat_1-del_lat_e

    if lon_1 == lon_2: #perfil pendiente vertical
        lon_pol_1=lon_1-del_lon_e
        lat_pol_1=lat_1

        lon_pol_2=lon_1+del_lon_e
        lat_pol_2=lat_1

        lon_pol_3=lon_1+del_lon_e
        lat_pol_3=lat_2

        lon_pol_4=lon_1-del_lon_e
        lat_pol_4=lat_2

    lon_pol_1=round (lon_pol_1, 3)
    lat_pol_1=round (lat_pol_1, 3)
    lon_pol_2=round (lon_pol_2, 3)
    lat_pol_2=round (lat_pol_2, 3)
    lon_pol_3=round (lon_pol_3, 3)
    lat_pol_3=round (lat_pol_3, 3)
    lon_pol_4=round (lon_pol_4, 3)
    lat_pol_4=round (lat_pol_4, 3)


    p1_pol=(lon_pol_1,lat_pol_1)
    p2_pol=(lon_pol_2,lat_pol_2)
    p3_pol=(lon_pol_3,lat_pol_3)
    p4_pol=(lon_pol_4,lat_pol_4)

    #print(p1_pol,p2_pol,p3_pol,p4_pol)
    return p1_pol,p2_pol,p3_pol,p4_pol


def grafica_perfil(lon_1,lat_1,lon_2,lat_2,r,FAM,exp_mag,alpha,ancho,color_prof,cmap,lp):
#def grafica_perfil(lon_1,lat_1,lon_2,lat_2,r,FAM,exp_mag,*args,**kwargs)
    
        

    p1_pol,p2_pol,p3_pol,p4_pol=construye_gepgrafia_perfil(ancho,lon_1,lat_1,lon_2,lat_2)
    perfil = Polygon([p1_pol,p2_pol,p3_pol,p4_pol])
    cords_perfil=perfil.exterior.coords[:]
    #plt.plot(*perfil.exterior.xy,"r",linewidth=0.5) #manera de pintar poligonos
    r_sort=sorted(r, key=itemgetter(5),reverse=True) #(4) pone los mas profundos abajo, (5) pone los mas grandes abajo. True o False para que clasifique de mayor a menor o visceverza

    lats=[]
    lons=[]
    mags_map=[]
    mags=[]    
    profs=[]
    profs_pos=[]
    elats=[]
    elons=[]
    eprofs=[]
    diss=[]
    profs_total_cat=[]

    r_perfil=[]
    for i in r_sort:
        lat=float(i[2])
        lon=float(i[3])
        profs_total_cat.append(-float(i[4]))
        p = Point(lon, lat)
        if p.within(perfil)==True:
            prof=float(i[4])
            profs.append(-prof)
            profs_pos.append(prof) #servira para el perfil 
            #dis=distancia_utm(lon_1,lat_1,lon,lat)
            proyec=proyeccion_sismo_a_perfil(lon_1,lat_1,lon_2,lat_2,lat,lon)
            #print(dis,dis2)
            diss.append(proyec)
            mag=float(i[5])        
            mags_map.append((((mag/8)*(8-FAM))+FAM)**exp_mag) 
            mags.append(mag)
            r_perfil.append(i)

   

    if len(profs) >0:
        longitud_perfil=distancia_utm(lon_1,lat_1,lon_2,lat_2)
        prof_perfil=abs(min(profs)) + 1 #1 km sobre el suelo, ojo min(prof) es la profunddad maxima y esta negativa, por eso el abs
        #prof_perfil=100 

        #----tamaño de grafica
        if longitud_perfil >= prof_perfil:

            tx=12
            ty=(6.2/4.8)*tx*prof_perfil/longitud_perfil
            #ty=int(25/(longitud_perfil/prof_perfil))

        if longitud_perfil < prof_perfil:
            ty=10
            tx=(4.8/6.2)*ty*longitud_perfil/prof_perfil
            #tx=int(20/(prof_perfil/longitud_perfil))
        
        
        fig = plt.figure(figsize=(tx, ty))
        #fig = plt.figure.gcf()

        #----------------
        
        ax = fig.add_subplot(111)
        ax.axis('equal')
        ax.xaxis.grid(True, which='major')
        ax.yaxis.grid(True, which='major')
        ax.set_xlim( 0, longitud_perfil)
        ax.hlines(0,ax.get_xlim()[0],ax.get_xlim()[1],"k",linewidth=1,linestyles="dashed")
        ax.set_ylim( min(profs), 2)
        #fig.set_size_inches(tx, ty)
        

        if color_prof=="discreto":

            colors, existe_prof=asigna_color_pepa(r_perfil)
            ax.scatter(diss,profs,color=colors,edgecolor='k',s=mags_map,zorder=10,alpha=alpha)
            
        if color_prof=="continuo":
            #adiciona los valores minimo y maxmimo de la pofundiad para los colores creados incluyan dichos valores, asi concidiran los colores en el mapa y en el perfil
            diss=diss+[0,0]
            profs=profs + [min(profs_total_cat),max(profs_total_cat)]
            profs_pos= profs_pos + [-max(profs_total_cat),-min(profs_total_cat)]
            mags_map=mags_map+[0,0]
            ax.scatter(diss,profs,c=profs_pos,edgecolor='k',s=mags_map,zorder=10,alpha=alpha,cmap=cmap)                  
                            
        
        #ax.scatter(diss,profs,edgecolor='k',s=mags_map,zorder=10,alpha=alpha)
        ymin, ymax = plt.ylim()
        xmin, xmax = plt.xlim()
        ax.text(xmin,0,lp, fontsize=25,horizontalalignment='center', color="k", style='italic')
        ax.text(xmax,0,lp+"'", fontsize=25,horizontalalignment='center', color="k", style='italic')
       

        ax.set_xlabel('\nDistancia (km)', color='k', fontsize=25)
        ax.set_ylabel('\nProfundidad (km)', color='k', fontsize=25)
        ax.xaxis.set_tick_params(labelsize=20)
        ax.yaxis.set_tick_params(labelsize=20)
        #plt.title("Perfil de sismicidad", color='k', fontsize=30)
        plt.savefig('Perfil_'+lp+'.png', bbox_inches='tight')
        plt.close()

        

    return cords_perfil,r_perfil

def genera_perfil(mapa,ax,lon_1,lat_1,lon_2,lat_2,r,FAM,exp_mag,alpha,ancho,color_prof,cmap,lp):

    cords_perfil, r_perfil = grafica_perfil(lon_1,lat_1,lon_2,lat_2,r,FAM,exp_mag,alpha,ancho,color_prof,cmap,lp)
    

    lats_perf=[]
    lons_perf=[]
    for cp in cords_perfil:
        lons_perf.append(cp[0])
        lats_perf.append(cp[1])
    
    x_perf,y_perf=mapa(lons_perf, lats_perf) 
    
    mapa.plot(x_perf, y_perf, '--k',linewidth=3,zorder=100)

    x_text,y_text=mapa([lon_1,lon_2], [lat_1,lat_2])
    ax.text(x_text[0],y_text[0], lp, fontsize=40, horizontalalignment='center', color="k", style='italic',zorder=100)
    ax.text(x_text[1],y_text[1], lp+"'", fontsize=40, horizontalalignment='center', color="k", style='italic',zorder=100)
    
    
    return mapa,ax

def adicion_contorno_poligono(mapa,ax,archivo_poligono):
    
    arc_pol=open(archivo_poligono,"r")
    lineas_pol=arc_pol.readlines()
    lons_pol=[]
    lats_pol=[]
    for lc in lineas_pol:
        vec_lc=lc.split()
        lons_pol.append(float(vec_lc[0]))
        lats_pol.append(float(vec_lc[1]))
    x_c,y_c=mapa(lons_pol, lats_pol)
    pol = mapa.plot(x_c, y_c, 'k', linewidth=3)
    '''
    for i in range(0,len(x_c)):
        ax.text(x_c[i],y_c[i], str(i), fontsize=40, horizontalalignment='center', color="k", style='italic')
    cir = mapa.plot(x_c, y_c, 'k', linewidth=3)     #---
    '''
    return mapa,ax


def mapa_crudo(r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,color_prof):

    
    print("Generando mapa crudo...")
    
    #-------
    lon_min,lon_max,lat_min,lat_max,lat_min_ini=configura_geografia(tipo_busqueda,radio,lon_centro,lat_centro,lat_max,lat_min,lon_max,lon_min,op,r)
    
    #---------------------
    fig, ax = plt.subplots(figsize=(33, 45)) #tamaño de la ventana del mapa    
    mapa = Basemap(lon_min,lat_min,lon_max,lat_max, epsg=3117, resolution=None,projection='merc') #grafica el rectangulo de lat y lon max y min
    
        
    #mapa_frp=mapa
    
    x_map,y_map=mapa([lon_min,lon_max],[lat_min,lat_max])

    try:
    
        if (lon_min < -82.0 or lat_max > 14.0) and (tipo_mapa=="Elevation/World_Hillshade" or tipo_mapa== "World_Shaded_Relief" or tipo_mapa=="Canvas/World_Light_Gray_Base" ):
            print("Obteniendo mapa base continental y oceanico...")
            mapa.arcgisimage(service=tipo_mapa, xpixels = 1800, verbose= False, alpha=1) # # en http://server.arcgisonline.com/ArcGIS/rest/services hay varios provedores de imagenes
            
            try:
                #print("h")    
                mapa.arcgisimage(service="Ocean/World_Ocean_Base", verbose= False, alpha=0.4) #notese que no tiene tamaño (xpixels) para que el lo construya por defecto, si se le pone 900 supera la resolucion del mapa y produce error
            except Exception as e:
                #continue
                print(e)
                print("Hubo problemas al generar el mapa base oceanico")
            
        else:
            
            mapa.arcgisimage(service=tipo_mapa, xpixels = 1800, verbose= False, alpha=1) # en http://server.arcgisonline.com/ArcGIS/rest/services hay varios provedores de imagenes
            print("Obteniendo mapa base ...")
            
        #ESRI_Imagery_World_2D (MapServer)     #ESRI_StreetMap_World_2D (MapServer)     #NatGeo_World_Map (MapServer)     #NGS_Topo_US_2D (MapServer)
        #Ocean_Basemap (MapServer)     #USA_Topo_Maps (MapServer)     #World_Imagery (MapServer)     #World_Physical_Map (MapServer)     #World_Shaded_Relief (MapServer)
        #World_Street_Map (MapServer)     #World_Terrain_Base (MapServer)     #World_Topo_Map (MapServer     #Canvas/World_Light_Gray_Base      #Elevation/World_Hillshade

    except Exception as e:
        #continue
        print(e)
        print("Hubo problemas al generar el mapa base, se intentara generar NatGeo_World_Map por defecto")

        try:
            mapa.arcgisimage(service="NatGeo_World_Map", xpixels = 900, verbose= True, alpha=0.4)
        except Exception as e:        
            print(e)
            print("Hubo problemas al generar el mapa base NatGeo_World_Map, se intentara graficar en fondo gris")
            mapa.drawlsmask(alpha=0.2)
    
    try:
        wms_server_sgc = "http://srvags.sgc.gov.co/arcgis/services/Mapa_Geologico_Colombia/Mapa_Geologico_Colombia/MapServer/WMSServer?request=GetCapabilities&service=WMS"
        
        '''
        1. fronteras     2. drenaje doble     3. lagunas     4. Bogota: grafica el poligono de bogota     5. Capitales Principales     6. drenaje sencillo
        7. Embalse     8. Departamentos     9. Ciénaga     10. Pueblos #pone los puntos de los pueblos, mas no sus anotaciones 11. Anotaciones Mapa Base
        12. A_Hidrografía     13. A_Limite maritimo      14. A_Límite Colombia    15. Pliegues    16. Fallas    17. Anotaciones Geología    18. A_Template
        19. A_Edificios volcánicos    20. A_Unidades cronoestratigráficas    21. A_Pliegues    22. A_Fallas    23. A_Índice volcanes    24. A_Elementos unidades cronoestratigráficas
        
        '''
        
        
        if tipo_mapa=="Elevation/World_Hillshade" or tipo_mapa=="World_Shaded_Relief" or tipo_mapa=="World_Imagery" or tipo_mapa=="NGS_Topo_US_2D" or tipo_mapa== "World_Physical_Map" or tipo_mapa=="World_Terrain_Base": #estos tipos de mapas no tienen departamentos ni municipios ni drenajes

            if incluir_fallas_sgc==True and incluir_fallas_movmasa==False:
                print("Obteniendo elementos geograficos y fallas del WMS del SGC ...")            
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18","16","22"],verbose=False,alpha=1, transparent=True, width=3000, height=3000)
                
            if incluir_fallas_sgc==False and incluir_fallas_movmasa==True:
                print("Obteniendo elementos geograficos del WMS del SGC y fallas del shape de fallas con riesgo potencial...")
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18"],verbose=False,alpha=1, transparent=True)
                graficar_fallas_rem_masa(mapa,ax,ancho_fallas,tam_tex_falla,en_negrilla)
                
                
            if incluir_fallas_sgc==True and incluir_fallas_movmasa==True:
                print("Obteniendo elementos geograficos y fallas del WMS del SGC y fallas del shape de fallas con riesgo potencial...")
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18","16","22"],verbose=False,alpha=1, transparent=True)
                graficar_fallas_rem_masa(mapa,ax,ancho_fallas,tam_tex_falla,en_negrilla)

            if incluir_fallas_sgc==False and incluir_fallas_movmasa==False:
                print("Obteniendo elementos geograficos del WMS del SGC...")
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18"],verbose=False,alpha=1, transparent=True)
                                
        else:
            if incluir_fallas_sgc==True and incluir_fallas_movmasa==False:
                print("Obteniendo elementos geograficos y fallas del WMS del SGC ...")            
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["16","22"],verbose=False,alpha=1, transparent=True )
                
            if incluir_fallas_sgc==False and incluir_fallas_movmasa==True:
                print("Obteniendo fallas del shape de fallas con riesgo potencial...")
                if tipo_mapa != "NatGeo_World_Map":
                    mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18"],verbose=False,alpha=1, transparent=True)
                graficar_fallas_rem_masa(mapa,ax,ancho_fallas,tam_tex_falla,en_negrilla)
                
            if incluir_fallas_sgc==True and incluir_fallas_movmasa==True:
                print("Obteniendo elementos geograficos y fallas del WMS del SGC y fallas del shape de fallas con riesgo potencial...")
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18","16","22"],verbose=False,alpha=1, transparent=True)
                graficar_fallas_rem_masa(mapa,ax,ancho_fallas,tam_tex_falla,en_negrilla)
            
            if incluir_fallas_sgc==False and incluir_fallas_movmasa==False:
                print("Obteniendo elementos geograficos del WMS del SGC...")
                mapa.wmsimage(wms_server_sgc, xpixels=5500, ypixels=5500,layers=["1","2","3","4","5","6","7","8","9","11","12","13","14","18"],verbose=False,alpha=1, transparent=True)    
        
    except Exception as e:
        #continue
        print(e)
        print("Hubo problemas al generar el mapa base")
        
    
    
    #FORMATO DE GRILLA DE LAT Y LON escribe los valores de lalitud y longitud cada delta_lat y delta_lon respectivamente
    tam_let_grill=28
    delta_lat=round((lat_max-lat_min)/5,1) #divide en 5 el intervalo de latitud y lo redondea a 1 decimal
    if delta_lat==0:
        delta_lat=0.1
    parallels = np.arange(lat_min,lat_max,delta_lat)
    mapa.drawparallels(parallels,labels=[True,False,False,False], dashes=[6,6], linewidth=0.1,fontsize=tam_let_grill, fmt='%0.4g')

    delta_lon=round((lon_max-lon_min)/5,1)
    if delta_lon==0:
        delta_lon=0.1
    meridians = np.arange(lon_min,lon_max,delta_lon)
    mapa.drawmeridians(meridians,labels=[False,True,True,False], dashes=[6,6], linewidth=0.1,fontsize=tam_let_grill, fmt='%0.5g')

    if tipo_busqueda=="c": 
        
        x_centro,y_centro=mapa(lon_centro,lat_centro) 
        
        if indica_centro_busqueda==True:

            evs = mapa.scatter(x_centro, y_centro, marker='*',color='k', edgecolor='w', s=1050, zorder=10) #indica el centro de la buqueda
            Circulo_busqueda = plt.Circle((x_centro,y_centro), 32*1000, color='blue',fill=False, linestyle="--",alpha=0.8)
            #ax.add_patch(Circulo_busqueda) #quitar esta linea si no se quiere resaltar el area de busqueda
    
    if incluir_rectangulo_leyenda==True:

        delYgrande=max(y_map)-min(y_map)
        x_rec, y_rec = mapa(lon_min, lat_min_ini)
        delYpequeno=y_rec-min(y_map)
        ancho_rec=(delYpequeno/delYgrande)-0.015 #si el rectangulo sale tapado por sismos, restarle un poco mas al ancho del rectangulo
        rectangulo = patches.Rectangle((0, 0), 1, ancho_rec,fill=True, transform=ax.transAxes, clip_on=False, facecolor='w',edgecolor='k',linewidth=2) #
        ax.add_patch(rectangulo)

    
    #trae de directorio el logo del servicio y leyendas de tamaños
    if color_prof=="continuo":
        zoom_leyenda=zoom_leyenda-0.15


    archivo_logo="/opt/rutinas/mapas_histogramas/catSeis_ests_leyens/leyenda2.png"
    logo = plt.imread(archivo_logo, format='png')
    imagebox = OffsetImage(logo, zoom=zoom_leyenda,transform=ax.transAxes) #modificar zoom para que cuadre
    x_im=min(x_map)# + ((1/10)*(max(x_map)-min(x_map)))
    y_im=min(y_map)# - (max(y_map)-min(y_map))/1000
    ab = AnnotationBbox(imagebox, (x_im,y_im), pad=0.5, box_alignment=(pos_ley_x,pos_ley_y)) #posicion del logo                         
    ax.add_artist(ab)
    

    #en las lineas de texto escribe sismicidad en tal año
    #tx1="SISMICIDAD REGISTRADA POR LA RSNC ENTRE: " +time_i +" y " + time_f +" (UTC)"
    tx1=" "
    ax.text(loc_x_txt,loc_y_txt, tx1, transform=ax.transAxes, fontsize=tam_text)
    
    #tamaño y posicion de la barra de escala
    esc_lon_min=lon_min+(loc_escala_x*(lon_max-lon_min)) #posicion de la barra
    esc_lon_max=lon_min+(loc_escala_x*(lon_max-lon_min))
    esc_lat_min=lat_min+(loc_escala_y*(lat_max-lat_min))
    esc_lat_max=lat_min+(loc_escala_y*(lat_max-lat_min))
    x_map,y_map=mapa([lon_min,lon_max],[lat_min,lat_max])
    long_barra=round((max(x_map)-min(x_map))/10000*3.5) #divide en 1000 para que de en km, luego en 10 para que de una decima parte del ancho del mapa, luego multiplica por 3.5 para obtener la barra ocupando el 35% del eje x (de la busqueda rectangular)
    long_barra= int(long_barra/10)*10 #para que la barra quede como 130 u no 133 por ejemplo
    
    # --------------------grafica Flecha dek Norte
    tam_letra_Norte=45
    porcentaje_tam_flecha=2.5 #porcentaje del tamaño de la flecha resecto a la longitud y del mapa
    if color_prof == "continuo":
        pos_x_flecha=0.9#posiciones en x y y de la flecha, estan normalizadas entre 0 y 1
    if color_prof == "discreto":
        pos_x_flecha=0.1#posiciones en x y y de la flecha, estan normalizadas entre 0 y 1
    pos_y_flecha=0.88
    overhang=1.8 #que tanto se deforma el triangulo de la flecha
    atLon = lon_min+(pos_x_flecha*(lon_max-lon_min))  # Longitude position for north arrow
    atLat = lat_min+(pos_y_flecha*(lat_max-lat_min)) # Latitude position for north arrow
    long_flecha=round( ((max(y_map)-min(y_map))/100) * porcentaje_tam_flecha ) #divide en 100 y multiplica por 2.5 para obtener la barra ocupando el 3% del eje y
    long_flecha= int(long_flecha/10)*10 #para que la barra quede como 130 u no 133 por ejemplo
    x, y = mapa(atLon, atLat)
    plt.arrow(x, y, 0, long_flecha, fc="k", ec="k", head_width=2.5*long_flecha, head_length=2*long_flecha, shape= 'right', head_starts_at_zero= False, alpha=0.8,linewidth=2, overhang=overhang,length_includes_head=True) 
    plt.arrow(x, y, 0, long_flecha, fc="None", ec="k", head_width=2.5*long_flecha, head_length=2*long_flecha, shape= 'left', head_starts_at_zero= False, alpha=0.8, linewidth=2, overhang=overhang,length_includes_head=True) 
    plt.text(x, y, "N", verticalalignment="top", horizontalalignment="center",fontsize=tam_letra_Norte)
     #---- fin grafica Flecha dek Norte

    if tipo_busqueda=="c":
        long_barra=radio
    #y configurar la escala como 3/10 de longitud del ancho del mapa

    if long_barra < 10: #corrige: al dar unas cordenadas muy apretadas en la busqueda, corrige el tamaño del mapa, para que sea mas grande y pueda traer la topografia, este minimo lo cosntruye para un long barra aprox de 8 km
        long_barra =10 

    mapa.drawmapscale(esc_lon_min, esc_lat_min, esc_lon_max, esc_lat_max, long_barra, barstyle='fancy', fontsize=22, zorder=40)
    
    #-----------------------
    return mapa, ax,lat_min,lat_max,lon_min,lon_max,long_barra,x_map,y_map



#def mapa_base(r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,alpha,exp_mag,color_prof,d_mag,y_mag,lat_pos_prof,incluir_mapa_ref,zoom_map_ref,tam_map_ref,caso_mapa,cmap):
def mapa_base(mapa, ax,long_barra,x_map,y_map,r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,alpha,exp_mag,color_prof,d_mag,y_mag,lat_pos_prof,incluir_mapa_ref,zoom_map_ref,tam_map_ref,caso_mapa,cmap):

    

    #mapa, ax,lat_min,lat_max,lon_min,lon_max,long_barra,x_map,y_map=mapa_crudo(r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,color_prof)
    
    if len(r)>0:
               

        r_sort=sorted(r, key=itemgetter(5),reverse=True) #(4) pone los mas profundos abajo, (5) pone los mas grandes abajo. True o False para que clasifique de mayor a menor o visceverza

        time=[]
        lats=[]
        lons=[]
        mags_map=[]
        mags=[]    
        profs=[]
        elats=[]
        elons=[]
        eprofs=[]
        rmss=[]
        gaps=[]
        
        
        for i in r_sort:

            fecha=str(i[0])
            horas=str(i[1])
            utc_fecha=UTCDateTime(fecha+horas).datetime #-timedelta(hours=5) #para hora local restar 5 horas
            time.append(utc_fecha)

            lat=float(i[2])
            lon=float(i[3])
            prof=float(i[4])                
            mag=float(i[5])
            lats.append(lat)
            lons.append(lon)
            profs.append(prof)        
            mags_map.append(( ( (mag/8)*(8-FAM)) +FAM )**exp_mag) 
            mags.append(mag)

            if op == "bd":

                if float(i[6]) > 10 or i[6]=="": #retira valores anomalos de la grafica
                    rmss.append(1) #se pone un valor para poder graficar, estos datos son anomalos, si se pone un valor grande, se pierde la grafica
                else:
                    rmss.append(float(i[6]))                
            
                if float(i[7]) > 360 or i[7]=="": #retira valores anomalos de la grafica
                    gaps.append(180) #se pone un valor para poder graficar, estos datos son anomalos, si se pone un valor grande, se pierde la grafica
                else:                    
                    gaps.append(float(i[7]))

                elats.append(float(i[8]))
                elons.append(float(i[9]))
                eprofs.append(float(i[10]))                
            
            if op == "ap" and len(i)>10: #en el elemnto 6 del archivo propio va el epicentro, en el elemento 11 esta eprof 

                if i[7].isdigit()==False:
                    rmss.append(0)
                elif float(i[7]) < 10: #retira valores anomalos de la grafica
                    rmss.append(float(i[7]))
                if i[8].isdigit()==False:
                    gaps.append(0)
                elif float(i[8]) < 360: #retira valores anomalos de la grafica
                    gaps.append(float(i[8]))
                if i[9].isdigit()==False:
                    elats.append(0)
                else:
                    elats.append(float(i[9]))
                if i[10].isdigit()==False:
                    elons.append(0)
                else:
                    elons.append(float(i[10]))
                if i[11].isdigit()==False:
                    eprofs.append(0)
                else:
                    eprofs.append(float(i[11]))
            
        x, y = mapa(lons, lats) #hace una proyeccion a cordenadas planas con origen (0,0) arbitrario en la esquina inferior izquierda del mapa
        
                
        if caso_mapa=="SISMICIDAD":            
            mapa, ax = mapa_profundidad(r_sort,long_barra,lon_min,lon_max,loc_x_txt,pos_mag_x,loc_y_txt,mapa,x,y,lat_pos_prof,lat_min,alpha,profs,mags_map,ax,color_prof,tam_text,y_mag,cmap)
        
        if caso_mapa =="GAP":            
            mapa, ax= mapa_gap(y_mag,long_barra,mapa,ax,x,y,gaps,mags_map,alpha,loc_x_txt,pos_mag_x,loc_y_txt,tam_text,cmap)
        
        if caso_mapa =="RMS":
            mapa, ax= mapa_rms(y_mag,long_barra,mapa,ax,x,y,rmss,mags_map,alpha,loc_x_txt,pos_mag_x,loc_y_txt,tam_text,cmap)

        if caso_mapa =="TIEMPO":
            mapa, ax= mapa_time(y_mag,long_barra,mapa,ax,x,y,time,mags_map,alpha,loc_x_txt,pos_mag_x,loc_y_txt,tam_text,cmap)



        mapa=leyenda_magnitudes(lon_min,lon_max,lat_min,pos_mag_x,d_mag,long_barra,mags,FAM,exp_mag,mapa,dir_ley_mag,y_mag)
        #---------------------leyenda de magnitudes
        

    else:
        lons=[]
        lats=[]
        profs=[]
        mags=[]

    
    #return mapa, ax, lons, lats, profs, mags, lat_min,lat_max,lon_min,lon_max



def mapa_calor(mapa, ax,long_barra,x_map,y_map,r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,prof_min,prof_max,
                time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,
                tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,
                tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,color_prof):


    if len(r)>1:
        
        lats=[]
        lons=[]
        profs=[]                
        
        for i in r:
           

            lat=float(i[2])
            lon=float(i[3])
            prof=float(i[4])       
            
            lats.append(lat)
            lons.append(lon)
            profs.append(prof)        
            
        op2="acumulado" #"prom_dia", "acumulado"
        escala_log10=True #solo sirve para la version sismos acumulados

        #print("Generando mapa de calor por "+op2+"...")
        #mapa, ax,lat_min,lat_max,lon_min,lon_max,long_barra,x_map,y_map=mapa_crudo(r,lon_centro,lat_centro,lat_min,lat_max,lon_min,lon_max,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda-0.04,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,color_prof)
        
        profs=[p*(-1) for p in profs]

        suavi=1.4
        res=15 # cambiar segun nececidad , entre mayor sea, menor el tamaño de los bins, recomendado para mapas pequeños 5, para mapas grandes 10
        db = int(long_barra/res) 
        if db==0:
            db=1
        print("Tamaño del bin (km): ",db)
        
        #-----Esta predeterminado para que calcule densidades usando celdas que comienzan desde los maximos y minimos del mapa y no dependa de la localizacion de los eventos, Si se desea fijar 
        #el calculo entre los maximos y minimos de posision de los eventos activar la siguiente  bloque de codigo
        '''
        lon_min=min(lons)
        lon_max=max(lons)
        lat_min=min(lats)
        lat_max=max(lats)        
        prof_min=min(profs)
        prof_max=max(profs)
        '''
        prof_min=min(profs)-(0.1*(prof_max-prof_min))
        prof_max=max(profs)+(0.1*(prof_max-prof_min))        
        #-----------
        

        lon_bins=np.arange(lon_min,lon_max,round(db/110,2)) 
        lat_bins=np.arange(lat_min,lat_max,round(db/110,2))
        prof_bins=np.arange(prof_min,prof_max,round(db,2))  
        
        porcentaje_minimo=0.0005 #  densidades menores al a este porcentaje del maximo no se muestran en el mapa de calor

        if op2=="prom_dia":
            
            del_t_dia=(UTCDateTime(time_f)-UTCDateTime(time_i))/3600/24
            #print(del_t_dia) numero de dias

            density,lon_bins, lat_bins = np.histogram2d(lons, lats, [lon_bins, lat_bins])
            density=density.T/del_t_dia
            
            
            density_lonprof,lon_binsP, prof_bins = np.histogram2d(lons, profs, [lon_bins, prof_bins])
            density_lonprof=density_lonprof.T/del_t_dia
            density_latprof, lat_binsP, prof_bins = np.histogram2d(lats, profs, [lat_bins, prof_bins])
            density_latprof=density_latprof.T/del_t_dia
            
            density=gaussian_filter(density, suavi) #para suavisar el mapa de calor, el parametro despues del array  incrementa la suavizada
            density_lonprof=gaussian_filter(density_lonprof, 1) #para suavisar el mapa de calor, el parametro despues del array  incrementa la suavizada
            density_latprof=gaussian_filter(density_latprof, 1) #para suavisar el mapa de calor, el parametro despues del array  incrementa la suavizada

            dens_minima_a_mostrar=porcentaje_minimo*np.amax(density)
            v_min=dens_minima_a_mostrar
        
        if op2=="acumulado":

            density,lon_bins, lat_bins = np.histogram2d(lons, lats, [lon_bins, lat_bins])
            density=density.T
            
                            
            density_lonprof,lon_binsP, prof_bins = np.histogram2d(lons, profs, [lon_bins, prof_bins])
            density_lonprof=density_lonprof.T
            density_latprof, lat_binsP, prof_bins = np.histogram2d(lats, profs, [lat_bins, prof_bins])
            density_latprof=density_latprof.T
            dens_minima_a_mostrar=int(porcentaje_minimo*np.amax(density))+1 #densidades menores porcentaje minimo no se muestran

            density=gaussian_filter(density, suavi) #para suavisar el mapa de calor, el parametro despues del array  incrementa la suavizada
            density_lonprof=gaussian_filter(density_lonprof, 0.5) #para suavisar el mapa de calor, el parametro despues del array  incrementa la suavizada
            density_latprof=gaussian_filter(density_latprof, 0.5) #para suavisar el mapa de calor, el parametro despues del array  incrementa la suavizada

            if escala_log10==True:
                density=np.log10(density)
                density_lonprof=np.log10(density_lonprof)
                density_latprof=np.log10(density_latprof)
                dens_minima_a_mostrar=int(porcentaje_minimo*np.amax(density))
    
                
            v_min=dens_minima_a_mostrar

        lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
        xs, ys = mapa(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh

        #plt.pcolormesh(xs, ys, density, cmap="seismic",alpha=0.4,shading="gourand",linewidth=0)
        
        #truco para generar alpha cero a valores indeseados
        my_norm = matplotlib.colors.Normalize(vmin=v_min, vmax=np.amax(density), clip=False)
        #my_norm = matplotlib.colors.Normalize(vmin=np.amin(density), vmax=np.amax(density), clip=False)
        my_cmap = copy(cm.get_cmap('jet')) # generando paleta propia 
        #'YlOrRd' claro a rojo
        #"jet" azul a rojo
        #"cool" azul claro a lila
        #"viridis" morado a amarillo
        my_cmap.set_under(alpha=0) # los valores bajo el minimo de normalizacion les pone alpha cero 
        #my_cmap.set_over(alpha=0)   # los valores sobre el minimo de normalizacion les pone alpha cero 



        extent=[xs.min(),xs.max(),ys.min(),ys.max()]
        #plt.contour(density,extent=extent, cmap="seismic",linewidths=1,vmin=20, levels =100,alpha=0.8)
        alpha_calor=0.4
        plt.imshow(density, extent=extent, cmap= my_cmap, norm=my_norm,interpolation='nearest',alpha=alpha_calor,origin = 'lower')
        
        

        cbar = plt.colorbar(orientation='vertical', shrink=0.6, aspect=30, fraction=0.2,pad=0.01)
        
        if op2=="prom_dia":
            cbar.set_label("\nSismos/dia\n"+time_i+" / "+time_f,size=40)
        if op2=="acumulado":
            if escala_log10==True:
                #cbar.set_label("\nNumero de sismos (10^)\n"+time_i+" / "+time_f,size=40)
                cbar.set_label("\nNumero de sismos (10^)",size=40)
            else:
                #cbar.set_label("\nNumero de sismos\n"+time_i+" / "+time_f,size=40)
                cbar.set_label("\nNumero de sismos",size=40)
        cbar.ax.tick_params(labelsize=30)        

        plt.savefig('MAPA_CALOR.png', bbox_inches='tight')          
        
        #--------------graficando mapa calor de perfil-
        #
        my_norm_lon_prof = matplotlib.colors.Normalize(vmin=v_min, vmax=np.amax(density_lonprof), clip=False)
        my_norm_lat_prof = matplotlib.colors.Normalize(vmin=v_min, vmax=np.amax(density_latprof), clip=False)
        
        #x- prof
        

        fig, ax = plt.subplots(figsize=(33, 45)) #tamaño de la ventana del mapa
        lon_bins_2d, prof_bins_2d = np.meshgrid(lon_bins, prof_bins)
        #xs, ys = mapa(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh
        plt.pcolormesh(lon_bins_2d, prof_bins_2d, density_lonprof, cmap=my_cmap, norm=my_norm_lon_prof)
        plt.ylim(prof_min, prof_max)
        cbar = plt.colorbar(orientation='horizontal', shrink=0.8, aspect=20, fraction=0.2,pad=0.05)
        
        if op2=="prom_dia":
            cbar.set_label("Sismos/dia\n"+time_i+" / "+time_f,size=40)
        if op2=="acumulado":
            if escala_log10==True:
                #cbar.set_label("\nNumero de sismos (10^)\n"+time_i+" / "+time_f,size=40)
                cbar.set_label("\nNumero de sismos (10^)",size=40)
            else:
                #cbar.set_label("\nNumero de sismos\n"+time_i+" / "+time_f,size=40)
                cbar.set_label("\nNumero de sismos",size=40)
        
        plt.title("Lon vs Z\n", fontsize=50)
        plt.xlabel('lon', color='k', fontsize=50)
        plt.ylabel('Z', color='k', fontsize=50)
        plt.xticks(fontsize=30)    
        plt.yticks(fontsize=30)
        cbar.ax.tick_params(labelsize=30)
        
        plt.savefig('lon_vs_z.png', bbox_inches='tight')
        

        #y- prof
        fig, ax = plt.subplots(figsize=(33, 45)) #tamaño de la ventana del mapa
        lat_bins_2d, prof_bins_2d = np.meshgrid(lat_bins, prof_bins)
        #xs, ys = mapa(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh

        plt.pcolormesh(lat_bins_2d, prof_bins_2d, density_latprof, cmap=my_cmap, norm=my_norm_lat_prof)
        plt.ylim(prof_min, prof_max)
        cbar = plt.colorbar(orientation='horizontal', shrink=0.625, aspect=20, fraction=0.2,pad=0.05)
        if op2=="prom_dia":
            cbar.set_label("Sismos/dia\n"+time_i+" / "+time_f,size=40)
        if op2=="acumulado":
            if escala_log10==True:
                cbar.set_label("\nNumero de sismos (10^)\n"+time_i+" / "+time_f,size=40)
            else:
                cbar.set_label("\nNumero de sismos\n"+time_i+" / "+time_f,size=40)
        
        plt.title("Lat vs Z\n", fontsize=50)
        plt.xlabel('lat', color='k', fontsize=50)
        plt.ylabel('Z', color='k', fontsize=50)
        plt.xticks(fontsize=30)    
        plt.yticks(fontsize=30)
        cbar.ax.tick_params(labelsize=30)
        
        plt.savefig('lat_vs_z.png', bbox_inches='tight')

        
        

        return mapa_calor    
    else:
        print("El reporte de sismicidad tiene menos de 10 sismos, no hay para que hacer mapa de calor")


def grafica_estaciones(lat_min,lat_max,lon_min,lon_max,lat_centro,lon_centro, radio,tipo_busqueda, mapa,ax,TX,TY,tam_text_est,tam_text,tam_triangulo,pos_triag_x,pos_triag_y,pos_tx_est_x,pos_tx_est_y,tam_triang_leyen,ver_temporales,ver_permanentes,imprimir_cod_estacion):
    
    #recibe de la funcion mapa_base lat_min,lat_max,lon_min,lon_max que fueron establecidas para el mapa, 
    #en el que se considera el ajuste de esos limites de acuerdo a las condiciones iniciales de busqueda (ver el inicio de esa funcion)

    print("Graficando estaciones sismologicas...")                                    
    
    lat_min=lat_min+((lat_max-lat_min)/8)
    #reajusta al tamaño que excedio en la funcion mapa (para pintar las leyendas) aqui quita nuevamente ese rectangulo
    #adional de latitud para que no se sobrepongan las estaciones sobre las leyendas
    
    filtrar_est=False#poner estaciones que se esta interesado en graficar
    estaciones_interes=["CAP2","PTB","URE","BRJC","JAMC","ZAR"]
    #----------------------------------------------------------------------------
    
    archivo_est_local="/opt/rutinas/mapas_histogramas/catSeis_ests_leyens/estaciones.csv" #archivo obtenido de la base de datos. 
    est=open(archivo_est_local, 'r')
    ests_perm=[]    
    ests_temp=[]    

    for e in est:
        vec_e=e.split(",")        
        tipo_est=vec_e[3][0:3]

        if tipo_est=="per":            
            ests_perm.append([vec_e[0],vec_e[1],vec_e[2]])
                        
        if tipo_est=="tem":
            if "VMM" in vec_e[0] or "DRL" in vec_e[0] or "LL" in vec_e[0] or "PGA" in vec_e[0] or "ACH" in vec_e[0]: #grafica temporales del VMM, DRL, LL, PGA o ACH, cambiar segun el caso    
                ests_temp.append([vec_e[0],vec_e[1],vec_e[2]])

    
    names=[]
    lats=[]
    lons=[]

    if ver_permanentes==True:

        for e in ests_perm:
            nombre_est=e[0]
            lat_est=float(e[1])
            lon_est=float(e[2])

            if filtrar_est==True:
                if lat_min<lat_est<lat_max and lon_min<lon_est<lon_max and nombre_est in estaciones_interes:
                    lats.append(lat_est)
                    lons.append(lon_est)
                    names.append(nombre_est)
            else:
                if lat_min<lat_est<lat_max and lon_min<lon_est<lon_max:
                    lats.append(lat_est)
                    lons.append(lon_est)
                    names.append(nombre_est)

    if ver_temporales==True:
        
        for e in ests_temp:
            nombre_est=e[0]
            lat_est=float(e[1])
            lon_est=float(e[2])

            if filtrar_est==True:
                if lat_min<lat_est<lat_max and lon_min<lon_est<lon_max and nombre_est in estaciones_interes:
                    lats.append(lat_est)
                    lons.append(lon_est)
                    names.append(nombre_est)
            else:
                if lat_min<lat_est<lat_max and lon_min<lon_est<lon_max:
                    lats.append(lat_est)
                    lons.append(lon_est)
                    names.append(nombre_est)

    if len(lons) > 0 and len(lats) > 0:
        print("\nSe encontraron "+str(len(lats))+" Estaciones")
        x,y=mapa(lons, lats) 
        est = mapa.scatter(x, y, marker='^',color="#0000ff", edgecolor='k', s=tam_triangulo*1000, zorder=10)

        
        i=0
        for n in names:
            ax.text(x[i]+(TX*1000),y[i]+(TY*1000), n, fontsize=tam_text_est, horizontalalignment='center', color="k", style='italic',zorder=1000)
            i=i+1

        if imprimir_cod_estacion==True:    
            tx_est="ESTACION"             
            ax.text(pos_tx_est_x,pos_tx_est_y, tx_est, transform=ax.transAxes, fontsize=tam_text, horizontalalignment='center')
        
        ax.plot(pos_triag_x,pos_triag_y,"^",color="#0000ff", markeredgecolor='k',transform=ax.transAxes, markersize=tam_triang_leyen)
        
    else:
        if filtrar_est==True:
            print("No se encontro ninguna estacion (cambie las estaciones de interes) ")
        if filtrar_est==False:
            print("No se encontro ninguna estacion en el area seleccionada")    


def restituir_geografia(lat_min_o,lat_max_o,lon_min_o,lon_max_o):
    lat_min=lat_min_o
    lat_max=lat_max_o
    lon_min=lon_min_o
    lon_max=lon_max_o

    return lat_min,lat_max,lon_min,lon_max


def generacion_mapas(bool_mapas_deseados,mapas_posibles,perfiles,incluir_perfil,long_barra,x_map_c,y_map_c,r,lon_centro,lat_centro,
    lat_min_o,lat_max_o,lon_min_o,lon_max_o,prof_max,prof_min,time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,
    zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,
    FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,
    color_prof,ancho,lat_min_c,lat_max_c,lon_min_c,lon_max_c,alpha,exp_mag,d_mag,y_mag,lat_pos_prof,cmap):
                         

    i=0                
    for caso_mapa in mapas_posibles:

        if bool_mapas_deseados[i]==True:

            nombre_archivo_leer="ax_c"
            with open(nombre_archivo_leer, "rb") as fp:
                ax_c_c = pickle.load(fp)

            nombre_archivo_leer="mapa_c"
            with open(nombre_archivo_leer, "rb") as fp:
                mapa_c_c = pickle.load(fp)

            print("Obteniendo recursos para mapa de "+caso_mapa)
            
            if caso_mapa == "GAP" or caso_mapa == "RMS" or caso_mapa == "TIEMPO" or caso_mapa == "CALOR":
                    color_prof="continuo" #estos mapas estan configurados para color continuo, volver a definir la variable ayuda a posicionar mejor el mapa de refrencia


            if caso_mapa == "CALOR":#---- ---genera mapa de calor, acumulacion espacial de sismicidad

                mapa_calor(mapa_c_c, ax_c_c,long_barra,x_map_c,y_map_c,r,lon_centro,lat_centro,lat_min_c,lat_max_c,lon_min_c,lon_max_c,-prof_max,-prof_min,
                            time_i,time_f,loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,
                            tipo_mapa,incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,
                            tam_tex_falla,indica_centro_busqueda,op,en_negrilla,incluir_mapa_ref,zoom_map_ref,tam_map_ref,color_prof)
                plt.close()
                 
            else:

                if incluir_perfil==True and caso_mapa == "SISMICIDAD":
                    for p in perfiles:
                        lon_1=p[0]
                        lon_2=p[1]
                        lat_1=p[2]
                        lat_2=p[3]
                        lp=p[4]
                        mapa_c_c,ax_c_c=genera_perfil(mapa_c_c,ax_c_c,lon_1,lat_1,lon_2,lat_2,r,FAM,exp_mag,alpha,ancho,color_prof,cmap,lp)

                
                mapa_base(mapa_c_c, ax_c_c,long_barra,x_map_c,y_map_c,r,lon_centro,lat_centro,lat_min_c,lat_max_c,lon_min_c,lon_max_c,time_i,time_f,
                        loc_x_txt,loc_y_txt,loc_escala_x,loc_escala_y,tam_text,zoom_leyenda,pos_mag_x,pos_ley_x,pos_ley_y,dir_ley_mag,tipo_mapa,
                        incluir_fallas_sgc,incluir_fallas_movmasa,incluir_rectangulo_leyenda,FAM,tipo_busqueda,radio,ancho_fallas,tam_tex_falla,
                        indica_centro_busqueda,op,en_negrilla,alpha,exp_mag,color_prof,d_mag,y_mag,lat_pos_prof,incluir_mapa_ref,zoom_map_ref,
                        tam_map_ref,caso_mapa,cmap) #h leyenda de magnitud horizontal
                
                plt.savefig('MAPA_'+caso_mapa+'.png', bbox_inches='tight')  #siempre incluir esta linea al final del codigo cuando se llama la fucnion mapa_base para que guarde la imagen en la carpeta
            
            if incluir_mapa_ref==True:          
                agrega_mapa_ref('MAPA_'+caso_mapa+'.png','MAPA_REF.png',zoom_map_ref,color_prof)
            
            plt.close()

        i=i+1

