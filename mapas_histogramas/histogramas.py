import csv
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import matplotlib.ticker as tck
import numpy as np
import pandas as pd
from obspy import UTCDateTime #importante para el condicional de tiempo en la fuincion selec_catalogo_seisan
import MySQLdb
from operator import itemgetter #tiene la funcion sort que sevira para clasificar los sismos de acurdo a sus profundidades o magnitudes
from itertools import groupby
from obspy.geodetics import kilometer2degrees
import sys
from datetime import date, datetime, timedelta
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, DateFormatter, date2num, num2date
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import copy
import math
from collections import Counter
import logging
import utm
import geopy.distance
from windrose import WindroseAxes
import matplotlib.cm as cm


sns.set(style="darkgrid")


def impime_mags_en_histo(t,cuentas,r,tik,ax,periodo):

    if tik > 20:
        tik=20
    #-------genera un vector de magnitudes con cada tiempo (en esta opcion de conteo por dia), 
    #sera usado para imprimir en la grafica la magnitud cuando es menor a dos meses
    time_mags=[]
    for t_unico in t:    
        mags=[]
        for i in r:
            if periodo == "d":
                fecha=str(i[0])                        
                utc_fecha=UTCDateTime(fecha+"00:00:00").datetime    
            if periodo == "h":
                fecha=str(i[0])    
                hora=str(i[1].split(":")[0])                    
                utc_fecha=UTCDateTime(fecha+hora+":00:00").datetime    

            if t_unico == utc_fecha:
                mag=float(i[5])
                mags.append(mag)
        
        time_mags.append([t_unico,mags])
    #-------------
    registro_y=[]
    #delta=int(tik/5)
    delta=(0.06*max(cuentas))+1

    max_num_sis_imp=6
    
    c=0
    for t_m in time_mags:
        t_unico=t_m[0]
        mags=t_m[1]
        i=2
        
        if cuentas[c] <= max_num_sis_imp:
            
            for m in mags:
            
                if m < 4:                        
                    ax.text(t_unico, cuentas[c] + (delta*i),"M:"+str(m),horizontalalignment='center',fontsize=tik-3,rotation=60)

                if m >=4:                        
                    ax.text(t_unico, cuentas[c] + (delta*i),"*M:"+str(m),horizontalalignment='center',fontsize=tik-2,rotation=60,color="#ff0000")
                i=i+1
            
        if cuentas[c] > max_num_sis_imp:
            
            mags_sort=sorted(mags)

            for m in mags_sort[::-1][0:max_num_sis_imp]: #imprime los (max_num_sis_imp) eventos mas grandes de ese dia ()
                if  m < 4:                        
                    ax.text(t_unico, cuentas[c] + (delta*i),"M:"+str(m),horizontalalignment='center',fontsize=tik-3,rotation=60)
                    i=i+1
                if m >=4:                        
                    ax.text(t_unico, cuentas[c] + (delta*i),"*M:"+str(m),horizontalalignment='center',fontsize=tik-2,rotation=60,color="#ff0000")
                    i=i+1
            ax.text(t_unico, cuentas[c] + (delta*i),"...",horizontalalignment='center',fontsize=tik,rotation=60)
                    
                
        registro_y.append(cuentas[c] + (delta*i))            
        
        c=c+1
    
    
    ax.set_ylim( 0, 1.1*max(registro_y) )

def genera_cuerpo_grafica(delta_time,time_i,time_f):
    
    utc_time_i = UTCDateTime(time_i).datetime
    utc_time_f = UTCDateTime(time_f).datetime    

    fig = plt.figure(figsize=(21., 11.))
    ax = fig.add_subplot(111)
    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')    
    
    return ax,plt,utc_time_i,utc_time_f

def genera_ejes(ax,delta_time,utc_time_i,utc_time_f):

    if  delta_time < 5:
        interval_minor=int( round (0.5*delta_time ,0) )
        interval_major=int( round (0.15*delta_time ,0) )
        if interval_minor<1:
            interval_minor=1
        if interval_major<1:
           interval_major=1 

        delta_ini_fin_graf=int(1.1*delta_time)
        ax.xaxis.set_minor_locator( HourLocator(interval = interval_minor))
        ax.xaxis.set_minor_formatter( DateFormatter('%d-%H:%M') )
        ax.xaxis.set_major_locator( DayLocator(interval = interval_major))
        ax.xaxis.set_major_formatter( DateFormatter('\n\n\n\n%d-%b-%Y'))
        ini_graf=utc_time_i-timedelta(hours=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(hours=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('Tiempo (dia-hora/dia-mes-año)', color='k', fontsize=18)
    
    if  5 <= delta_time < 30:
        
        interval_major=int( round (0.4*delta_time ,0) )
        delta_ini_fin_graf=int(0.9*delta_time)
        ax.xaxis.set_minor_locator( DayLocator(interval = 1))
        ax.xaxis.set_minor_formatter( DateFormatter('%d-%b') )
        ax.xaxis.set_major_locator( DayLocator(interval = interval_major))
        ax.xaxis.set_major_formatter( DateFormatter('\n\n\n%d-%b-%Y\n'))
        ini_graf=utc_time_i-timedelta(hours=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(hours=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('Tiempo (dia-mes/dia-mes-año)', color='k', fontsize=18)

    if  30 <= delta_time < 50:
        interval_minor=int( round (0.03*delta_time ,0) )
        interval_major=int( round (0.5*delta_time ,0) )
        delta_ini_fin_graf=int(0.9*delta_time)
        ax.xaxis.set_minor_locator( DayLocator(interval = interval_minor))
        ax.xaxis.set_minor_formatter( DateFormatter('%d-%b') )
        ax.xaxis.set_major_locator( DayLocator(interval = interval_major))
        ax.xaxis.set_major_formatter( DateFormatter('\n\n\n%d-%b-%Y\n'))
        ini_graf=utc_time_i-timedelta(hours=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(hours=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('Tiempo (dia-mes/dia-mes-año)', color='k', fontsize=18)

    if  50 <= delta_time < 450:
        interval_minor=int(delta_time/62)+1
        interval_major=int(delta_time/320)+1
        delta_ini_fin_graf=int(delta_time/30)+1
        ax.xaxis.set_minor_locator( DayLocator(interval = interval_minor))
        ax.xaxis.set_minor_formatter( DateFormatter('%d-%b') )
        ax.xaxis.set_major_locator( MonthLocator(interval = interval_major) )
        ax.xaxis.set_major_formatter( DateFormatter('\n\n\n%b-%Y\n'))
        ini_graf=utc_time_i-timedelta(days=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(days=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('\nTiempo (dia/mes-año)', color='k', fontsize=18)

    if 450 <= delta_time < 6000: 
        interval_minor=int(delta_time/1700)+1
        interval_major=int(delta_time/12000)+1
        delta_ini_fin_graf=int(delta_time/30)+1
        ax.xaxis.set_minor_locator( MonthLocator(interval = interval_minor))
        ax.xaxis.set_minor_formatter( DateFormatter('%b') )
        ax.xaxis.set_major_locator( YearLocator(interval_major) )
        ax.xaxis.set_major_formatter( DateFormatter('\n\n%Y\n'))
        ini_graf=utc_time_i-timedelta(days=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(days=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('Tiempo (mes/año)', color='k', fontsize=18)

    if 6000 <= delta_time < 12000: #
        #print("aqui")
        interval_minor=int(delta_time/1600)+1
        interval_major=int(delta_time/10000)+1
        delta_ini_fin_graf=int(delta_time/30)+1
        ax.xaxis.set_minor_locator( MonthLocator(interval = interval_minor))
        ax.xaxis.set_minor_formatter( DateFormatter('%b') )
        ax.xaxis.set_major_locator( YearLocator(interval_major) )
        ax.xaxis.set_major_formatter( DateFormatter('\n\n%Y\n'))
        ini_graf=utc_time_i-timedelta(days=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(days=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('Tiempo (mes/año)', color='k', fontsize=18)

    if delta_time >= 12000: #
        interval_minor=int(delta_time/1550)+1
        interval_major=int(delta_time/5000)+1
        delta_ini_fin_graf=int(delta_time/30)+1
        ax.xaxis.set_minor_locator( MonthLocator(interval = interval_minor))
        ax.xaxis.set_minor_formatter( DateFormatter('%b') )
        ax.xaxis.set_major_locator( YearLocator(interval_major) )
        ax.xaxis.set_major_formatter( DateFormatter('\n\n%Y\n'))
        ini_graf=utc_time_i-timedelta(days=delta_ini_fin_graf)
        fin_graf=utc_time_f+timedelta(days=delta_ini_fin_graf)
        ax.set_xlim( ini_graf, fin_graf)
        ax.set_xlabel('Tiempo (mes/año)', color='k', fontsize=18)
        
        #ax.set_xlim( min(vec_time)-timedelta(days=30), max(vec_time)+timedelta(days=30))  
    return ini_graf,fin_graf

def histo_conteo_repeticion_mes(time,delta_time):

    meses=[]
    for t in time:        
        mes=t.month                
        meses.append(mes)

    fig = plt.figure(figsize=(21., 11.))
    ax = fig.add_subplot(111)
    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')  
    
    bins = np.arange(1,14, 1)
    df, bins = np.histogram(meses, bins)
    ax.vlines(bins[:-1], [0],df, color="#0000ff",lw=30,alpha=0.9)

    #lis_meses=["ene","feb","mar","abr","may","jun","jul","ago","sep","oct","nov","dic"]
    #ax.set_xticklabels(lis_meses)
    ax.set_xlabel("Mes del año", color='k', fontsize=20)
    ax.set_ylabel("Conteo", color='k', fontsize=20)
        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    plt.title("Numero de sismos/mes del año", color='k', fontsize=14)
    plt.savefig('Num_sismos_mes_del_an.png', bbox_inches='tight')
    plt.close()

def histo_conteo_repeticion_dia(time,delta_time):

    dias=[]
    for t in time:        
        dia=t.day                
        dias.append(dia)

    fig = plt.figure(figsize=(21., 11.))
    ax = fig.add_subplot(111)
    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')  
    
    bins = np.arange(1,33, 1)
    df, bins = np.histogram(dias, bins)
    ax.vlines(bins[:-1], [0],df, color="#0000ff",lw=12,alpha=0.9)
    
    ax.set_xlabel("Dia del mes", color='k', fontsize=20)
    ax.set_ylabel("Conteo", color='k', fontsize=20)
        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    plt.title("Numero de sismos/dia del mes", color='k', fontsize=14)
    plt.savefig('Num_sismos_dia_del_mes.png', bbox_inches='tight')
    plt.close()



def histo_conteo_horas(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time,r):
    #sismicidad hora

    periodo="h"
    conteo = df_tiempo.groupby([df_tiempo['tiempo'].dt.to_period(periodo)]).count() #obtiene un dataframe con el conteo dependiendo del periodo p: h:horas, d: dias, m: meses, y: años
    t=conteo.index.to_timestamp() #obtiene un vector de tiempos en formato datatime
    cuentas=(conteo['tiempo'].tolist()) #obtiene un vector con el conteo

    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)

    if delta_time==0:
        ab=30
    else:
        ab=int(1.4+(440/(24*delta_time)))
    
    if ab==0:
        ab=1
    
    tik=int(ab/2)
    if tik<14:
        tik=14

    
    ax.vlines(t, [0],cuentas, lw=ab, color="#0000ff") 
    ax.set_ylabel('No. de eventos', color='k', fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)

    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=16,ha="right")
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=18)

    ax.yaxis.grid(True, which='major')
        
    #se grafica un texto sobre cada barra del histograma que es el numero de sismos que representa la barra
    a=0;

    for i in cuentas:

        if i>0 and delta_time<6: #no imprime si el numero de sismos es cero para cierto año
            
            tx = str(i)
            ax.text(t[a], cuentas[a]+(max(cuentas)*0.01), tx, horizontalalignment='center',fontsize=tik) #sobre 50 para que imprima el texto muy cercano al tope de la barra
        a=a+1

    if delta_time < 10:
        impime_mags_en_histo(t,cuentas,r,tik,ax,periodo)
    

def histo_conteo_dias(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time,r):
    

    periodo="d"
    conteo = df_tiempo.groupby([df_tiempo['tiempo'].dt.to_period(periodo)]).count() #obtiene un dataframe con el conteo dependiendo del periodo p: h:horas, d: dias, m: meses, y: años
    t=conteo.index.to_timestamp() #obtiene un vector de tiempos en formato datatime
    cuentas=(conteo['tiempo'].tolist()) #obtiene un vector con el conteo

    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)
        
    if delta_time==0:
        ab=30
    else:
        ab=int ( round ( (1.4+(400/delta_time)), 1) )
    
    if ab==0:
        ab=1
    
    tik=int ( round (0.32*ab ,0) )
    if tik<14:
        tik=14

    
    ax.vlines(t, [0],cuentas, lw=ab, color="#0000ff") 
    ax.set_ylabel('No. de eventos', color='k', fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)

    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=16)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=18)

    ax.yaxis.grid(True, which='major')
        
    #se grafica un texto sobre cada barra del histograma que es el numero de sismos que representa la barra
    

    if delta_time<100: #si hay menos de 100 dias imprime las cuentas

        a=0
    
        for i in cuentas:

            if i>0: #no imprime si el numero de sismos es cero 
                
                tx = str(i)
                ax.text(t[a], cuentas[a]+(max(cuentas)*0.01), tx, horizontalalignment='center',fontsize=tik) #sobre 50 para que imprima el texto muy cercano al tope de la barra
            a=a+1

    
    if delta_time < 62: #si hay dos meses pinta las magnitudes en el histograma
        impime_mags_en_histo(t,cuentas,r,tik,ax,periodo)
    
    '''
    arc=open("sismicidad_diaria.csv","w")
    for i in range(0,len(t)):
        arc.writelines(str(t[i])+","+str(cuentas[i])+"\n")
    arc.close()
    '''

def histo_conteo_meses(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time):
    
    
    #sismicidad mes
    conteo = df_tiempo.groupby([df_tiempo['tiempo'].dt.to_period("m")]).count() #obtiene un dataframe con el conteo dependiendo del periodo p: h:horas, d: dias, m: meses, y: años
    t=conteo.index.to_timestamp() #obtiene un vector de tiempos en formato datatime
    cuentas=(conteo['tiempo'].tolist()) #obtiene un vector con el conteo
    
    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)
        
    if delta_time==0:
        ab=30
    else:
        ab=int(1.4+(440/(delta_time/30)))
    
    if ab==0:
        ab=1
    
    tik=int(ab/2)
    if tik<14:
        tik=14

    
    ax.vlines(t, [0],cuentas, lw=ab, color="#0000ff") 
    ax.set_ylabel('No. de eventos', color='k', fontsize=20)
    
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)

    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=16)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=18)

    ax.yaxis.grid(True, which='major')
        
    #se grafica un texto sobre cada barra del histograma que es el numero de sismos que representa la barra
    a=0;

    for i in cuentas:

        if i>0 and delta_time<1800: #no imprime si el numero de sismos es cero para cierto año
            
            tx = str(i)
            ax.text(t[a], cuentas[a]+(max(cuentas)*0.01), tx, horizontalalignment='center',fontsize=tik) #sobre 50 para que imprima el texto muy cercano al tope de la barra
        a=a+1    
    
def histo_conteo_ans(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time):
    #sismicidad mes
    conteo = df_tiempo.groupby([df_tiempo['tiempo'].dt.to_period("y")]).count() #obtiene un dataframe con el conteo dependiendo del periodo p: h:horas, d: dias, m: meses, y: años
    t=conteo.index.to_timestamp() #obtiene un vector de tiempos en formato datatime
    cuentas=(conteo['tiempo'].tolist()) #obtiene un vector con el conteo

    #ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)
    
    interval_major=1#int(delta_time/8000)+1
    if delta_time > 8000:
        delta_ini_fin_graf=450#int(delta_time/5)+1
    if delta_time <= 8000:
        delta_ini_fin_graf=300#int(delta_time/5)+1

    ax.xaxis.set_major_locator( YearLocator(interval_major) )
    ax.xaxis.set_major_formatter( DateFormatter('%Y') )    
    ini_graf=utc_time_i-timedelta(days=delta_ini_fin_graf)
    fin_graf=utc_time_f+timedelta(days=delta_ini_fin_graf)
    ax.set_xlim( ini_graf, fin_graf)
    ax.set_xlabel('Tiempo (año)', color='k', fontsize=18)
    
    if delta_time==0:
        ab=30
    else:
        ab=int(25+(220/(delta_time/365)))
    if ab==0:
        ab=1
    
    tik=int(ab/5)
    if tik<14:
        tik=14

    
    ax.vlines(t, [0],cuentas, lw=ab, color="#0000ff") #genra rectangulos con altura desde cero hasta el numero de sismos en y, y en x lo posiciona en la fecha correspondiente con un ancho de 50
    ax.set_ylabel('No. de eventos', color='k', fontsize=20)
    
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
        
    labels = ax.xaxis.get_majorticklabels()
    if delta_time > 4500: #si hay mas de 12 años en la grafica pone los tiks inclinados
        plt.setp(labels, rotation=60, fontsize=18)
    else:
        plt.setp(labels, rotation=0, fontsize=18)
    
    ax.yaxis.grid(True, which='major')
        
    #se grafica un texto sobre cada barra del histograma que es el numero de sismos que representa la barra
    a=0;
    for i in cuentas:

        if i>0 and delta_time<20000: #no imprime si el numero de sismos es cero para cierto año #no imprime si el numero de sismos es cero para cierto año
            tx = str(i)
            ax.text(t[a], cuentas[a]+(max(cuentas)*0.01), tx, horizontalalignment='center',fontsize=tik) #sobre 50 para que imprima el texto muy cercano al tope de la barra
        a=a+1    
    

def magnitud_time(time,ax,plt,utc_time_i,utc_time_f,delta_time,mags,mags_grafica):
        
    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)

    ax.scatter(time, mags, marker='o', s=mags_grafica, c="#0000ff", edgecolor='k',alpha=0.5)
    labels = ax.xaxis.get_minorticklabels()
    
    ax.set_ylabel("Magnitud", color='k', fontsize=20)
        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=18)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=20)
    plt.title("Magnitud vs Tiempo", color='k', fontsize=14)
    plt.savefig('Magnitud_Tiempo.png', bbox_inches='tight')
    plt.close()

def error_time(time,ax,plt,utc_time_i,utc_time_f,delta_time,elats,elons,eprofs,rmss):
        
    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)

    ax.plot(time, elats, c="#ff0000",alpha=0.4,label='E_lat')
    ax.plot(time, elons, c="#0000ff",alpha=0.4,label='E_lon')
    ax.plot(time, eprofs, c="#32cd32",alpha=0.4,label='E_prof')
    ax.scatter(time, elats, marker='o', s=25, c="#ff0000", edgecolor='k',alpha=0.6)
    ax.scatter(time, elons, marker='o', s=25, c="#0000ff", edgecolor='k',alpha=0.6)
    ax.scatter(time, eprofs, marker='o', s=25, c="#32cd32", edgecolor='k',alpha=0.6)
    labels = ax.xaxis.get_minorticklabels()
    
    ax.set_ylabel("Error (km)", color='k', fontsize=20)
    ax.legend(loc='upper right',fontsize=25)
        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=18)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=20)
    plt.title("Error vs Tiempo", color='k', fontsize=14)
    plt.savefig('Error_Tiempo.png', bbox_inches='tight')
    plt.close()

def prof_time(time,ax,plt,utc_time_i,utc_time_f,delta_time,profs,mags_grafica):
        
    profs_n=[]
    for p in profs:
        if p < 0:
            profs_n.append(0)
        else:
            profs_n.append(-p)

    profs=profs_n

    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)

    if len(profs) < 100:
        tam_p=200
        alpha=0.8
    if 100<= len(profs) < 200:
        tam_p=150
        alpha=0.8
    if 200<= len(profs) < 300:
        tam_p=100
        alpha=0.7
    if 300<= len(profs) < 500:
        tam_p=80
        alpha=0.7
    if 500<= len(profs) < 1000:
        tam_p=60
        alpha=0.6
    if 1000<= len(profs) < 10000:
        tam_p=40
        alpha=0.6
    if len(profs) >= 10000:
        tam_p=30
        alpha=0.5

    s=mags_grafica #elegir entre tam_p y mags_grafica

    ax.scatter(time, profs, marker='o', s=s, c="#0000ff", edgecolor='k',alpha=alpha)
    #ax.scatter(time, profs, marker='o', s=mags_grafica, c="#0000ff", edgecolor='k',alpha=0.5)
    labels = ax.xaxis.get_minorticklabels()
    
    ax.set_ylabel("Profundidad (km)", color='k', fontsize=20)
        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=18)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=20)
    plt.title("Profundidad vs Tiempo", color='k', fontsize=14)
    plt.savefig('Profundidad_Tiempo.png', bbox_inches='tight')
    plt.close()


def energia_sismica_acumulada(r,ax,delta_time,utc_time_i,utc_time_f):
    #graficando energia sismica acumulada
    time=[]
    mags=[]
    E_sismica=[]

    Es=0
    E_sis_acum=[]
    for i in r:
        fecha=str(i[0])
        horas=str(i[1])
        vec_fecha=fecha.split("-")
        vec_hora=horas.split(":")
        an=int(vec_fecha[0])
        mes=int(vec_fecha[1])
        dia=int(vec_fecha[2])
        t=datetime(an, mes, dia)
        time.append(t)
        mags.append(float(i[5]))
        E_sismica.append(10**float(i[5]))
        
        
    dic={"tiempo": time, "energia_sismica": E_sismica}
    df = pd.DataFrame(dic)
    p="d" #PARA HACER EL GROOUBY del conteo de sismos ELEGIR "d" dias, "w" semanas, o "m" meses
    conteo = df.groupby([df['tiempo'].dt.to_period(p)]).count() #obtiene un dataframe con el conteo dependiendo del periodo p: h:horas, d: dias, m: meses, y: años
    tiempo_sismos=conteo.index.to_timestamp() #obtiene un vector de tiempos en formato datatime
    cuentas=(conteo['tiempo'].tolist()) #obtiene un vector con el conteo

    group_energia_sismica = df.groupby([df['tiempo'].dt.to_period(p)])['energia_sismica'].sum().values

    energia_sismica_acumulada=[]
    e_acum=0
    for es in group_energia_sismica:
        e_acum=e_acum+es
        energia_sismica_acumulada.append(e_acum)

    fig = plt.figure(figsize=(21., 11.))
    ax = fig.add_subplot(111)
    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')
    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)

    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=18)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=20)

    plt.title("Energia Sismica Acumulada",fontsize=20)
    
    ax.plot(tiempo_sismos,energia_sismica_acumulada, color="#0000ff", label='E_sis_acum (Sum ~10^M)')

    ax.set_ylim([0,1.2*max(energia_sismica_acumulada)])
    ax.set_ylabel("Energia_sismica_acumulada (Sum ~10^M) ("+p+")", color='k', fontsize=20)
    ax.set_xlabel("tiempo", color='k', fontsize=20)
    
    ax.legend(loc='upper left',fontsize='large')
    plt.savefig('Energia_sismica_acumulada.png', bbox_inches='tight')
    plt.close()
    '''
    arc=open("Energia_sismica.csv","w")
    for i in range(0,len(tiempo_sismos)):
        arc.writelines(str(tiempo_sismos[i])+","+str(energia_sismica_acumulada[i])+"\n")
    arc.close()
    '''
def del_time(r,delta_time,utc_time_i,utc_time_f,op):

    FAM=1.1
    exp_mag=5.0
    del_v=20 #En leyenda de magnitud entre mayor sea este numero mas se acercan las pepas
    #posicion de leyenda 0  (izquierda, abajo) a 1 (derecha, arriba)
    pos_ley_mag_x=0.1
    pos_ley_mag_y=0.8
    ancho_rec=0.12 #ancho rectangulo de la leyenda

    
    c=[]
    for s in r:                
        c.append( [ UTCDateTime ( s[0]+" "+s[1]) ] + s[2:] ) #c es el mismo catalogo pero fecha y hora en la misma celda

    #b es la version odenada por fecha 
    b = sorted(c, key= itemgetter(0), reverse=False)
    num_sismos=np.arange(0,len(b),1)
    
    time=[]
    time_str=[]
    mags_grafica=[]
    mags=[]
    
    num_sis=[]
    list_mag_t_num=[]
    i=1
    b_n=[]
    for s in b:
        time.append(s[0].datetime)
        time_str.append(str(s[0].datetime))
        num_sis.append(i)
        mag=float(s[4])
                
        mags_grafica.append((((mag/8)*(8-FAM))+FAM)**exp_mag)
        mags.append(mag)

        list_mag_t_num.append([mag,s[0].datetime,i])        
        b_n.append(s+[i])
        i=i+1

    #---- saca los 15 de mayor magnitud y los reordena nuevamente en time
    num_max=5 #numero de sismos mas grandes a resaltar
    if len(b)>num_max:
        

        b_sm=sorted(b_n, key= itemgetter(4), reverse=False)
        b_max_m=sorted(b_sm[-num_max:], key= itemgetter(0), reverse=False)

        time_tc=[]
        time_str_tc=[]
        num_tc=[]
        mags_grafica_tc=[]
        mags_tc=[]
        
        
        for s in b_max_m:
            

            time_tc.append(s[0].datetime)
            time_str_tc.append(str(s[0].datetime))
            num_tc.append(s[-1])
            mag=float(s[4])
                    
            mags_grafica_tc.append((((mag/8)*(8-FAM))+FAM)**exp_mag)
            mags_tc.append(mag)

                       


    def transparencia(b):
        
        if len(b) < 100:
            alpha=0.7            
        if 100<= len(b) < 500:            
            alpha=0.6        
        if 500<= len(b) < 1000:            
            alpha=0.5            
        if 1000<= len(b) < 10000:
            alpha=0.4            
        if len(b) >= 10000:            
            alpha=0.3            
        
        return alpha

    alpha=transparencia(b)
    
    fig = plt.figure(figsize=(21., 11.))
    ax = fig.add_subplot(111)

    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)#
    
    ax.plot(time,num_sis, c="#0000ff",linewidth=0.8)    
    ax.scatter(time,num_sis, marker='o', s=mags_grafica, c="r", edgecolor='k',alpha=alpha)
    ax.scatter(time, num_sis, marker='o', s=10, c="k", edgecolor='r',alpha=0.8) 
    

    if len(b)>num_max:
        
        #ax.scatter(time_tc,time_str_tc, marker='*', s=mags_grafica_tc, c="#ffff00", edgecolor='k',alpha=alpha)
        ax.scatter(time_tc,num_tc, marker='*', s=mags_grafica_tc, c="#ffff00", edgecolor='k',alpha=alpha)
            
    
    #ax2 = ax.twiny()
    
    #-----leyenda mag
    mag_min=int(min(mags))
    mag_max=int(max(mags))+2 
    
    xmin, xmax=ini_graf,fin_graf    
    ymin, ymax = plt.ylim()
    del_x=(xmax-xmin).total_seconds()

    ax.text(xmin + timedelta(seconds=pos_ley_mag_x*del_x),pos_ley_mag_y*(ymax-ymin),"MAGNITUD",fontsize=16,horizontalalignment='center',verticalalignment='center')
    
    for i in range(mag_min,mag_max):
        mag_leyenda=(((i/8)*(8-FAM))+FAM)**exp_mag
        ax.scatter(xmin + timedelta(seconds=pos_ley_mag_x*del_x),( (pos_ley_mag_y-0.05)-(i/del_v) )*(ymax-ymin), marker='o', s=mag_leyenda, c="r", edgecolor='k',alpha=0.4)
        ax.text(xmin + timedelta(seconds=pos_ley_mag_x*del_x),((pos_ley_mag_y-0.05)-(i/del_v))*(ymax-ymin),str(i),fontsize=14,horizontalalignment='center',verticalalignment='center')
            
    # Create a Rectangle patch        
    ini_x=xmin + timedelta(seconds=(pos_ley_mag_x- (ancho_rec/2) )*del_x)    
    ini_y=(pos_ley_mag_y-( (i+3)/del_v))*(ymax-ymin)
    ancho=timedelta(seconds=ancho_rec*del_x)    
    alto=( (pos_ley_mag_y+0.1) *(ymax-ymin)) - ini_y
    rect = patches.Rectangle((ini_x, ini_y),ancho , alto, linewidth=1, edgecolor='k', facecolor='none')
    # Add the patch to the Axes
    ax.add_patch(rect)
    
    if len(b)>num_max:

        #-----texto 5 sismos mayores 
        #ax.text(xmin + timedelta(seconds=pos_ley_mag_x*del_x),pos_ley_mag_y*(ymax-ymin),"MAGNITUD",fontsize=16,horizontalalignment='center',verticalalignment='center')
        pos_ley_txt_x=0.67
        pos_ley_txt_y=0.48
        del_y_txt=0.033

        print ("----------Los "+str(num_max)+" mas grandes------")

        j=0
        for i in b_max_m:
            
            dia_sis=i[0].datetime.strftime("%Y-%m-%d")
            #dia_sis=i[1].datetime.strftime("%Y-%m-%d")
            mag_txt=str(i[4])
            #mag_txt=str(i[0])
            if op == "ap":
                epicentro=str(i[5])
                if len(epicentro)<3: #es posible que al venir de un archivo propio se reporte NA o alguna sigla en el epicentro. Los archivos porpios escritos por estamisma rutina ponen el epi despues de la mag
                    epicentro=""

            if op == "bd":
                epicentro=str(i[-2]) #cuando la opcion es "bd", la lista r escribe al final el epicentro

            txt_s=dia_sis+", M "+mag_txt+", "+epicentro
            ax.text(xmin + timedelta(seconds=pos_ley_txt_x*del_x), pos_ley_txt_y- ( (j*del_y_txt)*(ymax-ymin) ), txt_s,fontsize=14,horizontalalignment='left',verticalalignment='center',bbox=dict(facecolor = 'b', alpha=0.1, pad = 6))
            j=j+1    
            print (txt_s)

        ax.scatter(xmin + timedelta(seconds=(pos_ley_txt_x-0.02)*del_x), pos_ley_txt_y- ( ((j/2)*del_y_txt)*(ymax-ymin) ), marker='*', s=mags_grafica_tc[-1], c="#ffff00", edgecolor='k',alpha=alpha,zorder=100)
        

    #------
    
    ax.set_ylabel("Número de Sismos", color='k', fontsize=20)
    ax.yaxis.set_tick_params(labelsize=18)

    '''
    num_tiks_eje_y =30   
    if len(b)>num_tiks_eje_y:    
        base_tick=int( round ( (len(b)/num_tiks_eje_y) , 2) ) #se muestran al rededor de 40 ticks en y
    else:
        base_tick=1
    ax.yaxis.set_major_locator(tck.MultipleLocator(base=base_tick))    
    '''
    
    labels = ax.xaxis.get_minorticklabels()
    ax.xaxis.set_tick_params(labelsize=16)    
    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=18)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=20)
    
    plt.title("Número de sismos vs Tiempo\n", color='k', fontsize=20)
    plt.savefig('Tiempo_Numero.png', bbox_inches='tight')
    plt.close()

def distancia_utm(lon1,lat1,lon2,lat2):
    #-------- halla distancia epicentral entre 2 eventos en km utm
    coords_1 = (lat1,lon1)
    coords_2 = (lat2,lon2)
    dist_epi=round ( geopy.distance.geodesic(coords_1, coords_2).km , 3)

    return dist_epi

def calcula_distacias_angulos_desde_centro(r,lat_centro,lon_centro):

    
    #----primer cuadrante (cuadrantes manecillas reloj)

    FAM=1.1
    exp_mag=4.7

    c=[]

    for s in r:                
        c.append( [ UTCDateTime ( s[0]+" "+s[1]) ] + s[2:] ) #c es el mismo catalogo pero fecha y hora en la misma celda

    #b es la version odenada por fecha 
    b = sorted(c, key= itemgetter(0), reverse=False)
    num_sismos=np.arange(0,len(b),1)
    
    time=[]
    time_str=[]
    mags_grafica=[]
    mags=[]        
    distancias_epi=[]
    tethas=[]

    for s in b:
        lat=float(s[1])
        lon=float(s[2])
        mag=float(s[4])

        time.append(s[0].datetime)               
        mags_grafica.append((((mag/8)*(8-FAM))+FAM)**exp_mag)
        mags.append(mag)

        
        dist_epi=distancia_utm(lon,lat,lon_centro,lat_centro)
        
        if lat_centro != lat:

            dis_sis_x=distancia_utm(lon_centro,lat_centro,lon,lat_centro)
            dis_sis_y=distancia_utm(lon_centro,lat_centro,lon_centro,lat)
            alfa=math.atan(dis_sis_x/dis_sis_y)*(180/3.141592) #queda en grados
            #cuadrante 1
            if lat >= lat_centro and lon >= lon_centro:            
                tetha=alfa

            #cuadrante 2
            if lat <= lat_centro and lon >= lon_centro:
                tetha=180-alfa

            #cuadrante 3
            if lat <= lat_centro and lon <= lon_centro:
                tetha=180+alfa

            #cuadrante 4
            if lat >= lat_centro and lon <= lon_centro:
                tetha=360-alfa 
        else:
            if  lon <= lon_centro:
                tetha=270
            if  lon > lon_centro:
                tetha=90


        distancias_epi.append(dist_epi)
        tethas.append(tetha)
    

    return time,distancias_epi,mags_grafica,tethas,mags

def distancia_al_centro(r,time,ax,plt,utc_time_i,utc_time_f,delta_time,lat_centro,lon_centro,mags_grafica):

    time,distancias_epi,mags_grafica,tethas,mags=calcula_distacias_angulos_desde_centro(r,lat_centro,lon_centro)

    #--- grafica distancia al centro vrs time
    ini_graf,fin_graf=genera_ejes(ax,delta_time,utc_time_i,utc_time_f)

    if len(time) < 100:        
        alpha=0.7
    if 100<= len(time) < 200:        
        alpha=0.6
    if 200<= len(time) < 300:        
        alpha=0.6
    if 300<= len(time) < 500:        
        alpha=0.5
    if 500<= len(time) < 1000:        
        alpha=0.5
    if 1000<= len(time) < 10000:        
        alpha=0.4
    if len(time) >= 10000:        
        alpha=0.3

    ax.scatter(time, distancias_epi, marker='o', s=mags_grafica, c="#0000ff", edgecolor='k',alpha=alpha)
    ax.plot(time, distancias_epi, c="#0000ff",linewidth=0.8,alpha=0.1)
    labels = ax.xaxis.get_minorticklabels()
    
    ax.set_ylabel("Distancia Epicentral (km)\n(xo,yo): ("+str(lon_centro)+"° , "+str(lat_centro)+"°)", color='k', fontsize=20)
        
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=18)
    labels = ax.xaxis.get_minorticklabels()
    plt.setp(labels, rotation=60, fontsize=18)
    labels = ax.get_xticklabels() 
    plt.setp(labels, rotation=0, fontsize=20)
    plt.title("Distancia Epicentral (km) vs Tiempo", color='k', fontsize=14)
    plt.savefig('Distancia_Epi_Tiempo.png', bbox_inches='tight')
    plt.close()

    #--- grafica windrose

    fig = plt.figure(figsize=(21., 11.))
    ax = fig.add_subplot(111)
    ax = WindroseAxes.from_ax()
    plt.title("Número de eventos por Magnitud")
    variable=mags #distancias_epi #mags
    bins=np.arange( int(min (variable)), int(max(variable))+1.5, 0.5 )
    #ax.contourf(tethas, variable, bins=bins, normed=False, cmap=cm.hot, alpha=0.2)
    #ax.contour(tethas, variable, bins=bins, normed=False, colors='black',linewidth=0.3)
    ax.bar(tethas, variable, normed=False, opening=0.8, edgecolor='white',bins= bins,nsector=36,cmap=cm.jet,alpha=1)
    #ax.box(tethas, variable, normed=False, edgecolor='white',bins= bins,nsector=36,cmap=cm.jet,alpha=1)
    #ax.set_legend("")
    ax.legend(loc='upper right',fontsize=8,title="Magnitud")
    plt.savefig('diagrama_rosa.png', bbox_inches='tight')
    plt.close()


def histo_sismos(r,time_i,time_f,op,lat_centro,lon_centro,tipo_busqueda):

    if len(r) > 0:

        #-------ordena r adecuadamente por fecha
        #-convierte fecha en utcdatetime (uniendo fecha y hora)
        c=[]
        for s in r:             
            c.append( [ UTCDateTime ( s[0]+" "+s[1]) ] + s[2:] ) #c es el mismo catalogo pero fecha y hora en la misma celda

        #b es la version odenada por fecha 
        b = sorted(c, key= itemgetter(0), reverse=False)

        #vuelve a separar fecha y hora en casillas diferentes,
        r=[]
        for rf in b:
            
            fecha=str(rf[0].datetime)
                    
            vec_fecha=fecha.split(" ")
            an_me_di=vec_fecha[0]
            hora_utc=vec_fecha[1]           

            r.append([an_me_di,hora_utc]+rf[1:]) #arma una lista concatenando fecha, hora, (lat, lon (rf[1:3])), prof, y el resto (rf[4:])

        #

        print("Generando Histogramas...")

        time=[]
        time_destacados=[]    
        mags=[]
        mags_grafica=[]
        profs=[]
        lats=[]
        lons=[]
        elats=[]
        elons=[]
        eprofs=[]
        rmss=[]
        gaps=[]
        
        r_destacados=[]
        for i in r:

            fecha=str(i[0])
            horas=str(i[1])
            
            if "/" in fecha:
                fecha=feha.replace("/","-")
            if "/" in horas:
                horas=horas.replace("/",":")
            if "-" in horas:
                horas=horas.replace("-",":")
            
            utc_fecha=UTCDateTime(fecha+horas).datetime #-timedelta(hours=5) #para hora local restar 5 horas
            time.append(utc_fecha)
            
            mag=float(i[5])
            
            if mag>=4.0: #magnitud para destacados
                time_destacados.append(utc_fecha)
                r_destacados.append(i)

            mags.append(mag)
            FAM=1.8
            exp_mag=4
            mags_grafica.append((((mag/8)*(8-FAM))+FAM)**exp_mag)    
            profs.append(float(i[4]))            
            lats.append(float(i[2]))
            lons.append(float(i[3]))

            if op == "bd":

                if float(i[6]) < 10: #retira valores anomalos de la grafica
                    rmss.append(float(i[6]))
                if float(i[7]) < 360: #retira valores anomalos de la grafica
                    gaps.append(float(i[7]))

                elats.append(float(i[8]))
                elons.append(float(i[9]))
                eprofs.append(float(i[10]))

                
            
            if op == "ap" and len(i)>10: #en el elemnto 6 del archivo propio va el epicentro, en el elemento 11 esta eprof 

                #print (i[7],i[8],i[9],i[10],i[11],i[12])


                if float(i[7]) < 10: #retira valores anomalos de la grafica
                    rmss.append(float(i[7]))
                else:
                    rmss.append(0)

                if float(i[8]) < 360: #retira valores anomalos de la grafica
                    gaps.append(float(i[8]))
                else:
                    gaps.append(0)
                
                elats.append(float(i[9]))
                
                elons.append(float(i[10]))
                
                eprofs.append(float(i[11]))

        print("-------------Datos de la sismicidad------")
        print("El Reporte contiene "+str(len(r))+" eventos" )
        print("La sismicidad se localiza entre ("+str(min(lons))+","+str(max(lons))+") y ("+str(min(lats))+","+str(max(lats))+")" )
        print("Sismo mas antiguo: "+min(time).strftime("%Y-%m-%d")+", Sismo mas reciente: "+max(time).strftime("%Y-%m-%d"))
        print("Magnitud menor: "+str(min(mags))+", Magnitud mayor: "+str(max(mags)))
        print("Menor profundidad: "+str(min(profs))+", Mayor profundidad: "+str(max(profs)))       


        #----genera dataframe para generar conteos de sismidad en la unidad de tiempo
        dic={"tiempo": time}
        df_tiempo = pd.DataFrame(dic)

        dic_destacados={"tiempo": time_destacados}
        df_tiempo_destacados = pd.DataFrame(dic_destacados)
        utc_time_i = UTCDateTime(time_i).datetime
        utc_time_f = UTCDateTime(time_f).datetime
        delta_time=(utc_time_f-utc_time_i).total_seconds()/3600/24 #en dias
        #print(delta_time)
        #(df_tiempo['tiempo'].max()-df_tiempo['tiempo'].min()).total_seconds()/3600/24 #en dias
        #con este tomara desicion de que histogramas hacer

        if 2 < delta_time < 6: #hara histograma horas 
            #hoas
            ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
            histo_conteo_horas(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time,r)
            plt.title('Sismicidad por Hora\nCatalogo del SGC', fontsize=19)
            plt.savefig('Sismicidad_Hora.png', bbox_inches='tight')
            plt.close()
            
            
            if len (time_destacados)>0:
                ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
                histo_conteo_horas(df_tiempo_destacados,ax,plt,utc_time_i,utc_time_f,delta_time,r_destacados)
                plt.title('Sismicidad por Hora Destacada\nCatalogo del SGC', fontsize=19)
                plt.savefig('Sismicidad_Hora_Destacada.png', bbox_inches='tight')
                plt.close()
        
        if delta_time >= 3: #hara histograma dias

            #dias
            ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
            histo_conteo_dias(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time,r)
            plt.title('Sismicidad Diaria\nCatalogo del SGC', fontsize=19)
            plt.savefig('Sismicidad_Diaria.png', bbox_inches='tight')
            plt.close()
                        
            if len (time_destacados)>0:
                ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
                histo_conteo_dias(df_tiempo_destacados,ax,plt,utc_time_i,utc_time_f,delta_time,r_destacados)
                plt.title('Sismicidad Diaria Destacada\nCatalogo del SGC', fontsize=19)
                plt.savefig('Sismicidad_Diaria_Destacada.png', bbox_inches='tight')
                plt.close()
            #fin_dias
        
        if delta_time > 180: #hara histograma meses

            #hace grafica cuando han ocurrido sismos de acuerdo al dia del mes 
            histo_conteo_repeticion_dia(time,delta_time)

            #meses
            ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
            histo_conteo_meses(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time)
            plt.title('Sismicidad mensual\nCatalogo del SGC', fontsize=19)
            plt.savefig('Sismicidad_Mensual.png', bbox_inches='tight')
            plt.close()
                        
            if len (time_destacados)>0:
                ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
                histo_conteo_meses(df_tiempo_destacados,ax,plt,utc_time_i,utc_time_f,delta_time)
                plt.title('Sismicidad Mensual Destacados\nCatalogo del SGC', fontsize=19)
                plt.savefig('Sismicidad_Mensual_Destacada.png', bbox_inches='tight')
                plt.close()
            #fin_meses
        
        if delta_time > 1080: #hara histograma años (lo hace si hay un poco mas de tres años de datos)
            
            #hace grafica cuando han ocurrido sismos de acuerdo al mes del año 
            histo_conteo_repeticion_mes(time,delta_time) 

            #años
            ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
            histo_conteo_ans(df_tiempo,ax,plt,utc_time_i,utc_time_f,delta_time)
            plt.title('Sismicidad anual\nCatalogo del SGC', fontsize=19)
            plt.savefig('Sismicidad_Anual.png', bbox_inches='tight')
            plt.close()    
            
            if len (time_destacados)>0:
                ax,plt,utc_time_i,utc_time_f = genera_cuerpo_grafica(delta_time,time_i,time_f)
                histo_conteo_ans(df_tiempo_destacados,ax,plt,utc_time_i,utc_time_f,delta_time)
                plt.title('Sismicidad Anual Destacada\nCatalogo del SGC', fontsize=19)
                plt.savefig('Sismicidad_anual_destacada.png', bbox_inches='tight')
                plt.close()

        #-------------------------------------------------------------------------------
        
        # ---------------------GRAFICAs en TIME (magnitud, profundidad y energia)

        ax,plt,utc_time_i,utc_time_f=genera_cuerpo_grafica(delta_time,time_i,time_f)        
        magnitud_time(time,ax,plt,utc_time_i,utc_time_f,delta_time,mags,mags_grafica)

        ax,plt,utc_time_i,utc_time_f=genera_cuerpo_grafica(delta_time,time_i,time_f)
        prof_time(time,ax,plt,utc_time_i,utc_time_f,delta_time,profs,mags_grafica)



        if tipo_busqueda=="c":
            ax,plt,utc_time_i,utc_time_f=genera_cuerpo_grafica(delta_time,time_i,time_f)
            distancia_al_centro(r,time,ax,plt,utc_time_i,utc_time_f,delta_time,lat_centro,lon_centro,mags_grafica)

        del_time(r,delta_time,utc_time_i,utc_time_f,op)

        #grafica energia acumulada
        energia_sismica_acumulada(r,ax,delta_time,utc_time_i,utc_time_f)

    

        #------------------graficando numero de sismos por profundidad------------------

        fig = plt.figure(figsize=(21., 11.))
        ax1 = fig.add_subplot(111)
        ax1.set_xlim( min(profs),1.1*max(profs))
        
        if int(max(profs)/20) != 0:
            bins=np.arange(0,max(profs),int(max(profs)/20))
        if int(max(profs)/20) == 0:
            bins=np.arange(0,max(profs),1)
        #N, bins, patches = ax1.hist(profs, bins=[0,30,60,90,120,150,180], alpha=0.9, color="#0000ff",edgecolor = "#e5e5ff", linewidth=3)
        N, bins, patches = ax1.hist(profs, bins, alpha=0.9, color="#0000ff",edgecolor = "#e5e5ff", linewidth=3)

        ax1.set_ylabel('No. de eventos', color='k', fontsize=20)
        ax1.set_xlabel(r'Profundidad (km)', color='k', fontsize=20)

        ax1.xaxis.grid(True, which='major')
        ax1.yaxis.grid(True, which='major')

        ax1.xaxis.set_tick_params(labelsize=20)
        ax1.yaxis.set_tick_params(labelsize=20)

        plt.title('Profundidad vs Nº de Sismos\nCatalogo del SGC', fontsize=19)
        plt.savefig('Profundidad vs Numero de Sismos.png', bbox_inches='tight')
        plt.close()

            
        #------------------graficando numero de sismos por magnitud------------------

        fig = plt.figure(figsize=(21., 11.))
        ax1 = fig.add_subplot(111)
        ax1.set_xlim( 0, 1.1*max(mags))
        
        #bins=np.arange(0,max(mags),round((max(mags)/20),1))
        bins=np.arange(0,max(mags) +1,0.5)
        ax1.hist(mags, bins=bins, alpha=0.9, color="#0000ff",edgecolor = "#e5e5ff", linewidth=3)

        ax1.set_ylabel('No. de eventos', color='k', fontsize=20)
        ax1.set_xlabel(r'Magnitud', color='k', fontsize=20)

        ax1.xaxis.set_tick_params(labelsize=20)
        ax1.yaxis.set_tick_params(labelsize=20)

        plt.title('Magnitud vs Nº de Sismos\nCatalogo del SGC', fontsize=19)
        plt.savefig('Magnitud_Numero_Sismos.png', bbox_inches='tight')
        plt.close()        
        

        #--------------grafica errores, rms y gap-------------

        #--------------grafica errores
        if op == "bd" or  (op == "ap" and len(i)>10):

            #---------------error en time
            ax,plt,utc_time_i,utc_time_f=genera_cuerpo_grafica(delta_time,time_i,time_f)
            error_time(time,ax,plt,utc_time_i,utc_time_f,delta_time,elats,elons,eprofs,rmss)

            # --------------errores en frecuencia
            fig = plt.figure(figsize=(21., 11.))

            def labels(eje,ax):
                ax.set_ylabel('No. de eventos', color='k', fontsize=20)
                ax.set_xlabel('Error_'+eje, color='k', fontsize=20)
                ax.xaxis.grid(True, which='major')
                ax.yaxis.grid(True, which='major')

                ax.xaxis.set_tick_params(labelsize=20)
                ax.yaxis.set_tick_params(labelsize=20)

                ax.legend(fontsize=20)
            
            max_err=max (max(elats),max(elons),max(eprofs))
            bins=np.arange(0,max_err +1,0.5) 
            ax1 = fig.add_subplot(311)
            ax1.set_xlim( 0, 1.1*max_err)            
            ax1.hist(elats, bins=bins, alpha=0.9, color="b",edgecolor = "b", linewidth=1,label='e_lat')
            ax1.set_title('Error vs Nº de Sismos\nCatalogo del SGC', fontsize=19)  
            labels("lat",ax1)

            ax2 = fig.add_subplot(312)
            ax2.set_xlim( 0, 1.1*max_err)             
            ax2.hist(elons, bins=bins, alpha=0.9, color="r",edgecolor = "r", linewidth=1,label='e_lon')  
            labels("lon",ax2)

            ax3 = fig.add_subplot(313)
            ax3.set_xlim( 0, 1.1*max_err)             
            ax3.hist(eprofs, bins=bins, alpha=0.9, color="g",edgecolor = "g", linewidth=1,label='e_prof')  
            labels("prof",ax3)
            
            plt.subplots_adjust(hspace=0.45)

            
            plt.savefig('Error_Numero_Sismos.png', bbox_inches='tight')
            plt.close()

            #--------------grafica rms-------------

            fig = plt.figure(figsize=(21., 11.))
            ax1 = fig.add_subplot(111)
            ax1.set_xlim( 0, 1.1*max(rmss))
            
            #bins=np.arange(0,max(mags),round((max(mags)/20),1))
            #print(max(rmss))
            bins=np.arange(0,max(rmss) +0.1,0.2)
            ax1.hist(rmss, bins=bins, alpha=0.9, color="#0000ff",edgecolor = "#e5e5ff", linewidth=3)
            
            ax1.set_ylabel('No. de eventos', color='k', fontsize=20)
            ax1.set_xlabel('RMS', color='k', fontsize=20)
            ax1.xaxis.grid(True, which='major')
            ax1.yaxis.grid(True, which='major')

            ax1.xaxis.set_tick_params(labelsize=20)
            ax1.yaxis.set_tick_params(labelsize=20)

            plt.title('RMS vs Nº de Sismos\nCatalogo del SGC', fontsize=19)
            plt.savefig('RMS_Numero_Sismos.png', bbox_inches='tight')
            plt.close()

            #--------------grafica GAP-------------

            fig = plt.figure(figsize=(21., 11.))
            ax1 = fig.add_subplot(111)
            ax1.set_xlim( 0, 1.1*max(gaps))
            
            #bins=np.arange(0,max(mags),round((max(mags)/20),1))
            bins=np.arange(0,max(gaps) +10,10)
            ax1.hist(gaps, bins=bins, alpha=0.9, color="#0000ff",edgecolor = "#e5e5ff", linewidth=3)
            
            ax1.set_ylabel('No. de eventos', color='k', fontsize=20)
            ax1.set_xlabel('GAP', color='k', fontsize=20)
            ax1.xaxis.grid(True, which='major')
            ax1.yaxis.grid(True, which='major')

            ax1.xaxis.set_tick_params(labelsize=20)
            ax1.yaxis.set_tick_params(labelsize=20)

            plt.title('GAP vs Nº de Sismos\nCatalogo del SGC', fontsize=19)
            plt.savefig('GAP_Numero_Sismos.png', bbox_inches='tight')
            plt.close()
        
        

        
        #---------------------------------------------------
            
    
        
        




    else:
        print("El reporte de sismicidad esta vacio, no se pueden generar histogramas")



