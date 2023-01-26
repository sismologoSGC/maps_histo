import simplekml


def r2kml (r):

	print("Generando mapa en Google Earth 'Mapa_Google_Earth.kml ...")
	kml=simplekml.Kml()

	style1 = simplekml.Style() #creates shared style for all points
	style1.iconstyle.color = 'ff0000ff' #magenta
	style1.iconstyle.icon.href ='http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png' #can change to any desired icon URL

	style2 = simplekml.Style() #creates shared style for all points
	style2.iconstyle.color = 'ff0045ff' #magenta
	style2.iconstyle.icon.href ='http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png' #can change to any desired icon URL

	style3 = simplekml.Style() #creates shared style for all points
	style3.iconstyle.color = 'ff00ffff' #magenta
	style3.iconstyle.icon.href ='http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png' #can change to any desired icon URL

	style4 = simplekml.Style() #creates shared style for all points
	style4.iconstyle.color = 'ff008000' #magenta
	style4.iconstyle.icon.href ='http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png' #can change to any desired icon URL

	style5 = simplekml.Style() #creates shared style for all points
	style5.iconstyle.color = 'ffff0000' #magenta
	style5.iconstyle.icon.href ='http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png' #can change to any desired icon URL

	for i in r:

		lat=float(i[2])
		lon=float(i[3])
		prof=float(i[4])
			
			#en i[5] esta ml y en i[6] esta mw (en caso que este), el codigo busca si hay valor de mw, si lo hay, toma este como la magnitud mag, de lo contrario toma a ml=mag
		mag=float(i[5])

		descripcion=f"CATALOGO SGC\nFecha: {i[0]} - {i[1]} (UTC) \nM = {mag} - Profundidad = {prof} \nRMS = {i[7]} - GAP = {i[8]}"
	

		if prof<=30:
			
			pnt=kml.newpoint(coords=[(lon,lat)],description=descripcion)
			pnt.style = style1 # magent
		if 30 < prof <=60:
	  		
	  		pnt=kml.newpoint(coords=[(lon,lat)],description=descripcion)
	  		pnt.style = style2
		if 60 < prof <=120:
	  		
	  		pnt=kml.newpoint(coords=[(lon,lat)],description=descripcion)
	  		pnt.style = style3
		if 120 < prof <=300:
	  		
	  		pnt=kml.newpoint(coords=[(lon,lat)],description=descripcion)
	  		pnt.style = style4
		if prof >300:
	  		
	  		pnt=kml.newpoint(coords=[(lon,lat)],description=descripcion)
	  		pnt.style = style5
		pnt.style.iconstyle.scale = mag*0.2
		


	kml.save('Mapa_Google_Earth.kml')
	