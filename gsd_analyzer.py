# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GSDAnalyzer
                                 A QGIS plugin
 This plugin is to compute and spatially visualize the local GSD for a mapping mission at constant altitude with respect to take off.
                              -------------------
        begin                : 2019-01-17
        git sha              : $Format:%H$
        copyright            : (C) 2019 by Jonas Hener - Avy BV
        email                : jonas@avy.eu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt4.QtGui import QAction, QIcon, QFileDialog
from qgis.core import *
from qgis.core import QgsMapLayerRegistry
from osgeo import gdal, ogr
from qgis.gui import QgsMessageBar
import struct
import sys
import processing
import numpy as np
import math
from osgeo.gdalconst import *
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from gsd_analyzer_dialog import GSDAnalyzerDialog
import os.path



class GSDAnalyzer:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgisInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'GSDAnalyzer_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)
    
        # Create the dialog (after translation) and keep reference
        self.dlg = GSDAnalyzerDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&GSD Analyzer')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'GSDAnalyzer')
        self.toolbar.setObjectName(u'GSDAnalyzer')

        #Clear Output File text and connect PushButton
        self.dlg.OutLine.clear()
        self.dlg.OutButton.clicked.connect(self.select_output_file)

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('GSDAnalyzer', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/GSDAnalyzer/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Analyze the GSD of your mapping mission'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&GSD Analyzer'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar
            
            
    def select_output_file(self):
        
        filename = QFileDialog.getSaveFileName(self.dlg, "Select output file ","", '*.tif')
        self.dlg.OutLine.setText(filename)


    def getPointGeo(self, myLayers):
        #assign Vector Layer from dlg
        vectorLayerIndex = self.dlg.VectorBox.currentIndex()
        vectorLayer = myLayers[vectorLayerIndex]
        features = vectorLayer.getFeatures()
        for feature in features:
            geom = feature.geometry()
            # show some information about the feature
            if geom.type() == QGis.Point:
                x = geom.asPoint()
                geoX = int(x[0])
                geoY = int(x[1])
                return geoX, geoY
            else:
                print "No Point Geometry Feature"


    def geoToInd(self, geoX, geoY, rasterGT):
        #only works for geotransforms with no rotation.
        indX = int((geoX - rasterGT[0]) / rasterGT[1]) #x pixel
        indY = int((geoY - rasterGT[3]) / rasterGT[5]) #y pixel
    
        return indY, indX #reverse order intentional, don't ask why...


    def analysis(self, DEM, a, p, hfov, vfov, noData = -99): #auxfile as argument for debugging
        
        DEM[DEM<0] = -1
        
        deltaH = a * np.ones(DEM.shape) - DEM
        deltaH[deltaH<0] = 0
        blank = np.zeros(DEM.shape)
        
        h_e = hfov*math.pi/360.     #deliberately by 2
        v_e = vfov*math.pi/360.     #deliberately by 2

        for i in range(DEM.shape[0]):
            for j in range(DEM.shape[1]):

                h = deltaH[i, j]
                
                if h > 1 and DEM[i, j] > 0:
                    hPrint = 2*h*math.tan(h_e) #sensor footprint width in m
                    vPrint = 2*h*math.tan(v_e) #sensor footprint length in m
                    fPrint = vPrint*hPrint # sensor footprint area in m^2
                    GSD = math.sqrt(fPrint * 100 * 100 / p) #m^2 to cm^2 by pixels
                    blank[i, j] = GSD
                    #auxfile.write("A, h = "+ str(h) + " , DEM = " + str(DEM[i,j]) + " , fPrint = " + str(fPrint) + " , GSD = " + str(blank[i, j]) + "\n")
                
                else:
                    blank[i, j] = noData
                    deltaH[i, j] = noData
                    #auxfile.write("B, h = "+ str(h) + " , DEM = " + str(DEM[i,j]) +  "\n")
        

        return blank, deltaH


    def writeOut(self, openRaster, filename, npOut, noData = -99):
        # create the output image1 from frame of open raster DEM and output array
        driver = openRaster.GetDriver()
        rasterGT = openRaster.GetGeoTransform()
        rows = openRaster.RasterYSize
        cols = openRaster.RasterXSize
                
        #print driver
        outRaster = driver.Create(filename, cols, rows, 1, GDT_Int32)
        if outRaster is None:
            print 'Could not create Tif File'
            sys.exit(1)
        
        outBand = outRaster.GetRasterBand(1)
        #outData = np.zeros(npBand1.shape) #for when cropping DEM is allowed
        #outData[loX:hiX, loY:hiY] = out   #for when cropping DEM is allowed
        outData = npOut
        
        # write the data
        outBand.WriteArray(outData, 0, 0)
        
        # flush data to disk, set the NoData value and calculate stats
        outBand.FlushCache()
        outBand.SetNoDataValue(noData)
        
        # georeference the image and set the projection
        outRaster.SetGeoTransform(openRaster.GetGeoTransform())
        outRaster.SetProjection(openRaster.GetProjection())
        
        del outData

    
    def run(self):
        
        #gets the layers loaded in QGIS and adds it to the comboBox object from the plugin dialog
        iface = self.iface
        #clear combos
        self.dlg.RasterBox.clear();self.dlg.VectorBox.clear()

        #take layers from canvas and sort
        myLayers = iface.mapCanvas().layers()
    
        for i in range(len(myLayers)):
            myLayer = myLayers[i]

            #add all to RasterBox
            self.dlg.RasterBox.addItem(myLayer.name(),myLayer.id())
            
            #add all to VectorBox
            self.dlg.VectorBox.addItem(myLayer.name(),myLayer.id())
        
        # TODO: Sort for valid inputs and sync ComboBox index to mapCanvas
        
        #-----------------------------------------------------------------------------#

        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
           
            #for debugging and logging
            filename = self.dlg.OutLine.text()
            parts = filename.split(".")
            auxfilename = parts[0] + "-aux.txt"
            #auxfile = open(auxfilename, 'w')
            
            #for logging metadata, mainly analysis parameters
            metafilename = parts[0] + "-0.txt"
            metafile = open(metafilename, 'w')
            
            #assign Raster Layer from dlg
            rasterLayerIndex = self.dlg.RasterBox.currentIndex()
            rasterLayer = myLayers[rasterLayerIndex]
            
            #-------------------------------------------------------------------------#
            
            #open raster, get dims and band1
            gdal.AllRegister()
            provider = rasterLayer.dataProvider()
            openRaster = gdal.Open(str(provider.dataSourceUri()), gdal.GA_Update)
            
            if openRaster is None:
                print 'Could not open image file'
                sys.exit(1)
            
            rasterGT = openRaster.GetGeoTransform()
            rows = openRaster.RasterYSize
            cols = openRaster.RasterXSize
            band1 = openRaster.GetRasterBand(1)
            npBand1 = band1.ReadAsArray(0,0,cols,rows)
            
            #-------------------------------------------------------------------------#
            
            #get Geo Coordinates of Observer
            geoX, geoY = self.getPointGeo(myLayers)
            #Convert from map to pixel coordinates
            indX, indY = self.geoToInd(geoX, geoY, rasterGT)
            
            #determine Observer Height
            obsZ = npBand1[indX, indY]
            #auxfile.write("obsZ: " + str(obsZ))
            altGCS = int(self.dlg.AltLine.text())
            #auxfile.write("altGCS: " + str(altGCS))
            alt = obsZ+altGCS
            #auxfile.write("alt: " + str(alt))

            #get sensor characteristics:
            pix = float(self.dlg.PixelLine.text())*1000000 #pixelcount, from MP to P
            HFoV = float(self.dlg.HLine.text())
            VFoV = float(self.dlg.VLine.text())
       
            
            #np.savetxt("/Users/jonashener/Desktop/yomama-aux0" + '.csv', npBand1, delimiter=",")

            #-------------------------------------------------------------------------#
            
            out1, out2 = self.analysis(npBand1, alt, pix, HFoV, VFoV)   #auxfile as argument for debugging

            #-------------------------------------------------------------------------#
            
            # create the output image1
            filename1 = parts[0] + "-1.tif"
            self.writeOut(openRaster, filename1, out1)
            
            # create the output image2
            filename2 = parts[0] + "-2.tif"
            self.writeOut(openRaster, filename2, out2)


            metafile.write("input raster:  " + str(self.dlg.RasterBox.currentText()) + "\n" + "input vector:  " + str(self.dlg.VectorBox.currentText()) + "\n" + "flight altitude w.r.t GCS:  " + str(altGCS) + " m" + "\n" + "GCS height:  " + str(obsZ) + " m" + "\n" + "pixelcount:  " + str(pix/1000000.) + " MP" + "\n" + "hfov:  " + str(HFoV) + " deg" + "\n" + "vfov:  " + str(VFoV) + " deg" + "\n")
            
            

