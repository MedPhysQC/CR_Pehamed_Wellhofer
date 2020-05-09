#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This code is an analysis module for WAD-QC 2.0: a server for automated 
# analysis of medical images for quality control.
#
# The WAD-QC Software can be found on 
# https://bitbucket.org/MedPhysNL/wadqc/wiki/Home
# 
#
# Changelog:
#   20200508: dropping support for python2; dropping support for WAD-QC 1; toimage no longer exists in scipy.misc
#   20190426: Fix for matplotlib>3
#   20161220: Removed class variables; removed testing stuff
#   20161216: added use_mustbeinverted param
#   20160802: sync with pywad1.0
#   20160622: removed adding limits (now part of analyzer)
#   20160620: remove quantity and units
#
# ./QCXRay_wadwrapper.py -c Config/cr_philips_umcu_series.json -d TestSet/StudyPehamed -r results_pehamed.json

__version__ = '20200508'
__author__ = 'aschilham'

import os
# this will fail unless wad_qc is already installed
from wad_qc.module import pyWADinput
from wad_qc.modulelibs import wadwrapper_lib

if not 'MPLCONFIGDIR' in os.environ:
    import pkg_resources
    try:
        #only for matplotlib < 3 should we use the tmp work around, but it should be applied before importing matplotlib
        matplotlib_version = [int(v) for v in pkg_resources.get_distribution("matplotlib").version.split('.')]
        if matplotlib_version[0]<3:
            os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
    except:
        os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
        
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

try:
    import pydicom as dicom
except ImportError:
    import dicom
import datetime
import QCXRay_lib

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

def logTag():
    return "[QCXRay_wadwrapper] "

# helper functions
"""
    roomWKZ1 = Room("WKZ1",outvalue=1023,tablesid=1150,wallsid=2000, tablepid=65, wallpid=50,phantom=lit.stWellhofer)
    <params>
      <roomname>WKZ1</roomname>
      <phantom>wellhofer</phantom>
      <tablesidmm>1150</tablesidmm>
      <tablepidmm>65</tablepidmm>
      <wallsidmm>2000</wallsidmm>
      <wallpidmm>50</wallpidmm>
      <outvalue>1023</outvalue>
      
      <sensitivities>
        <threshold date="20100101" value="35" />
      </sensitivities>
      
      <sdthreshold>40</sdthreshold>
    </params>
"""
def _getRoomDefinition(params):
    # Use the params in the config file to construct an Scanner object
    try:
        # a name for identification
        roomname = params['roomname']

        # phantom name (only pehamed or wellhofer)
        phantoms_supported = ['pehamed','wellhofer']
        phantom = params['phantom']
        if not phantom in phantoms_supported:
            raise ValueError(logTag()+' unsupported phantom %s'%phantom)

        # load the locations of markers on the linepair pattern. if these are not given, use the hardcoded values
        linepair_type = params['linepair_type']
        if not linepair_type in ['typ38']:
            raise ValueError('Incorrect linepair type %s'%linepair_type)

        linepairmarkers = {}
        try:
            mnames = ['xymm1.8','xymm0.6','xymm1.4','xymm4.6']
                
            for mname in mnames:
                marker  = params[mname]
                vals = [float(v) for v in  marker.split(';')]
                linepairmarkers[mname] = [vals[0],vals[1]]
        except:
            print(logTag()+' exact locations of markers on linepair pattern not supplied by config. Using empirical values; please check if these are valid here.')
            
        # Source to Detector distance and Patient to Detector distance for wall and table (both in mm)
        tablepidmm  = float(params['tablepidmm'])
        wallpidmm   = float(params['wallpidmm'])

        outvalue    = -1 # not supplied
        wallsidmm   = -1 # not supplied
        try: # only for FCR
            wallsidmm   = float(params['wallsidmm'])
        except:
            pass
        
        tablesidmm  = -1 # not supplied
        try: # only for FCR
            tablesidmm  = float(params['tablesidmm'])
        except:
            pass

        outvalue    = -1 # not supplied
        try: # only for FCR
            # pixelvalue that defines 'outside phantom' use '-1' to calculate from four cornerpoints
            outvalue    = int(params['outvalue'])
        except:
            pass

        # for fcr systems there is no dicom tag to indicate wall or table, but a hack on SD or Sensitivity is possible
        try:
            thresholdlist = []
            sensitivitydatavalue = params['sensitivitydatavalue']
            for sn in sensitivitydatavalue.split('|'):
                vals = sn.split(';')
                thresholdlist.append([int(vals[0]),int(vals[1])])
            return QCXRay_lib.Room(roomname, outvalue=outvalue,
                                   tablesid=tablesidmm, wallsid=wallsidmm, 
                                   tablepid=tablepidmm, wallpid=wallpidmm,
                                   phantom=phantom, sens_threshold = thresholdlist,
                                   linepairmarkers=linepairmarkers)
        except:
            pass

        # no sensitivity threshold, so try if threshOnSD exists
        try:
            sdthreshold = float(params["sdthreshold"])
            return QCXRay_lib.Room(roomname, outvalue=outvalue,
                                   tablesid=tablesidmm, wallsid=wallsidmm, 
                                   tablepid=tablepidmm, wallpid=wallpidmm,
                                   phantom=phantom, sdthresh = sdthreshold,
                                   linepairmarkers=linepairmarkers)
        except:
            pass

        try:
            use_mustbeinverted = params['use_mustbeinverted']
            if use_mustbeinverted.lower() == 'true':
                mustbeinverted = True
            elif use_mustbeinverted.lower() == 'false':
                mustbeinverted = False
            else:
                raise ValueError('Unknown value %s for param use_mustbeinverted'%use_mustbeinverted)
        except:
            mustbeinverted = None

        # no artificial thresholds present or needed
        return QCXRay_lib.Room(roomname, outvalue=outvalue,
                               tablesid=tablesidmm, wallsid=wallsidmm, 
                               tablepid=tablepidmm, wallpid=wallpidmm,
                               phantom=phantom,linepairmarkers=linepairmarkers,
                               mustbeinverted=mustbeinverted)
    except AttributeError as e:
        raise ValueError(logTag()+" missing room definition parameter!"+str(e))


def override_settings(cs, params):
    """
    Look for 'use_' params in to force behaviour of module
    """
    return
    try:
        use_mustbeinverted = params.find('use_mustbeinverted').text
        if use_mustbeinverted.lower() == 'true':
            cs.mustbeinverted = True
        elif use_mustbeinverted.lower() == 'false':
            cs.mustbeinverted = False
        else:
            raise ValueError('Unknown value %s for param use_mustbeinverted'%use_mustbeinverted)
    except:
        pass


###### Series wrappers
def qc_series(data, results, action):
    """
    QCXRay_UMCU checks:
        Horizontal uniformity
        XRayEdges
        LowContrast
        DynamicRange
        MTF

    Workflow:
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
        6. Build artefact picture thumbnail
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    inputfile = data.series_filelist[0]  # give me a filename

    ## 2. Check data format
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(inputfile,headers_only=False,logTag=logTag())

    ## 3. Build and populate qcstructure
    remark = ""
    qclib = QCXRay_lib.XRayQC()
    room = _getRoomDefinition(params)
    cs = QCXRay_lib.XRayStruct(dcmInfile,pixeldataIn,room)
    cs.verbose = False # do not produce detailed logging
    override_settings(cs, params)

    ## 4. Run tests
    error,msg = qclib.QC(cs)

    ## 5. Build xml output
    ## Struct now contains all the results and we can write these to the WAD IQ database
    stand = qclib.TableOrWall(cs)
    idname = '_'+stand

    labvals = qclib.ReportEntries(cs)
    tmpdict={}
    for elem in labvals:
        #labvals.append( {'name':'label','value':0, 'quantity':'columnname','level':'1:default, 2: detail','pos':missing or a number} )
        varname = elem['name']+str(idname)
        results.addFloat(varname, elem['value'])

    ## 6. Build artefact picture thumbnail
    filename = 'test'+idname+'.jpg' # Use jpg if a thumbnail is desired

    qclib.saveAnnotatedImage(cs,filename)
    varname = 'AnnotatedImage'+idname
    results.addObject(varname, filename)

def acqdatetime_series(data, results, action):
    """
    Read acqdatetime from dicomheaders and write to IQC database

    Workflow:
        1. Read only headers
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    dt = wadwrapper_lib.acqdatetime_series(dcmInfile)

    results.addDateTime('AcquisitionDateTime', dt) 


def header_series(data, results, action):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Read only headers
        2. Run tests
        3. Build xml output
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    info = 'qcwad'

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    ## 2. Run tests
    qclib = QCXRay_lib.XRayQC()
    room = _getRoomDefinition(params)

    ## Table or Wall? from distances and sensitivity; for well defined protocols to be defined in DESCRIPTION field
    cs = QCXRay_lib.XRayStruct(dcmInfile,None,room)
    cs.verbose = False # do not produce detailed logging
    override_settings(cs, params)

    dicominfo = qclib.DICOMInfo(cs,info)
    idname = '_'+qclib.TableOrWall(cs)

    ## 3. Build xml output
    floatlist = [
        'Exposure (mAs)',
        'DistanceSourceToDetector (mm)',
        'ExposureTime (ms)',
        'ImageAreaDoseProduct',
        'Sensitivity',
        'kVp'
    ]
    offset = -25
    varname = 'pluginversion'+idname
    results.addString('pluginversion'+idname, str(qclib.qcversion)) 
    for elem in dicominfo:
        varname = elem['name']+str(idname)
        if elem['name'] in floatlist:
            results.addFloat(varname, elem['value'])
        else:
            results.addString(varname, str(elem['value'])[:min(len(str(elem['value'])),100)])
    
    varname = 'room'+idname
    results.addString(varname, cs.forceRoom.name) 
    varname = 'stand'+idname
    results.addString(varname, qclib.TableOrWall(cs))

if __name__ == "__main__":
    data, results, config = pyWADinput()

    # read runtime parameters for module
    for name,action in config['actions'].items():
        if name == 'acqdatetime':
            acqdatetime_series(data, results, action)

        elif name == 'header_series':
            header_series(data, results, action)
        
        elif name == 'qc_series':
            qc_series(data, results, action)

    #results.limits["minlowhighmax"]["mydynamicresult"] = [1,2,3,4]

    results.write()
