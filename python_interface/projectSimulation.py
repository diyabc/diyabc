# -*- coding: utf-8 -*-

import time
import os
import subprocess
from project import Project
#from subprocess import Popen, PIPE, STDOUT 
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import uic
#from uis.project_ui import *
from setHistoricalModelSimulation import SetHistoricalModelSimulation
from setGenDataSimulation import SetGeneticDataSimulation
#from mutationModel.setMutationModelMsat import SetMutationModelMsat
#from mutationModel.setMutationModelSequences import SetMutationModelSequences
#from summaryStatistics.setSummaryStatisticsMsat import SetSummaryStatisticsMsat
#from summaryStatistics.setSummaryStatisticsSeq import SetSummaryStatisticsSeq
import os.path
import output
from utils.cbgpUtils import log

class ProjectSimulation(Project):
    """ classe qui représente un projet de simulation
    """
    def __init__(self,name,dir=None,parent=None):
        super(ProjectSimulation,self).__init__(name,dir,parent)

        self.dico_loc_nb = None
        self.sexRatio = None

        self.ui.projNameLabel.setText("Data file generic name :")
        self.ui.label.setText("Target directory :")
        self.ui.dirEdit.setText("%s"%self.dir)
        self.ui.label_4.hide()
        self.ui.browseDataFileButton.hide()
        self.ui.dataFileEdit.hide()
        self.ui.label_10.hide()
        self.ui.dataFileInfoLabel.hide()

        self.ui.label_8.hide()
        self.ui.nbSetsDoneEdit.hide()

        self.ui.groupBox_6.show()
        self.ui.groupBox_7.show()
        self.ui.groupBox_8.show()

        self.ui.removeTab(1)
        self.ui.setTabText(0,QString("Simulate data sets"))

        self.ui.setHistoricalButton.setDisabled(False)
        self.ui.setGeneticButton.setDisabled(False)

        self.hist_model_win = SetHistoricalModelSimulation(self)
        self.hist_model_win.hide()

        self.locusNumberFrame = uic.loadUi("uis/setLocusNumber.ui")
        self.locusNumberFrame.parent = self
        QObject.connect(self.locusNumberFrame.okButton,SIGNAL("clicked()"),self.checkSampleNSetGenetic)
        QObject.connect(self.locusNumberFrame.mxEdit,SIGNAL("textChanged(QString)"),self.checkXYValues)
        QObject.connect(self.locusNumberFrame.myEdit,SIGNAL("textChanged(QString)"),self.checkXYValues)
        QObject.connect(self.locusNumberFrame.sxEdit,SIGNAL("textChanged(QString)"),self.checkXYValues)
        QObject.connect(self.locusNumberFrame.syEdit,SIGNAL("textChanged(QString)"),self.checkXYValues)
        QObject.connect(self.locusNumberFrame.mhEdit,SIGNAL("textChanged(QString)"),self.haploidChanged)
        QObject.connect(self.locusNumberFrame.shEdit,SIGNAL("textChanged(QString)"),self.haploidChanged)
        self.locusNumberFrame.mxEdit.setText('1')
        self.locusNumberFrame.mxEdit.setText('0')

        self.gen_data_win = SetGeneticDataSimulation(self)
        self.gen_data_win.hide()

        #self.setHistValid(False)
        #self.setGenValid(False)
        #self.connect(self.ui.runReftableButton, SIGNAL("clicked()"),self.runSimulation)

    def stopUiGenReftable(self):
        self.ui.runReftableButton.setText("Run computations")
        self.ui.runReftableButton.setDisabled(False)
        self.ui.progressBar.hide()

    def haploidChanged(self,srt):
        try:
            mh = int(str(self.locusNumberFrame.mhEdit.text()))
            sh = int(str(self.locusNumberFrame.shEdit.text()))
            if sh > 0 or mh > 0:
                for ed in [self.locusNumberFrame.sxEdit,self.locusNumberFrame.mxEdit,self.locusNumberFrame.myEdit,self.locusNumberFrame.syEdit,self.locusNumberFrame.maEdit,self.locusNumberFrame.saEdit]:
                    ed.setDisabled(True)
                    ed.setText("0")
            else:
                for ed in [self.locusNumberFrame.sxEdit,self.locusNumberFrame.mxEdit,self.locusNumberFrame.myEdit,self.locusNumberFrame.syEdit,self.locusNumberFrame.maEdit,self.locusNumberFrame.saEdit]:
                    ed.setDisabled(False)
        except Exception,e:
            pass


    def checkXYValues(self,val):
        try:
            mx = int(str(self.locusNumberFrame.mxEdit.text()))
            my = int(str(self.locusNumberFrame.myEdit.text()))
            sx = int(str(self.locusNumberFrame.sxEdit.text()))
            sy = int(str(self.locusNumberFrame.syEdit.text()))
            if mx <= 0 and my <= 0 and sx <= 0 and sy <= 0:
                self.locusNumberFrame.sexRatioEdit.setDisabled(True)
            else:
                self.locusNumberFrame.sexRatioEdit.setDisabled(False)
        except Exception,e:
            self.locusNumberFrame.sexRatioEdit.setDisabled(True)


    def save(self):
        pass

    def setGenetic(self):
        """ passe sur l'onglet correspondant
        """
        log(1,"Entering in Genetic Data Setting")
        self.ui.refTableStack.addWidget(self.locusNumberFrame)
        self.ui.refTableStack.setCurrentWidget(self.locusNumberFrame)
        self.sexRatio = None
        self.setGenValid(False)

    def checkSampleNSetGenetic(self):
        dico_loc_nb = {}
        editList = self.locusNumberFrame.findChildren(QLineEdit)
        editList.remove(self.locusNumberFrame.sexRatioEdit)
        nbXY = 0
        nbHap = 0
        nbDip = 0
        nbpos = 0
        try:
            for le in editList:
                dico_loc_nb[str(le.objectName())[:2]] = int(le.text())
                if int(le.text()) < 0:
                    output.notify(self,"Number of sample error","Number of sample must be positive integers")
                    return
                elif int(le.text()) > 0:
                    nbpos += 1
                    if "x" in str(le.objectName())[:2] or "y" in str(le.objectName())[:2]:
                        nbXY += 1
                    elif "d" in str(le.objectName())[:2]:
                        nbDip += 1
                    elif "h" in str(le.objectName())[:2]:
                        nbHap += 1
            if nbpos == 0:
                output.notify(self,"Number of sample error","You must have at leat one sample")
                return
        except Exception,e:
            output.notify(self,"Number of sample error","Input error : \n%s"%e)
            return
        if nbXY > 0:
            try:
                val = float(str(self.locusNumberFrame.sexRatioEdit.text()))
                if val < 0 or val > 100:
                    output.notify(self,"Value error","Sex ratio value must be in [0,1]")
                    return
                self.sexRatio = val
            except Exception,e:
                output.notify(self,"Value error","Sex ratio value must be in [0,1]")
                return

        # verification : a-t-on modifié le nb de loci ?
        if self.dico_loc_nb != None:
            changed = False
            for key in self.dico_loc_nb:
                if self.dico_loc_nb[key] != dico_loc_nb[key]:
                    changed = True
                    break
        else:
            changed = True

        self.dico_loc_nb = dico_loc_nb

        log(3,"Numbers of locus : %s"%dico_loc_nb)
        if changed:
            # on vide les gen data
            self.gen_data_win = SetGeneticDataSimulation(self)
            self.gen_data_win.hide()
            self.gen_data_win.fillLocusTable(dico_loc_nb)

        self.ui.refTableStack.removeWidget(self.ui.refTableStack.currentWidget())
        self.ui.refTableStack.addWidget(self.gen_data_win)
        self.ui.refTableStack.setCurrentWidget(self.gen_data_win)


    def returnToMe(self):
        self.ui.refTableStack.removeWidget(self.ui.refTableStack.currentWidget())
        self.ui.refTableStack.setCurrentIndex(0)

    def writeHeaderSim(self):
        #if self.verifyRefTableValid():
        # première ligne
        sexRatioTxt = ""
        if self.sexRatio != None:
            sexRatioTxt = self.sexRatio
        nb_rec_edit = str(self.ui.nbSetsReqEdit.text())
        if nb_rec_edit.isdigit() and int(nb_rec_edit) > 0:
            nb_rec = int(nb_rec_edit)
        else:
            output.notify(self,"Value error","Required number of simulated data sets must be a positive integer")
            return
        print "%s %s %s"%(self.name,nb_rec,sexRatioTxt)
        print self.hist_model_win.getConf()
        print ""
        print self.gen_data_win.getConf().replace(u'\xb5','u')
        fdest = open("%s/headersim.txt"%self.dir,"w")
        fdest.write("%s %s %s\n"%(self.name,nb_rec,sexRatioTxt))
        fdest.write("%s\n\n"%self.hist_model_win.getConf())
        fdest.write("%s"%self.gen_data_win.getConf().replace(u'\xb5','u'))
        fdest.close()

    @pyqtSignature("")
    def on_btnStart_clicked(self):
        self.writeHeaderSim()
        try:
            nb_to_gen = int(self.ui.nbSetsReqEdit.text())
        except Exception,e:
            output.notify(self,"value error","Check the value of required number of data sets\n\n%s"%e)
            return
        self.th = SimulationThread(self,nb_to_gen)
        self.th.connect(self.th,SIGNAL("increment"),self.incProgress)
        self.th.connect(self.th,SIGNAL("simulationProblem"),self.simulationProblem)
        self.th.connect(self.th,SIGNAL("simulationLog"),self.refTableLog)
        #self.ui.progressBar.connect (self, SIGNAL("canceled()"),self.th,SLOT("cancel()"))
        self.th.start()

    def simulationProblem(self):
        output.notify(self,"Simulation problem","Something happened during the simulation :\n %s"%(self.th.problem))
        self.stopSimulation()

    def stopSimulation(self):
        if self.th != None:
            self.th.terminate()
            self.th.killProcess()
            self.th = None
        if os.path.exists("%s/simulation.out"%(self.dir)):
            os.remove("%s/simulation.out"%(self.dir))
        log(1,"Simulation stopped")


class SimulationThread(QThread):
    """ thread de traitement qui met à jour la progressBar en fonction de l'avancée de
    la génération de la reftable
    """
    def __init__(self,parent,nb_to_gen):
        super(SimulationThread,self).__init__(parent)
        self.parent = parent
        self.nb_to_gen = nb_to_gen
        self.processus = None

        self.logmsg = ""
        self.loglvl = 3

    def log(self,lvl,msg):
        """ evite de manipuler les objets Qt dans un thread non principal
        """
        self.loglvl = lvl
        self.logmsg = msg
        self.emit(SIGNAL("simulationLog"))

    def killProcess(self):
        self.log(3,"Attempting to kill simulation process")
        if self.processus != None:
            if self.processus.poll() == None:
                self.processus.kill()
                self.log(3,"Killing simulation process (pid:%s) DONE"%(self.processus.pid))

    def run (self):
        # lance l'executable
        #outfile = os.path.expanduser("~/.diyabc/general.out")
        outfile = "%s/simulation.out"%self.parent.dir
        if os.path.exists(outfile):
            os.remove(outfile)
        if os.path.exists("%s/progress.txt"%self.parent.dir):
            os.remove("%s/progress.txt"%self.parent.dir)
        fg = open(outfile,"w")
        try:
            self.log(2,"Running the executable for the simulation")
            exPath = self.parent.parent.preferences_win.getExecutablePath()
            nbMaxThread = self.parent.parent.preferences_win.getMaxThreadNumber()
            cmd_args_list = [exPath,"-p", "%s/"%self.parent.dir, "-k", "-m", "-t", "%s"%nbMaxThread]
            #print " ".join(cmd_args_list)
            self.log(3,"Command launched : %s"%" ".join(cmd_args_list))
            p = subprocess.Popen(cmd_args_list, stdout=fg, stdin=subprocess.PIPE, stderr=subprocess.STDOUT) 
            self.processus = p
        except Exception,e:
            #print "Cannot find the executable of the computation program %s"%e
            self.problem = "Problem during program launch \n%s"%e
            self.emit(SIGNAL("simulationProblem"))
            #output.notify(self.parent(),"computation problem","Cannot find the executable of the computation program")
            fg.close()
            return

        # boucle toutes les secondes pour verifier les valeurs dans le fichier
        finished = False
        while not finished:
            time.sleep(1)
            # verification de l'arret du programme
            if p.poll() != None:
                # lecture 
                if os.path.exists("%s/progress.txt"%(self.parent.dir)):
                    #print 'open'
                    f = open("%s/progress.txt"%(self.parent.dir),"r")
                    lines = f.readlines()
                    f.close()
                    finished = True
                    self.log(2,"Simulation terminated normaly")
                else:
                    fg.close()
                    fout = open(outfile,'r')
                    lastline = fout.readlines()
                    if len(lastline) > 0:
                        lastline = lastline[-1]
                    fout.close()
                    self.problem = "Simulation program exited anormaly\n%s"%lastline
                    self.emit(SIGNAL("simulationProblem"))
                    return
                    
                fg.close()
        fg.close()
