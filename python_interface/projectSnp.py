# -*- coding: utf-8 -*-

import os
import codecs
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from projectReftable import ProjectReftable
from summaryStatistics.setSummaryStatisticsSnp import SetSummaryStatisticsSnp
from utils.data import DataSnp
import os.path
import output
from utils.cbgpUtils import log

class ProjectSnp(ProjectReftable):
    """ classe qui représente un projet de simulation
    """
    def __init__(self,name,dir=None,parent=None):
        super(ProjectSnp,self).__init__(name,dir,parent)

        self.ui.nbMicrosatLabel.hide()
        self.ui.nbSequencesLabel.hide()
        # TODO suppr cette ligne
        self.setGenValid(False)
        #QObject.connect(self.ui.setSumSnpButton,SIGNAL("clicked()"),self.setSumStat)
        self.ui.label_15.setText("Summary statistics")

        self.typesOrdered = ["A","H","X","Y","M"]

        self.ui.frame_11.show()

    def setGenetic(self):
        """ initie la définition des summary statistics
        """
        log(1,"Entering in Summary statistics")
        ty = str(self.typeCombo.currentText())
        self.ui.refTableStack.addWidget(self.sum_stat_wins[ty])
        self.ui.refTableStack.setCurrentWidget(self.sum_stat_wins[ty])
        self.setGenValid(False)

    def getNbSumStats(self):
        nb = 0
        for ty in self.sum_stat_wins.keys():
            (nstat,stat_txt) = self.sum_stat_wins[ty].getSumConf()
            nb += int(nstat)
        return nb

    def updateNbStats(self):
        nb = self.getNbSumStats()
        if nb > 1:
            plur='s'
            self.setGenValid(True)
        else:
            plur=''
        self.ui.nbSumStatsLabel.setText("Total : %s summary statistic%s"%(nb,plur))


    #def setGenValid(self,valid):
    #    """ met à jour l'état des genetic data
    #    et change l'icone du bouton en fonction de la validité
    #    """
    #    self.gen_state_valid = valid
    #    self.verifyRefTableValid()
    #    if valid:
    #        self.ui.setSumSnpButton.setIcon(QIcon("docs/icons/ok.png"))
    #    else:
    #        self.ui.setSumSnpButton.setIcon(QIcon("docs/icons/redcross.png"))

    def checkAll(self):
        """ vérification du modèle historique et mutationnel
        cette fonction est appelée au chargement du projet pour restituer l'etat du projet
        """
        log(2,"Checking validity of Historical Model and Sum Stats")
        # historical model : 
        self.hist_model_win.definePriors(silent=True)
        if self.hist_model_win.checkAll(silent=True):
            self.setHistValid(True)
            self.hist_model_win.majProjectGui()
        # mutation model : plus facile d'utiliser directement la validation
        for ty in self.sum_stat_wins.keys():
            self.sum_stat_wins[ty].validate(silent=True)

    def returnToMe(self):
        self.ui.refTableStack.removeWidget(self.ui.refTableStack.currentWidget())
        self.ui.refTableStack.setCurrentIndex(0)

    def loadSumStatsConf(self):
        if os.path.exists(self.dir):
            if os.path.exists("%s/%s"%(self.dir,self.parent.gen_conf_name)):
                f = codecs.open("%s/%s"%(self.dir,self.parent.gen_conf_name),"r","utf-8")
                lines = f.readlines()
                f.close()
                nl = 1
                diconb = {}
                dicoFrom = {}
                dicoGrTy = {}
                dicoStatLines = {}
                while (lines[nl].strip() != ""):
                    ty = lines[nl].split('<')[1].split('>')[0]
                    numGr = int(lines[nl].strip().split(' ')[2][1:])
                    diconb[ty] = lines[nl].strip().split(' ')[0]
                    dicoFrom[ty] = lines[nl].strip().split(' ')[4]
                    dicoGrTy[numGr] = ty
                    self.sum_stat_wins[ty].ui.takenEdit.setText(diconb[ty])
                    self.sum_stat_wins[ty].ui.fromEdit.setText(dicoFrom[ty])
                    nl += 1
                nl += 2
                while nl < len(lines) and lines[nl].strip() != "":
                    numGroup = int(lines[nl].strip().split(' ')[1][1:])
                    ty = dicoGrTy[numGroup]
                    statlines = []
                    nl += 1
                    while nl < len(lines) and "group" not in lines[nl]:
                        statlines.append( lines[nl] )
                        nl += 1
                    dicoStatLines[ty] = statlines

                for ty in dicoStatLines.keys():
                    self.sum_stat_wins[ty].setSumConf(dicoStatLines[ty])

    def loadFromDir(self):
        """ charge les infos à partir du répertoire self.dir
        """
        log(2,"Launching load procedures")
        # GUI
        self.ui.dirEdit.setText(self.dir)
        self.ui.browseDataFileButton.setDisabled(True)
        self.ui.browseDataFileButton.hide()
        self.ui.browseDirButton.hide()
        #self.ui.groupBox.show()
        self.ui.groupBox_6.show()
        self.ui.groupBox_7.show()
        self.ui.groupBox_8.show()
        self.ui.setHistoricalButton.setDisabled(False)
        self.ui.setGeneticButton.setDisabled(False)

        # lecture du meta project
        if self.loadMyConf():
            # lecture de conf.hist.tmp
            self.hist_model_win.loadHistoricalConf()
            self.loadSumStatsConf()
            self.loadAnalysis()
        else:
            raise Exception("Impossible to read the project configuration")
            output.notify(self,"Load error","Impossible to read the project configuration")

    def loadMyConf(self):
        """ lit le fichier conf.tmp pour charger le fichier de données
        """
        log(2,"Reading '%s' to load datafile"%self.parent.main_conf_name)
        if os.path.exists(self.ui.dirEdit.text()+"/%s"%self.parent.main_conf_name):
            f = open("%s/%s"%(self.dir,self.parent.main_conf_name),"r")
            lines = f.readlines()
            self.dataFileName = lines[0].strip()
            self.ui.dataFileEdit.setText(lines[0].strip())
            # lecture du dataFile pour les infos de Gui Projet
            if self.loadDataFile("%s/%s"%(self.dir,lines[0].strip())):
                return True
            ## comme on a lu le datafile, on peut remplir le tableau de locus dans setGeneticData
            #self.gen_data_win.fillLocusTableFromData()
        return False

    def loadDataFile(self,name):
        """ Charge le fichier de données passé en paramètre. Cette fonction est appelée lors
        de l'ouverture d'un projet existant et lors du choix du fichier de données pour un nouveau projet
        """
        log(2,"Loading datafile '%s'"%name)

        try:
            self.data = DataSnp(name)
            typestr = ""
            for ty in self.data.ntypeloc.keys():
                typestr += " %s : %s,"%(ty,self.data.ntypeloc[ty])
            typestr = typestr[:-1]
            self.ui.dataFileInfoLabel.setText("%s loci (SNP)\n%s individuals in %s samples\n%s" % (self.data.nloc,self.data.nindtot,self.data.nsample,typestr))
            self.ui.dataFileEdit.setText(name)
            self.dataFileSource = name
            self.ui.browseDirButton.setDisabled(False)

        except Exception,e:
            keep = ""
            if self.ui.dataFileEdit.text() != "":
                #keep = "\n\nKeeping previous selected file"
                keep = "\n\nThe file was not loaded, nothing was changed"
            output.notify(self,"Data file error","%s%s"%(e,keep))
            return False

        # on declare les sumstats apres avoir chargé le datafile car c'est nécessaire
        # feinte pour que le parent.parent renvoie au projet
        self.dummy = QFrame()
        self.dummy.parent = self
        self.sum_stat_wins = {}
        for ty in self.data.ntypeloc.keys():
            self.sum_stat_wins[ty] = SetSummaryStatisticsSnp(parent=self.dummy,numGroup=ty)
            self.sum_stat_wins[ty].hide()

        # selection du type de snp pour sumstats
        self.typeCombo = QComboBox(self)
        for ty in self.typesOrdered:
            if ty in self.data.ntypeloc.keys():
                self.typeCombo.addItem(ty)
        self.ui.horizontalLayout_6.addWidget(QLabel("for locus type :"))
        self.ui.horizontalLayout_6.addWidget(self.typeCombo)

        return True

    def save(self):
        """ sauvegarde du projet -> mainconf, histconf, genconf, theadconf
        Si le gen et hist sont valides, on génère le header
        """
        #print "je me save"
        log(2,"Saving project '%s'"%self.dir)
        self.parent.showStatus("Saving project %s"%self.name)

        if self.dir != None and self.dataFileName != "":
            # save meta project
            if os.path.exists(self.dir+"/%s"%self.parent.main_conf_name):
                os.remove("%s/%s"%(self.dir,self.parent.main_conf_name))

            f = codecs.open(self.dir+"/%s"%self.parent.main_conf_name,'w',"utf-8")
            f.write("%s\n"%self.dataFileName)
            # recup du nombre de params (depuis historical model et les mutational qui varient)
            nb_param = self.hist_model_win.getNbParam()
            #nb_param += self.gen_data_win.getNbParam()
            #nb_sum_stats = self.gen_data_win.getNbSumStats()
            f.write("%s parameters and %s summary statistics\n\n"%(nb_param,self.getNbSumStats()))
            f.close()

            # save hist conf
            self.hist_model_win.writeHistoricalConfFromGui()
            # save gen conf
            self.writeGeneticConfFromGui()
            # save th conf et production du reftable header
            if self.gen_state_valid and self.hist_state_valid:
                self.writeThConf()
                self.writeRefTableHeader()
            self.saveAnalysis()
            self.parent.clearStatus()
            self.parent.showStatus("Project %s successfully saved"%self.name,2000)
        else:
            output.notify(self,"Saving is impossible","Project %s is not yet completly created"%self.name)
            self.parent.clearStatus()

    def writeThConf(self):
        """ ecrit le header du tableau de resultat qui sera produit par les calculs
        il contient, les paramètres historicaux,  les summary stats
        """
        log(2,"Writing last part of the header (the parameter table header) in %s"%self.parent.table_header_conf_name)
        hist_params_txt = self.hist_model_win.getParamTableHeader()
        sum_stats_txt = self.getSumStatsTableHeader()
        if os.path.exists(self.dir+"/%s"%self.parent.table_header_conf_name):
            os.remove("%s/%s"%(self.dir,self.parent.table_header_conf_name))
        f = codecs.open(self.dir+"/%s"%self.parent.table_header_conf_name,'w',"utf-8")
        f.write("scenario%s%s"%(hist_params_txt,sum_stats_txt))
        f.close()
        
    def getSumStatsTableHeader(self):
        """ retourne la partie sumstats du table header
        """
        result = u""
        numGroup = 1
        for ty in self.typesOrdered:
            if ty in self.sum_stat_wins.keys():
                sums_txt = self.sum_stat_wins[ty].getSumStatsTableHeader(numGroup)
                result += sums_txt
                numGroup += 1
        return result

    def writeGeneticConfFromGui(self):
        locidesc = ""
        statsdesc = ""
        numGr = 1
        totNbStat = 0
        for ty in self.typesOrdered:
            if ty in self.sum_stat_wins.keys():
                (nstat,statstr) = self.sum_stat_wins[ty].getSumConf()
                totNbStat += nstat
                if nstat > 0:
                    locidesc += "%s <%s> G%s from %s\n"%(self.sum_stat_wins[ty].ui.takenEdit.text(),ty,numGr,self.sum_stat_wins[ty].ui.fromEdit.text())
                    statsdesc += "group G%s (%s)\n%s"%(numGr,nstat,statstr)
                    numGr += 1
        res = "loci description (%s)\n"%(numGr - 1)
        res += locidesc

        res += "\ngroup summary statistics (%s)\n"%(totNbStat)
        res += statsdesc
        res += "\n"

        print res
        if os.path.exists(self.dir+"/%s"%self.parent.gen_conf_name):
            os.remove("%s/%s" % (self.dir,self.parent.gen_conf_name))

        f = codecs.open(self.dir+"/%s"%self.parent.gen_conf_name,'w',"utf-8")
        f.write(res)
        f.close()

    def freezeGenData(self,yesno=True):
        # TODO
        pass