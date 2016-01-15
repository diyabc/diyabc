# -*- coding: utf-8 -*-

import sys
import time
import os
import itertools
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui

class SetSummaryStatistics(QFrame):
    """ ecran de selection des summary statistics
    superclasse abstraite regroupant les méthodes communes aux setsumstats msat et seq
    les deux classes filles n'ont pas la même interface graphique, l'import de l'UI n'est donc pas
    effectué ici
    """
    def __init__(self,parent=None,box_group=None,numGroup=0):
        super(SetSummaryStatistics,self).__init__(parent)
        self.parent=parent
        self.box_group = box_group
        self.numGroup = numGroup
        self.oneSampleList = []
        self.twoSampleList = []
        self.admixSampleList = []

        self.createWidgets()

        # remplissage des sample
        for i in range(self.parent.parent.data.nsample):
            self.addOneSample(i+1)

        if self.parent.parent.data.nsample >= 2:
            for i in range(self.parent.parent.data.nsample):
                j = i+1
                while j < self.parent.parent.data.nsample:
                    self.addTwoSample(i+1,j+1)
                    j+=1
            if self.parent.parent.data.nsample >= 3:
                pass
            else:
                self.ui.admixFrame.hide()
        else:
            self.ui.twoFrame.hide()
            self.ui.admixFrame.hide()

        for i in range(self.parent.parent.data.nsample):
            self.ui.adm1Combo.addItem(str(i+1))
            self.ui.adm2Combo.addItem(str(i+1))
            self.ui.adm3Combo.addItem(str(i+1))

    def allPressed(self):
        """ clic sur bouton 'all' ou 'none'. analyse les informations pour
        savoir ce qu'il faut cocher ou décocher
        """
        button_name = str(self.sender().objectName())
        checkname = button_name.replace('Button','').replace('all','').replace('none','').lower()
        if checkname == "ci2":
            self.checkAll(('all' in button_name),"%sRightCheck"%checkname)
            self.checkAll(('all' in button_name),"%sLeftCheck"%checkname)
        else:
            self.checkAll(('all' in button_name),"%sCheck"%checkname)

    def checkAll(self,yesno,objname):
        """ coche ou décoche toute une ligne
        """
        for chkbox in self.findChildren(QCheckBox,objname):
            chkbox.setChecked(yesno)

    def addAdmixSampleGui(self,num1,num2,num3):
        """ methode d'ajout d'un bloc dans 'three sample admix sum stats'
        """
        frame_12 = QtGui.QFrame(self.ui.admixFrame)
        frame_12.setMinimumSize(QtCore.QSize(55, 0))
        frame_12.setMaximumSize(QtCore.QSize(55, 16777215))
        frame_12.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_12.setFrameShadow(QtGui.QFrame.Raised)
        frame_12.setObjectName("frame_12")
        verticalLayout_6 = QtGui.QVBoxLayout(frame_12)
        verticalLayout_6.setSpacing(1)
        #verticalLayout_6.setAlignment(Qt.AlignHCenter)
        verticalLayout_6.setContentsMargins(1, 1, 1, 1)
        verticalLayout_6.setObjectName("verticalLayout_6")
        threeSampleLabel = QtGui.QLabel("Samp\n%i&%i&%i"%(num1,num2,num3),frame_12)
        font = QtGui.QFont()
        font.setPointSize(8)
        threeSampleLabel.setFont(font)
        threeSampleLabel.setObjectName("threeSampleLabel")
        threeSampleLabel.setMinimumSize(QtCore.QSize(0, 48))
        threeSampleLabel.setMaximumSize(QtCore.QSize(16777215, 48))
        verticalLayout_6.addWidget(threeSampleLabel)
        rmButton = QtGui.QPushButton("remove",frame_12)
        font = QtGui.QFont()
        font.setPointSize(8)
        rmButton.setFont(font)
        rmButton.setObjectName("rmAdmixButton")
        verticalLayout_6.addWidget(rmButton)
        horizontalLayout_17 = QtGui.QHBoxLayout()
        horizontalLayout_17.setObjectName("horizontalLayout_17")
        spacerItem20 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_17.addItem(spacerItem20)
        ml3Check = QtGui.QCheckBox(frame_12)
        ml3Check.setMinimumSize(QtCore.QSize(20, 0))
        ml3Check.setMaximumSize(QtCore.QSize(20, 16777215))
        ml3Check.setObjectName("ml3Check")
        ml3Check.setChecked(True)
        ml3Check.setDisabled(True)
        horizontalLayout_17.addWidget(ml3Check)
        spacerItem21 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_17.addItem(spacerItem21)
        verticalLayout_6.addLayout(horizontalLayout_17)
        self.ui.horizontalLayout_2.addWidget(frame_12)

        self.admixSampleList.append(frame_12)

        QObject.connect(rmButton,SIGNAL("clicked()"),self.rmAdmix)

    def rmAdmix(self):
        """ supprime un admix sample
        """
        box = self.sender().parent()
        box.hide()
        self.admixSampleList.remove(box)

    def addAdmix(self):
        """ ajoute un admix sample en fonction des valeurs des combo box
        """
        val1 = str(self.ui.adm1Combo.currentText())
        val2 = str(self.ui.adm2Combo.currentText())
        val3 = str(self.ui.adm3Combo.currentText())
        # verif de la coherence
        if val1 == val2 or val2 == val3 or val3 == val1:
            QMessageBox.information(self,"Value error","Each value must be different to add an admixture summary statistic")
        else:
            # verif qu'il n'existe pas d'équivalent
            ok = True
            for box in self.admixSampleList:
                txt = str(box.findChild(QLabel,"threeSampleLabel").text()).replace('Samp\n','')
                valbox = (txt.split('&')[0],txt.split('&')[1],txt.split('&')[2])
                if valbox == (val1,val2,val3):
                    ok = False
                    break
            if not ok:
                QMessageBox.information(self,"Value error","This admixture summary statistic is already set")
            else:
                self.addAdmixSampleGui(int(val1),int(val2),int(val3))

    def allAdmix(self):
        """ add all admix sample combinations
        """
        self.noneAdmix()
        for i in range(self.parent.parent.data.nsample):
            s = set(range(self.parent.parent.data.nsample))
            s.remove(i)
            for j in itertools.combinations(s, 2) :
                self.addAdmixSampleGui(i+1, j[0]+1, j[1]+1)

    def noneAdmix(self):
        """ remove all admix samples
        """
        while len(self.admixSampleList) > 0 :
            for box in self.admixSampleList:
                box.hide()
                self.admixSampleList.remove(box)

    def getSumConf(self):
        """ retourne un tuple contenant le nombre de statistiques demandés
        et la configuration à écrire
        """
        (nstat,dico_stats) = self.getStats()
        conf_txt = ""
        #for k in dico_stats.keys():
        for k in self.statList:
            if len(dico_stats[k]) > 0:
                conf_txt += "%s "%k
                for lsample in dico_stats[k]:
                    conf_txt += "%s "%lsample
                conf_txt += "\n"
        return (nstat,conf_txt)

    def exit(self):
        ## reactivation des onglets
        #self.parent.parent.setTabEnabled(self.parent.parent.indexOf(self.parent),True)
        #self.parent.parent.removeTab(self.parent.parent.indexOf(self))
        #self.parent.parent.setCurrentIndex(self.parent.parent.indexOf(self.parent))
        self.parent.parent.ui.refTableStack.removeWidget(self)
        self.parent.parent.ui.refTableStack.setCurrentWidget(self.parent)
        self.parent.majProjectGui(ss=self.parent.getNbSumStats())

    def validate(self,silent=False):
        self.exit()
        if self.box_group != None:
            (nstat,stat_txt) = self.getSumConf()
            ## on met le bouton en police normale pour signaler qu'il est valide
            #set_sum_button = self.box_group.findChild(QPushButton,"setSumButton")
            #fontt = set_sum_button.font()
            #if nstat > 0:
            #    fontt.setBold(False)
            #    set_sum_button.setStyleSheet("border-color: #000000")
            #else:
            #    fontt.setBold(True)
            #    set_sum_button.setStyleSheet("background-color: #EFB1B3")
            #set_sum_button.setFont(fontt)


