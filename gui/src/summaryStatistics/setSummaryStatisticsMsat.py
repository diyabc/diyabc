# -*- coding: utf-8 -*-

import sys, codecs
import time
import os
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
from PyQt4 import uic
#from uis.setSummaryStatisticsMsat_ui import Ui_Frame
from setSummaryStatistics import SetSummaryStatistics
import output
import variables
from variables import UIPATH
from utils.cbgpUtils import getFsEncoding

formSetSummaryStatisticsMsat,baseSetSummaryStatisticsMsat = uic.loadUiType((u"%s/setSummaryStatisticsMsat.ui"%UIPATH).encode(getFsEncoding(logLevel=False)))

class SetSummaryStatisticsMsat(SetSummaryStatistics,formSetSummaryStatisticsMsat,baseSetSummaryStatisticsMsat):
    """ ecran de selection des summary statistics pour les microsats
    """
    def __init__(self,parent=None,box_group=None,numGroup=0):
        super(SetSummaryStatisticsMsat,self).__init__(parent,box_group,numGroup)

        self.statList = ["NAL","HET","VAR","MGW","N2P","H2P","V2P","FST","LIK","DAS","DM2","AML"]

    def createWidgets(self):
        self.ui=self
        self.ui.setupUi(self)

        QObject.connect(self.ui.exitButton,SIGNAL("clicked()"),self.exit)
        QObject.connect(self.ui.clearButton,SIGNAL("clicked()"),self.clear)
        QObject.connect(self.ui.okButton,SIGNAL("clicked()"),self.validate)

        # all buttons
        QObject.connect(self.ui.allMnoaButton,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allMgdButton,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allMsvButton,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allMgwButton,SIGNAL("clicked()"),self.allPressed)
        # none buttons
        QObject.connect(self.ui.noneMnoaButton,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneMgdButton,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneMsvButton,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneMgwButton,SIGNAL("clicked()"),self.allPressed)

        # all buttons 2
        QObject.connect(self.ui.allMnoa2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allMgd2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allMsv2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allFst2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allCi2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allSad2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.allDd2Button,SIGNAL("clicked()"),self.allPressed)
        #none buttons 2
        QObject.connect(self.ui.noneMnoa2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneMgd2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneMsv2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneFst2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneCi2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneSad2Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneDd2Button,SIGNAL("clicked()"),self.allPressed)
        # all and none for admixtures
        QObject.connect(self.ui.allMl3Button,SIGNAL("clicked()"),self.allPressed)
        QObject.connect(self.ui.noneMl3Button,SIGNAL("clicked()"),self.allPressed)

        QObject.connect(self.ui.addAdmixButton,SIGNAL("clicked()"),self.addAdmix)
        QObject.connect(self.ui.allAdmixButton,SIGNAL("clicked()"),self.allAdmix)
        QObject.connect(self.ui.noneAdmixButton,SIGNAL("clicked()"),self.noneAdmix)

        self.ui.horizontalLayout_4.setAlignment(Qt.AlignLeft)
        self.ui.horizontalLayout_2.setAlignment(Qt.AlignLeft)
        self.ui.horizontalLayout_5.setAlignment(Qt.AlignLeft)

        self.ui.frame_12.hide()
        self.ui.frame_11.hide()
        self.ui.frame_9.hide()
        self.ui.allMl3Button.hide()
        self.ui.noneMl3Button.hide()

    def addOneSample(self,num):
        """ methode d'ajout d'un bloc dans 'one sample sum stats'
        """
        frame_9 = QtGui.QFrame(self.ui.oneFrame)
        frame_9.setMinimumSize(QtCore.QSize(60, 0))
        frame_9.setMaximumSize(QtCore.QSize(60, 16777215))
        frame_9.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_9.setFrameShadow(QtGui.QFrame.Raised)
        frame_9.setObjectName("frame_9")
        verticalLayout_4 = QtGui.QVBoxLayout(frame_9)
        verticalLayout_4.setObjectName("verticalLayout_4")
        verticalLayout_4.setSpacing(1)
        verticalLayout_4.setContentsMargins(1, 1, 1, 1)
        oneSampleLabel = QtGui.QLabel("Samp %i"%num,frame_9)
        font = QtGui.QFont()
        font.setPointSize(8)
        oneSampleLabel.setFont(font)
        oneSampleLabel.setMinimumSize(QtCore.QSize(0, 14))
        oneSampleLabel.setMaximumSize(QtCore.QSize(50000, 14))
        oneSampleLabel.setAlignment(QtCore.Qt.AlignCenter)
        oneSampleLabel.setObjectName("oneSampleLabel")
        verticalLayout_4.addWidget(oneSampleLabel)
        horizontalLayout_3 = QtGui.QHBoxLayout()
        horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_3.addItem(spacerItem)
        mnoaCheck = QtGui.QCheckBox(frame_9)
        mnoaCheck.setMinimumSize(QtCore.QSize(20, 14))
        mnoaCheck.setMaximumSize(QtCore.QSize(20, 14))
        mnoaCheck.setLayoutDirection(QtCore.Qt.LeftToRight)
        mnoaCheck.setObjectName("mnoaCheck")
        horizontalLayout_3.addWidget(mnoaCheck)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_3.addItem(spacerItem1)
        verticalLayout_4.addLayout(horizontalLayout_3)
        horizontalLayout_7 = QtGui.QHBoxLayout()
        horizontalLayout_7.setObjectName("horizontalLayout_7")
        spacerItem2 = QtGui.QSpacerItem(18, 13, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_7.addItem(spacerItem2)
        mgdCheck = QtGui.QCheckBox(frame_9)
        mgdCheck.setMinimumSize(QtCore.QSize(20, 14))
        mgdCheck.setMaximumSize(QtCore.QSize(20, 14))
        mgdCheck.setObjectName("mgdCheck")
        horizontalLayout_7.addWidget(mgdCheck)
        spacerItem3 = QtGui.QSpacerItem(18, 13, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_7.addItem(spacerItem3)
        verticalLayout_4.addLayout(horizontalLayout_7)
        horizontalLayout_8 = QtGui.QHBoxLayout()
        horizontalLayout_8.setObjectName("horizontalLayout_8")
        spacerItem4 = QtGui.QSpacerItem(18, 13, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_8.addItem(spacerItem4)
        msvCheck = QtGui.QCheckBox(frame_9)
        msvCheck.setMinimumSize(QtCore.QSize(20, 14))
        msvCheck.setMaximumSize(QtCore.QSize(20, 14))
        msvCheck.setObjectName("msvCheck")
        horizontalLayout_8.addWidget(msvCheck)
        spacerItem5 = QtGui.QSpacerItem(18, 13, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_8.addItem(spacerItem5)
        verticalLayout_4.addLayout(horizontalLayout_8)
        horizontalLayout_10 = QtGui.QHBoxLayout()
        horizontalLayout_10.setObjectName("horizontalLayout_10")
        spacerItem6 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_10.addItem(spacerItem6)
        mgwCheck = QtGui.QCheckBox(frame_9)
        mgwCheck.setMinimumSize(QtCore.QSize(20, 14))
        mgwCheck.setMaximumSize(QtCore.QSize(20, 14))
        mgwCheck.setObjectName("mgwCheck")
        horizontalLayout_10.addWidget(mgwCheck)
        spacerItem7 = QtGui.QSpacerItem(18, 13, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_10.addItem(spacerItem7)
        verticalLayout_4.addLayout(horizontalLayout_10)
        self.ui.horizontalLayout_4.addWidget(frame_9)

        self.oneSampleList.append(frame_9)

    def addTwoSample(self,num1,num2):
        """ methode d'ajout d'un bloc dans 'two samples sum stats'
        """
        frame_11 = QtGui.QFrame(self.ui.twoFrame)
        frame_11.setMinimumSize(QtCore.QSize(70, 0))
        frame_11.setMaximumSize(QtCore.QSize(70, 16777215))
        frame_11.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_11.setFrameShadow(QtGui.QFrame.Raised)
        frame_11.setObjectName("frame_11")
        verticalLayout_5 = QtGui.QVBoxLayout(frame_11)
        verticalLayout_5.setSpacing(1)
        verticalLayout_5.setContentsMargins(1, 1, 1, 1)
        verticalLayout_5.setObjectName("verticalLayout_5")
        twoSampleLabel = QtGui.QLabel("Samp %i&%i"%(num1,num2),frame_11)
        font = QtGui.QFont()
        font.setPointSize(8)
        twoSampleLabel.setFont(font)
        twoSampleLabel.setMinimumSize(QtCore.QSize(0, 14))
        twoSampleLabel.setMaximumSize(QtCore.QSize(16777215, 14))
        twoSampleLabel.setObjectName("twoSampleLabel")
        verticalLayout_5.addWidget(twoSampleLabel)
        horizontalLayout_16 = QtGui.QHBoxLayout()
        horizontalLayout_16.setObjectName("horizontalLayout_16")
        spacerItem8 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_16.addItem(spacerItem8)
        mnoa2Check = QtGui.QCheckBox(frame_11)
        mnoa2Check.setMinimumSize(QtCore.QSize(20, 14))
        mnoa2Check.setMaximumSize(QtCore.QSize(20, 14))
        mnoa2Check.setObjectName("mnoa2Check")
        horizontalLayout_16.addWidget(mnoa2Check)
        spacerItem9 = QtGui.QSpacerItem(17, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_16.addItem(spacerItem9)
        verticalLayout_5.addLayout(horizontalLayout_16)
        horizontalLayout_15 = QtGui.QHBoxLayout()
        horizontalLayout_15.setObjectName("horizontalLayout_15")
        spacerItem10 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_15.addItem(spacerItem10)
        mgd2Check = QtGui.QCheckBox(frame_11)
        mgd2Check.setMinimumSize(QtCore.QSize(20, 14))
        mgd2Check.setMaximumSize(QtCore.QSize(20, 14))
        mgd2Check.setObjectName("mgd2Check")
        horizontalLayout_15.addWidget(mgd2Check)
        spacerItem11 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_15.addItem(spacerItem11)
        verticalLayout_5.addLayout(horizontalLayout_15)
        horizontalLayout_14 = QtGui.QHBoxLayout()
        horizontalLayout_14.setObjectName("horizontalLayout_14")
        spacerItem12 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_14.addItem(spacerItem12)
        msv2Check = QtGui.QCheckBox(frame_11)
        msv2Check.setMinimumSize(QtCore.QSize(20, 14))
        msv2Check.setMaximumSize(QtCore.QSize(20, 14))
        msv2Check.setObjectName("msv2Check")
        horizontalLayout_14.addWidget(msv2Check)
        spacerItem13 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_14.addItem(spacerItem13)
        verticalLayout_5.addLayout(horizontalLayout_14)
        horizontalLayout_13 = QtGui.QHBoxLayout()
        horizontalLayout_13.setObjectName("horizontalLayout_13")
        spacerItem14 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_13.addItem(spacerItem14)
        fst2Check = QtGui.QCheckBox(frame_11)
        fst2Check.setMinimumSize(QtCore.QSize(20, 14))
        fst2Check.setMaximumSize(QtCore.QSize(20, 14))
        fst2Check.setObjectName("fst2Check")
        horizontalLayout_13.addWidget(fst2Check)
        spacerItem15 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_13.addItem(spacerItem15)
        verticalLayout_5.addLayout(horizontalLayout_13)
        horizontalLayout_6 = QtGui.QHBoxLayout()
        horizontalLayout_6.setObjectName("horizontalLayout_6")
        ci2LeftCheck = QtGui.QCheckBox(frame_11)
        ci2LeftCheck.setMinimumSize(QtCore.QSize(20, 14))
        ci2LeftCheck.setMaximumSize(QtCore.QSize(20, 14))
        ci2LeftCheck.setObjectName("ci2LeftCheck")
        horizontalLayout_6.addWidget(ci2LeftCheck)
        ci2RightCheck = QtGui.QCheckBox(frame_11)
        ci2RightCheck.setMinimumSize(QtCore.QSize(20, 14))
        ci2RightCheck.setMaximumSize(QtCore.QSize(20, 14))
        ci2RightCheck.setObjectName("ci2RightCheck")
        horizontalLayout_6.addWidget(ci2RightCheck)
        verticalLayout_5.addLayout(horizontalLayout_6)
        horizontalLayout_12 = QtGui.QHBoxLayout()
        horizontalLayout_12.setObjectName("horizontalLayout_12")
        spacerItem16 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_12.addItem(spacerItem16)
        sad2Check = QtGui.QCheckBox(frame_11)
        sad2Check.setMinimumSize(QtCore.QSize(20, 14))
        sad2Check.setMaximumSize(QtCore.QSize(20, 14))
        sad2Check.setObjectName("sad2Check")
        horizontalLayout_12.addWidget(sad2Check)
        spacerItem17 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_12.addItem(spacerItem17)
        verticalLayout_5.addLayout(horizontalLayout_12)
        horizontalLayout_11 = QtGui.QHBoxLayout()
        horizontalLayout_11.setObjectName("horizontalLayout_11")
        spacerItem18 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_11.addItem(spacerItem18)
        dd2Check = QtGui.QCheckBox(frame_11)
        dd2Check.setMinimumSize(QtCore.QSize(20, 14))
        dd2Check.setMaximumSize(QtCore.QSize(20, 14))
        dd2Check.setObjectName("dd2Check")
        horizontalLayout_11.addWidget(dd2Check)
        spacerItem19 = QtGui.QSpacerItem(18, 18, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        horizontalLayout_11.addItem(spacerItem19)
        verticalLayout_5.addLayout(horizontalLayout_11)
        self.ui.horizontalLayout_5.addWidget(frame_11)

        self.twoSampleList.append(frame_11)

    def getStats(self):
        """ retourne les stats sous la forme d'un dico indexé par les noms de ligne
        """
        nstat = 0
        # dico indexé par les noms de ligne (NAL,FST...) qui donne la liste des samples concernés
        dico_stats = {}
        dico_stats["NAL"] = []
        dico_stats["HET"] = []
        dico_stats["VAR"] = []
        dico_stats["MGW"] = []

        dico_stats["N2P"] = []
        dico_stats["H2P"] = []
        dico_stats["V2P"] = []
        dico_stats["FST"] = []
        dico_stats["LIK"] = []
        dico_stats["DAS"] = []
        dico_stats["DM2"] = []

        dico_stats["AML"] = []
        for box in self.oneSampleList:
            lab = str(box.findChild(QLabel,"oneSampleLabel").text()).split()[1]
            if box.findChild(QCheckBox,"mnoaCheck").isChecked():
                dico_stats["NAL"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"mgdCheck").isChecked():
                dico_stats["HET"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"msvCheck").isChecked():
                dico_stats["VAR"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"mgwCheck").isChecked():
                dico_stats["MGW"].append(lab)
                nstat += 1
        for box in self.twoSampleList:
            lab = str(box.findChild(QLabel,"twoSampleLabel").text()).split()[1]
            if box.findChild(QCheckBox,"mnoa2Check").isChecked():
                dico_stats["N2P"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"mgd2Check").isChecked():
                dico_stats["H2P"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"msv2Check").isChecked():
                dico_stats["V2P"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"fst2Check").isChecked():
                dico_stats["FST"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"ci2LeftCheck").isChecked():
                dico_stats["LIK"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"ci2RightCheck").isChecked():
                lab_inverted = "%s&%s"%(lab.split('&')[1],lab.split('&')[0])
                dico_stats["LIK"].append(lab_inverted)
                nstat += 1
            if box.findChild(QCheckBox,"sad2Check").isChecked():
                dico_stats["DAS"].append(lab)
                nstat += 1
            if box.findChild(QCheckBox,"dd2Check").isChecked():
                dico_stats["DM2"].append(lab)
                nstat += 1
        for box in self.admixSampleList:
            lab = str(box.findChild(QLabel,"threeSampleLabel").text()).split('\n')[1]
            if box.findChild(QCheckBox,"ml3Check").isChecked():
                dico_stats["AML"].append(lab)
                nstat +=1
        if dico_stats.has_key("LIK"):
            dico_stats["LIK"] = self.sortLik(dico_stats["LIK"])
        return (nstat,dico_stats)

    def sortLik(self,l):
        res = []
        for i,e in enumerate(l):
            el = int(e.split('&')[0])
            er = int(e.split('&')[1])
            idest = 0
            while idest < len(res) and el > int(res[idest].split('&')[0]):
                #print "elem : %s, plus loin que %s l"%(e,res[idest])
                idest += 1
            while idest < len(res) and er > int(res[idest].split('&')[1]) and el == int(res[idest].split('&')[0]):
                #print "elem : %s, plus loin que %s r"%(e,res[idest])
                idest += 1
            #print "j'insere %s a %s"%(e,idest)
            res.insert(idest,e)
        return res

    def setSumConf(self,lines):
        """ grace aux lignes du fichier de conf, remet les sum stats
        """
        confToStat = {"NAL" : "mnoa","HET" : "mgd", "VAR": "msv","MGW" : "mgw", "N2P" : "mnoa2", "H2P" : "mgd2", "V2P" : "msv2", "FST" : "fst2", "DAS" : "sad2", "DM2" : "dd2", "AML" : "ml3", "LIK" : "" }

        # construction du dico de stats (le même que dans getSumConf)
        dico_stats = {}
        for line in lines:
            sline = line.strip()
            if sline.strip() != "":
                t = sline.split()[0]
                dico_stats[t] = []
                for sample in sline.split()[1:]:
                    dico_stats[t].append(sample)
        #print "dico stats :",dico_stats

        # pour chaque ligne (de case à cocher)
        for k in dico_stats.keys():
            name_chk_box = confToStat[k]
            # pour chaque colonne chaque sample dans one Sample
            if k in ["NAL","HET","VAR","MGW"]:
                for box in self.oneSampleList:
                    num_sample = str(box.findChild(QLabel,"oneSampleLabel").text()).split()[1].strip()
                    if num_sample in dico_stats[k]:
                        # on coche
                        box.findChild(QCheckBox,"%sCheck"%name_chk_box).setChecked(True)
            elif k in ["H2P","V2P","FST","DAS","DM2","N2P"]:
                for box in self.twoSampleList:
                    num_sample = str(box.findChild(QLabel,"twoSampleLabel").text()).split()[1].strip()
                    if num_sample in dico_stats[k]:
                        # on coche
                        box.findChild(QCheckBox,"%sCheck"%name_chk_box).setChecked(True)
            elif k == "LIK":
                for lik in dico_stats[k]:
                    num1 = int(lik.strip().split('&')[0])
                    num2 = int(lik.strip().split('&')[1])
                    ind_box = 0
                    while ind_box < len(self.twoSampleList) and str(self.twoSampleList[ind_box].findChild(QLabel,"twoSampleLabel").text()).strip().split()[1] != lik.strip() and\
                            str(self.twoSampleList[ind_box].findChild(QLabel,"twoSampleLabel").text()).strip().split()[1] != "%s&%s"%(num2,num1):
                        ind_box += 1
                    if num1 > num2:
                        self.twoSampleList[ind_box].findChild(QCheckBox,"ci2RightCheck").setChecked(True)
                    else:
                        self.twoSampleList[ind_box].findChild(QCheckBox,"ci2LeftCheck").setChecked(True)

            elif k == "AML":
                for aml in dico_stats[k]:
                    num1 = aml.strip().split('&')[0]
                    num2 = aml.strip().split('&')[1]
                    num3 = aml.strip().split('&')[2]
                    self.addAdmixSampleGui(int(num1),int(num2),int(num3))
                    self.admixSampleList[-1].findChild(QCheckBox,"%sCheck"%name_chk_box).setChecked(True)



    def clear(self):
        self.parent.clearSummaryStats(self.box_group)

    def getSumStatsTableHeader(self):
        """ retourne l'en-tête pour le tableau à écrire dans conf.th
        """
        result = u""
        if self.numGroup == 0:
            gnumber = self.parent.groupList.index(self.box_group)+1
        else:
            gnumber = self.numGroup
        dico = {"NAL" : [],
        "HET" : [],
        "VAR" : [],
        "MGW" : [],
        "N2P" : [],
        "H2P" : [],
        "V2P" : [],
        "FST" : [],
        "LIK" : [],
        "DAS" : [],
        "DM2" : [],
        "AML" : []}
        for box in self.oneSampleList:
            lab = str(box.findChild(QLabel,"oneSampleLabel").text()).split()[1]
            if box.findChild(QCheckBox,"mnoaCheck").isChecked():
                stat_name = "NAL_%s_%s"%(gnumber,lab)
                dico["NAL"].append(stat_name)
            if box.findChild(QCheckBox,"mgdCheck").isChecked():
                stat_name = "HET_%s_%s"%(gnumber,lab)
                dico["HET"].append(stat_name)
            if box.findChild(QCheckBox,"msvCheck").isChecked():
                stat_name = "VAR_%s_%s"%(gnumber,lab)
                dico["VAR"].append(stat_name)
            if box.findChild(QCheckBox,"mgwCheck").isChecked():
                stat_name = "MGW_%s_%s"%(gnumber,lab)
                dico["MGW"].append(stat_name)
        for box in self.twoSampleList:
            lab = str(box.findChild(QLabel,"twoSampleLabel").text()).split()[1]
            if box.findChild(QCheckBox,"mnoa2Check").isChecked():
                stat_name = "N2P_%s_%s"%(gnumber,lab)
                dico["N2P"].append(stat_name)
            if box.findChild(QCheckBox,"mgd2Check").isChecked():
                stat_name = "H2P_%s_%s"%(gnumber,lab)
                dico["H2P"].append(stat_name)
            if box.findChild(QCheckBox,"msv2Check").isChecked():
                stat_name = "V2P_%s_%s"%(gnumber,lab)
                dico["V2P"].append(stat_name)
            if box.findChild(QCheckBox,"fst2Check").isChecked():
                stat_name = "FST_%s_%s"%(gnumber,lab)
                dico["FST"].append(stat_name)
            if box.findChild(QCheckBox,"ci2LeftCheck").isChecked():
                stat_name = "%s"%(lab)
                dico["LIK"].append(stat_name)
            if box.findChild(QCheckBox,"ci2RightCheck").isChecked():
                lab_inverted = "%s&%s"%(lab.split('&')[1],lab.split('&')[0])
                stat_name = "%s"%(lab_inverted)
                dico["LIK"].append(stat_name)
            if box.findChild(QCheckBox,"sad2Check").isChecked():
                stat_name = "DAS_%s_%s"%(gnumber,lab)
                dico["DAS"].append(stat_name)
            if box.findChild(QCheckBox,"dd2Check").isChecked():
                stat_name = "DM2_%s_%s"%(gnumber,lab)
                dico["DM2"].append(stat_name)
        for box in self.admixSampleList:
            lab = str(box.findChild(QLabel,"threeSampleLabel").text()).split('\n')[1]
            if box.findChild(QCheckBox,"ml3Check").isChecked():
                stat_name = "AML_%s_%s"%(gnumber,lab)
                dico["AML"].append(stat_name)

        dico["LIK"] = self.sortLik(dico["LIK"])
        a=[]
        for e in dico["LIK"]:
            a.append("LIK_%s_%s"%(gnumber,e))
        dico["LIK"] = a

        for t in ["NAL","HET","VAR","MGW","N2P","H2P","V2P","FST","LIK","DAS","DM2","AML"]:
            l = dico[t]
            for s in l:
                result += output.centerHeader(s,14)

        return result

