# Author: Thomas Toft Lindkvist
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.rc("font", family=["Helvetica", "Arial"]) # skifter skrifttype
plt.rc("axes", labelsize=18)   # skriftstørrelse af `xlabel` og `ylabel`
plt.rc("xtick", labelsize=16, top=True, direction="in")  # skriftstørrelse af ticks og viser ticks øverst
plt.rc("ytick", labelsize=16, right=True, direction="in")
plt.rc("axes", titlesize=18)
plt.rc("legend", fontsize=16)

import numpy as np

import PyQt5.Qt
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout, QSlider, QGraphicsView, QHBoxLayout, QGroupBox, QFormLayout, QSpinBox, QLabel, QCheckBox, QLineEdit
from PyQt5.QtCore import QCoreApplication

from PyQt5.QtGui import QFont, QDoubleValidator

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from scipy.special import mathieu_a, mathieu_b

# Units
amu = 1.6605e-27
ec = 1.60217662e-19

## Setup Constants
RF = 300e3
Omega = 2*np.pi*RF
r0 = 0.0092
z0 = 0.0088

## TODO - Functionality
    # - Load parameters from file
    # - Check stability of high and low mass - Color code for each DC
    # - Plot export image
    # - Export current settings

## TODO - cleanup

lfont = QFont('Arial', 20)
mfont = QFont('Arial', 16)
sfont = QFont('Arial', 12)

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        self.setParent(parent)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):

        super(MainWindow, self).__init__(*args, **kwargs)
        self.setMouseTracking(True)

        self.qs = [np.linspace(0, 1.23978, 100), np.linspace(0.780906, 1.351218, 100), np.linspace(0, 0.780906), np.linspace(1.23978, 1.351218)]
        self.mathieuCurves = [-mathieu_a(0, self.qs[0]), -mathieu_b(1, self.qs[1]), mathieu_a(0, self.qs[2]/2), mathieu_b(1, self.qs[3]/2)]
        
        self.crossings = [0.780906, 1.23978, 1.351218]
        self.RF_q = self.crossings[0]
        self.plot_lims = [[0, 160], [40, -20]]

        # Low, High, Target
        self.cs = ['lightblue', 'lightpink', 'black']
        self.ms = np.array([68*amu, 88*amu, 78*amu])
        self.charge = 1*ec

        self.centralwidget = QWidget()
        self.setCentralWidget(self.centralwidget)
        self.resize(1000, 600)
        self.setWindowTitle("Mathieu Stability Diagram")
        self.create_UI()

        self.updateMass()
        self.show()
    def create_UI(self):
        self.DPI = 100
        self.canvas = MplCanvas(self, width=6, height=4, dpi=self.DPI)

        # RF Slider
        self.RF_slider_scaling = 1/500
        self.RF_slider = QSlider(self)
        self.RF_slider.setMaximum(500)
        self.RF_slider.setMinimum(0)
        self.RF_slider.setValue((self.crossings[0]/self.crossings[1])/self.RF_slider_scaling + 1) # Integer division and round up
        self.RF_slider.setOrientation(PyQt5.QtCore.Qt.Horizontal)
        self.RF_slider.valueChanged.connect(self.eventRFSliderMove)

        # DC Sliders
        self.DC_slider_scaling = 1/200

        self.DC_slider_high = QSlider(self)
        self.DC_slider_high.setMinimum(0)
        self.DC_slider_high.setMaximum(200)
        self.DC_slider_high.setValue(0)
        self.DC_slider_high.setOrientation(PyQt5.QtCore.Qt.Horizontal)
        self.DC_slider_high.valueChanged.connect(self.eventDCSliderMove)

        self.DC_slider_low = QSlider(self)
        self.DC_slider_low.setMinimum(0)
        self.DC_slider_low.setMaximum(200)
        self.DC_slider_low.setValue(0)
        self.DC_slider_low.setOrientation(PyQt5.QtCore.Qt.Horizontal)
        self.DC_slider_low.valueChanged.connect(self.eventDCSliderMove)

        # Check follow graph
        self.check_follow = QCheckBox(self)
        self.check_follow.setCheckState(QtCore.Qt.Checked)
        self.check_follow.stateChanged.connect(self.eventOnCheckChange)

        # Voltage readout
        # self.label_voltage = QLabel(self)
        # self.label_voltage.setText("RF: x\nDC1: y\DC2: z")
        # self.label_voltage.setFont(lfont)
        self.field_RF_voltage = QLineEdit(self)
        self.field_RF_voltage.setValidator(QDoubleValidator(0.0, 2000.0, 1))
        self.field_RF_voltage.setFont(lfont)
        self.field_RF_voltage.setAlignment(QtCore.Qt.AlignRight)
        self.field_RF_voltage.returnPressed.connect(self.eventManualRF)

        self.field_DC1_voltage = QLineEdit(self)
        self.field_DC1_voltage.setValidator(QDoubleValidator(0.0, 2000.0, 1))
        self.field_DC1_voltage.setFont(lfont)
        self.field_DC1_voltage.setAlignment(QtCore.Qt.AlignRight)
        self.field_DC1_voltage.returnPressed.connect(self.eventManualDC(1))

        self.field_DC2_voltage = QLineEdit(self)
        self.field_DC2_voltage.setValidator(QDoubleValidator(0.0, 2000.0, 1))
        self.field_DC2_voltage.setFont(lfont)
        self.field_DC2_voltage.setAlignment(QtCore.Qt.AlignRight)
        self.field_DC2_voltage.returnPressed.connect(self.eventManualDC(2))

        # Mass selection
        self.field_low_mass = QSpinBox(self)
        self.field_low_mass.setMaximum(10000)
        self.field_low_mass.setValue(int(self.ms[0]/amu))
        self.field_high_mass = QSpinBox(self)
        self.field_high_mass.setMaximum(10000)
        self.field_high_mass.setValue(int(self.ms[1]/amu))
        self.field_target_mass = QSpinBox(self)
        self.field_target_mass.setMaximum(10000)
        self.field_target_mass.setValue(int(self.ms[2]/amu))

        # Charge selection
        self.field_charge = QSpinBox(self)
        self.field_charge.setValue(int(self.charge/ec))


        # Update button
        mass_button = QPushButton('Update particles', self)
        mass_button.clicked.connect(self.updateMass)
        mass_button.setFont(sfont)

        # Settings print
        self.label_settings = QLabel(self)
        self.label_settings.setText(f"r_0: {r0*1000:0.1f} mm\nz_0: {z0*1000:0.1f} mm\nRF: {RF/1e3:0.0f} kHz")
        self.label_settings.setFont(sfont)
        
        # Layout setup
        self.hbox = QHBoxLayout(self.centralwidget) 
        self.hbox.addWidget(self.canvas)
        self.horizontalGroupBox = QGroupBox("Settings")
        self.hbox.addWidget(self.horizontalGroupBox)

        self.main_vertical_layout = QFormLayout()
        self.horizontalGroupBox.setLayout(self.main_vertical_layout)
        
        # Settings setup
        self.main_vertical_layout.addRow("High mass", self.field_high_mass)
        self.main_vertical_layout.addRow("Target mass", self.field_target_mass)
        self.main_vertical_layout.addRow("Low mass", self.field_low_mass)
        self.main_vertical_layout.addRow("Charge", self.field_charge)

        self.main_vertical_layout.addRow(mass_button)

        self.main_vertical_layout.addRow("Follow graph", self.check_follow)

        self.main_vertical_layout.addRow("RF V-Tuning", self.RF_slider)
        self.main_vertical_layout.addRow("DC1 V-Tuning", self.DC_slider_high)
        self.main_vertical_layout.addRow("DC2 V-Tuning", self.DC_slider_low)
        self.main_vertical_layout.addRow("RF V: ", self.field_RF_voltage)
        self.main_vertical_layout.addRow("DC 1: ", self.field_DC1_voltage)
        self.main_vertical_layout.addRow("DC 2: ", self.field_DC2_voltage)
        self.main_vertical_layout.addRow(self.label_settings)

        # Set fonts(large) for labels
        for w in (self.field_RF_voltage, self.field_DC1_voltage, self.field_DC2_voltage):
            self.main_vertical_layout.labelForField(w).setFont(lfont)

        # Set fonts(small) for the labels 
        for w in (self.RF_slider, self.DC_slider_high, self.DC_slider_low, self.check_follow):
            self.main_vertical_layout.labelForField(w).setFont(sfont)

        # Set fonts (small) for the input fields and labels
        for w in (self.field_charge, self.field_high_mass, self.field_low_mass, self.field_target_mass):
            self.main_vertical_layout.labelForField(w).setFont(sfont)
            w.setFont(sfont)

    def eventOnCheckChange(self):
        self.eventManualRF()
    def eventRFSliderMove(self):
        
        q_to_V = (self.ms[2]*Omega**2*(r0**2 + 2*z0**2)/(8*self.charge))
        self.RF_q = self.crossings[1]*self.RF_slider.value()*self.RF_slider_scaling
        self.updateDCRange()

        self.updatePlot()
        self.updateReadout()
    def eventDCSliderMove(self):
        self.updatePlot()
        self.updateReadout()
    def eventManualRF(self):
        V_to_q = 8*self.charge/(self.ms[2]*Omega**2*(r0**2 + 2*z0**2))
        self.RF_slider.setValue((float(self.field_RF_voltage.text())*V_to_q/self.crossings[1])/self.RF_slider_scaling + 1)
        self.eventRFSliderMove()
    def eventManualDC(self, n):
        def specific():
            # Convert voltage to percentage
            if(self.check_follow.checkState() == QtCore.Qt.Checked):
                print(self.DC_Range)
                q_to_V = (self.ms[2]*Omega**2*(r0**2 + 2*z0**2)/(8*self.charge))
                
                DC_voltage_high_m = 0
                if(self.RF_q > self.crossings[0]):
                    DC_voltage_high_m = -mathieu_b(1, self.RF_q)*q_to_V/2
                else:
                    DC_voltage_high_m = mathieu_a(0, self.RF_q/2)*q_to_V

                DC_voltage_low_m = -mathieu_a(0, self.RF_q)*q_to_V/2
                DC_middle = (DC_voltage_low_m+DC_voltage_high_m)/2

                if(n == 1):
                    self.DC_slider_high.setValue((1 - abs(float(self.field_DC1_voltage.text()) - DC_middle)/(0.5*self.DC_Range))/self.DC_slider_scaling + 1)
                elif(n == 2):
                    print(DC_middle)
                    self.DC_slider_low.setValue((1 - abs(float(self.field_DC2_voltage.text()) - DC_middle)/(0.5*self.DC_Range))/self.DC_slider_scaling + 1)
            else:
                if(n == 1):
                    self.DC_slider_high.setValue((abs(float(self.field_DC1_voltage.text()))/(0.25*self.DC_Range))/self.DC_slider_scaling + 1)
                elif(n == 2):
                    self.DC_slider_low.setValue((abs(float(self.field_DC2_voltage.text()))/(0.75*self.DC_Range))/self.DC_slider_scaling + 1)
            self.eventDCSliderMove()
        return specific
    def updateDCRange(self):
        q_to_V = (self.ms[2]*Omega**2*(r0**2 + 2*z0**2)/(8*self.charge))
        if(self.check_follow.checkState() == QtCore.Qt.Checked):
            DC_voltage_high_m = 0
            if(self.RF_q > self.crossings[0]):
                DC_voltage_high_m = -mathieu_b(1, self.RF_q)*q_to_V/2
            else:
                DC_voltage_high_m = mathieu_a(0, self.RF_q/2)*q_to_V

            DC_voltage_low_m = -mathieu_a(0, self.RF_q)*q_to_V/2
            self.DC_Range = DC_voltage_low_m-DC_voltage_high_m
        else:
            DC_voltage_high_m = -mathieu_b(1, self.crossings[0])*q_to_V/2
            DC_voltage_low_m = -mathieu_a(0, self.crossings[2])*q_to_V/2
            self.DC_Range = DC_voltage_low_m-DC_voltage_high_m

    def updateMass(self):
        self.ms[0] = self.field_low_mass.value()*amu
        self.ms[1] = self.field_high_mass.value()*amu
        self.ms[2] = self.field_target_mass.value()*amu
        self.charge = self.field_charge.value()*ec
        
        self.updateDCRange()

        self.updatePlot()
        self.updateReadout()
    def updateReadout(self):
        q_to_V = (self.ms[2]*Omega**2*(r0**2 + 2*z0**2)/(8*self.charge))
        # self.RF_q = self.crossings[1]*self.RF_slider.value()*self.RF_slider_scaling
        
        DC_voltage_high_m = 0
        DC_voltage_low_m = 0
        if(self.check_follow.checkState() == QtCore.Qt.Checked):
            if(self.RF_q > self.crossings[0]):
                DC_voltage_high_m = -mathieu_b(1, self.RF_q)*q_to_V/2 + 0.5*self.DC_slider_high.value()*self.DC_slider_scaling*self.DC_Range
            else:
                DC_voltage_high_m = mathieu_a(0, self.RF_q/2)*q_to_V + 0.5*self.DC_slider_high.value()*self.DC_slider_scaling*self.DC_Range

            DC_voltage_low_m = -mathieu_a(0, self.RF_q)*q_to_V/2 - 0.5*self.DC_slider_low.value()*self.DC_slider_scaling*self.DC_Range
        else: 
            DC_voltage_high_m = -0.25*self.DC_slider_high.value()*self.DC_slider_scaling*self.DC_Range
            DC_voltage_low_m = 0.75*self.DC_slider_low.value()*self.DC_slider_scaling*self.DC_Range

        self.field_RF_voltage.setText(f"{self.RF_q*q_to_V:0.2f}")
        self.field_DC1_voltage.setText(f"{DC_voltage_high_m:0.2f}")
        self.field_DC2_voltage.setText(f"{DC_voltage_low_m:0.2f}")

        # TODO: Check stability of high and low


    # TODO: Port updatePort to PyQtGraph
    def updatePlotPyQT(self):
        return 0
    def updatePlot(self):
        self.canvas.axes.cla()
        q_to_Vs = (self.ms*Omega**2*(r0**2 + 2*z0**2)/(8*self.charge))
        #Plotting
        for q_to_V, c in zip(q_to_Vs, self.cs):
            # Z stable solutions between
            self.canvas.axes.plot(self.qs[0]*q_to_V, self.mathieuCurves[0]*q_to_V/2, c=c, lw=1.5)
            self.canvas.axes.plot(self.qs[1]*q_to_V, self.mathieuCurves[1]*q_to_V/2, c=c, lw=1.5)

            # R stable solutions between
            self.canvas.axes.plot(self.qs[2]*q_to_V, self.mathieuCurves[2]*q_to_V, c=c, lw=1.5)
            self.canvas.axes.plot(self.qs[3]*q_to_V, self.mathieuCurves[3]*q_to_V, c=c, lw=1.5)

        # self.RF_q = self.crossings[1]*self.RF_slider.value()*self.RF_slider_scaling
        
        # Plot DC crosses
        if(self.check_follow.checkState() == QtCore.Qt.Checked):
            DC_voltage_high_m = 0
            if(self.RF_q > self.crossings[0]):
                DC_voltage_high_m = -mathieu_b(1, self.RF_q)*q_to_Vs[2]/2 + 0.5*self.DC_slider_high.value()*self.DC_slider_scaling*self.DC_Range
            else:
                DC_voltage_high_m = mathieu_a(0, self.RF_q/2)*q_to_Vs[2] + 0.5*self.DC_slider_high.value()*self.DC_slider_scaling*self.DC_Range
            
            self.canvas.axes.plot(q_to_Vs[2]*(self.RF_q), DC_voltage_high_m, 'kx')
            self.canvas.axes.plot(q_to_Vs[2]*(self.RF_q), -mathieu_a(0, self.RF_q)*q_to_Vs[2]/2 - 0.5*self.DC_slider_low.value()*self.DC_slider_scaling*self.DC_Range, 'kx')
        else:
            self.canvas.axes.plot(q_to_Vs[2]*(self.RF_q), -0.25*self.DC_slider_high.value()*self.DC_slider_scaling*self.DC_Range, 'kx')
            self.canvas.axes.plot(q_to_Vs[2]*(self.RF_q), 0.75*self.DC_slider_low.value()*self.DC_slider_scaling*self.DC_Range, 'kx')

        # Set axes
        self.canvas.axes.set_xlim(self.plot_lims[0])
        self.canvas.axes.set_ylim(self.plot_lims[1])
        self.canvas.axes.set_xlabel("RF Voltage")
        self.canvas.axes.set_ylabel("DC Voltage")
        self.canvas.axes.grid()
        self.canvas.draw()
    
    # Zoom functionality
    def wheelEvent(self, event):
        pos = [event.pos().x(), event.pos().y()]
        canvasSize = self.canvas.get_width_height()
        if(0 <= pos[0] and pos[0] <= canvasSize[0] and 0 <= pos[1] and pos[1] <= canvasSize[1]): # If the mouse is on the canvas, zoom
            x_frac = pos[0]/canvasSize[0]
            y_frac = pos[1]/canvasSize[1]
            
            sign = np.sign(event.angleDelta().y())
            width = self.plot_lims[0][1]-self.plot_lims[0][0]
            height = self.plot_lims[1][0]-self.plot_lims[1][1]
            self.plot_lims[0][0] += 0.1*width*sign*x_frac
            self.plot_lims[0][1] -= 0.1*width*sign*(1-x_frac)
            self.plot_lims[1][0] -= 0.1*height*sign*(1-y_frac)
            self.plot_lims[1][1] += 0.1*height*sign*y_frac
            self.updatePlot()

        return super().wheelEvent(event)

    # On window resize update plot size
    def resizeEvent(self, event):
        size = [event.size().width(), event.size().height()]
        self.canvas.figure.set_size_inches(size[0]*0.7/self.DPI, size[1]/self.DPI)
        self.updatePlot()
        self.hbox.addWidget(self.canvas)
        self.hbox.addWidget(self.horizontalGroupBox)
        return super().resizeEvent(event)

if __name__ == '__main__':
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    else:
        print('QApplication instance already exists: %s' % str(app))
    app.aboutToQuit.connect(app.deleteLater)

    window = MainWindow()

    app.exec_()