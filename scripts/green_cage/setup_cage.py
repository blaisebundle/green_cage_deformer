from functools import partial

import pymel.core as pm
import maya.OpenMayaUI as mui

from PySide2 import QtWidgets, QtCore, QtGui
import shiboken2


def get_maya_window():
    ptr = mui.MQtUtil.mainWindow()
    if ptr:
        return shiboken2.wrapInstance(long(ptr), QtWidgets.QWidget)
    else:
        return None


def apply_cage(target_mesh=None, cage_mesh=None):
    if not target_mesh or not cage_mesh:
        currently_selected = pm.selected()    
        if len(currently_selected) != 2:        
            warning('Select a target and cage mesh!')
            return None
        target_mesh, cage_mesh = currently_selected
    else:
        target_mesh = pm.PyNode(target_mesh)
        cage_mesh = pm.PyNode(cage_mesh)
        
    if not pm.pluginInfo('green_cage_deformer', loaded=True, query=True):
        pm.loadPlugin('green_cage_deformer')

    pm.select(target_mesh, r=True)
    cage_deformer = pm.deformer(type='green_cage_deformer')[0]
    cage_mesh.worldMesh[0] >> cage_deformer.cageMesh

    return cage_deformer


class SelectMeshWidget(QtWidgets.QWidget):
    def __init__(self, name='', parent=None):
        super(SelectMeshWidget, self).__init__(parent)
        self.name = name
        self._setup_ui()
        self._setup_signals()

    def _setup_ui(self):
        main_layout = QtWidgets.QHBoxLayout()

        self.label = QtWidgets.QLabel(self.name)
        self.label.setMinimumWidth(80)
        self.mesh_le = QtWidgets.QLineEdit()
        self.mesh_le.setReadOnly(True)
        self.select_btn = QtWidgets.QPushButton('<<')

        main_layout.addWidget(self.label)
        main_layout.addWidget(self.mesh_le)
        main_layout.addWidget(self.select_btn)

        self.setLayout(main_layout)

    def _setup_signals(self):
        self.select_btn.clicked.connect(self.add_to_mesh_le)

    def add_to_mesh_le(self):
        currently_selected = pm.selected()
        self.mesh_le.setText(currently_selected[0].name())

    def get_name(self):
        return str(self.mesh_le.text())


class GCDeformerUI(QtWidgets.QDialog):
    _widgets = []
    def __init__(self, parent=None):
        super(GCDeformerUI, self).__init__(parent)
        self._setup_ui()
        self._setup_signals()
        self._widgets.append(self)

    def _setup_ui(self):
        self.setWindowTitle('Green Cage Deformer UI')
        main_layout = QtWidgets.QVBoxLayout()
        self.target_widget = SelectMeshWidget(name='Target Mesh')
        self.cage_widget = SelectMeshWidget(name='Cage Mesh')
        self.setup_cage_btn = QtWidgets.QPushButton('Setup Cage')
        
        main_layout.addWidget(self.target_widget)
        main_layout.addWidget(self.cage_widget)
        main_layout.addWidget(self.setup_cage_btn)
        self.setLayout(main_layout)

    def _setup_signals(self):
        self.setup_cage_btn.clicked.connect(self.setup_cage)

    def setup_cage(self):
        target_mesh = self.target_widget.get_name()
        cage_mesh = self.cage_widget.get_name()
        apply_cage(target_mesh=target_mesh, cage_mesh=cage_mesh)

    @classmethod
    def show_window(cls):
        if cls._widgets:
            cls._widgets[0].close()
            cls._widgets[0].deleteLater()
            cls._widgets.pop()


        ui = cls(parent=get_maya_window())
        ui.show()
