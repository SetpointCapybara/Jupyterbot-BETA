from utils import *
import numpy as np
import robot
from box import *
from cylinder import *
from meshmaterial import *
from model3d import *
from links import *


def _create_kukakr5(htm, name, color, opacity):

    if not Utils.is_a_matrix(htm, 4, 4):
        raise Exception("The optional parameter 'htm' should be a 4x4 homogeneous transformation matrix")

    if not (str(type(name)) == "<class 'str'>"):
        raise Exception("The optional parameter 'name' should be a string")

    if not Utils.is_a_color(color):
        raise Exception("The parameter 'color' should be a color")

    if (not Utils.is_a_number(opacity)) or opacity < 0 or opacity > 1:
        raise Exception("The parameter 'opacity' should be a float between 0 and 1")


    link_info = [[0,0,0,0,0,0],
                [0.335, 0.000, 0.000, -0.405, 0.000, 0.080],  # "d" translation in z
                [-1.570, 0.000, 1.570, -1.570, -1.570, 0],  # "alfa" rotation in x
                [0.075, 0.365, 0.090, 0.000, 0.000, 0],  # "a" translation in x
                [0.000, 0.000, 0.000, 0.000, 0.000, 0.000]];

    n = 6

    # Collision model
    col_model = [[], [], [], [], [], []]

    htm_0_0 = np.array([[7.963e-04, 1.000e+00, 1.067e-22, -7.500e-02],
                        [-7.963e-04, 6.341e-07, -1.000e+00, 1.700e-01],
                        [-1.000e+00, 7.963e-04, 7.963e-04, -1.354e-04],
                        [0.000e+00, 0.000e+00, 0.000e+00, 1.000e+00]])

    htm_0_1 = np.array([[7.671e-08, 1.000e+00, 7.963e-04, -2.391e-05],
                        [1.000e+00, 6.341e-07, -8.927e-04, -4.976e-03],
                        [-8.927e-04, 7.963e-04, -1.000e+00, 3.006e-02],
                        [0.000e+00, 0.000e+00, 0.000e+00, 1.000e+00]])

    htm_1_0 = np.array([[7.970e-04, 7.957e-04, 1.000e+00, -2.001e-01],
                        [7.957e-04, 1.000e+00, -7.963e-04, 1.977e-02],
                        [-1.000e+00, 7.963e-04, 7.963e-04, 1.202e-01],
                        [0.000e+00, 0.000e+00, 0.000e+00, 1.000e+00]])

    htm_1_1 = np.array([[-1.000e+00, 7.957e-04, 8.933e-04, 9.968e-03],
                        [7.964e-04, 1.000e+00, 7.956e-04, -3.305e-04],
                        [-8.927e-04, 7.963e-04, -1.000e+00, 4.036e-02],
                        [0.000e+00, 0.000e+00, 0.000e+00, 1.000e+00]])

    htm_2_0 = np.array([[7.970e-04, 7.957e-04, 1.000e+00, 1.787e-04],
                        [-1.000e+00, 1.593e-03, 7.957e-04, 7.801e-04],
                        [-1.592e-03, -1.000e+00, 7.970e-04, -2.246e-01],
                        [0.000e+00, 0.000e+00, 0.000e+00, 1.000e+00]])

    col_model[0].append(Cylinder(htm=htm_0_0, name=name + "_C0_0", radius=0.12, height=0.33, color="red"))
    col_model[0].append(Cylinder(htm=htm_0_1, name=name + "_C0_1", radius=0.095, height=0.30, color="red"))

    col_model[1].append(Box(htm=htm_1_0, name=name + "_C1_0", width=0.1, height=0.5, depth=0.16, color="green"))
    col_model[1].append(Cylinder(htm=htm_1_1, name=name + "_C1_1", radius=0.095, height=0.28, color="green"))

    col_model[2].append(Box(htm=htm_2_0, name=name + "_C2_0", width=0.143, height=0.12, depth=0.45, color="blue"))


    #Create 3d objects


    base_3d_obj = Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Base.obj', 0.001,
                          Utils.rotz(3.14) @ Utils.rotx(3.14/2),
                          MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color="#242526", opacity=opacity))

    link_3d_obj = []

    link_3d_obj.append(
        Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Axis1.obj', 0.001,
                Utils.trn([0,0,0.203]) @ Utils.rotx(3.14/2),
                MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color=color, opacity=opacity))
    )

    link_3d_obj.append(
        Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Axis2.obj', 0.001,
                Utils.trn([0,0,0.1]) @ Utils.rotz(-3.14 / 2+3.14/13)  @ Utils.rotx(3.14/2),
                MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color=color, opacity=opacity))
    )

    link_3d_obj.append(
        Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Axis3.obj', 0.001,
                Utils.rotx(3.14/2) @ Utils.rotz(-3.14/2),
                MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color=color, opacity=opacity))
    )

    link_3d_obj.append(
        Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Axis4.obj', 0.001,
                Utils.trn([0,0,-0.218]) @ Utils.rotx(3.14/2) @ Utils.rotz(-3.14/2),
                MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color=color, opacity=opacity))
    )

    link_3d_obj.append(
        Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Axis5.obj', 0.001,
                Utils.rotx(3.14/2),
                MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color=color, opacity=opacity))
    )

    link_3d_obj.append(
        Model3D('https://raw.githubusercontent.com/SetpointCapybara/kukakr5/main/models/Axis6.obj', 0.001,
                Utils.trn([0,0,-0.012]) @ Utils.rotx(3.14/2),
                MeshMaterial(metalness=0.7, clearcoat=1, roughness=0.5, normal_scale=[0.5, 0.5], color="#242526", opacity=opacity))
    )

    #Create links

    links = []
    for i in range(n):
        links.append(Link(i, link_info[0][i], link_info[1][i], link_info[2][i], link_info[3][i], link_info[4][i], link_3d_obj[i]))

        for j in range(len(col_model[i])):
            links[i].attach_col_object(col_model[i][j], col_model[i][j].htm)

    #Define initial configuration
    q0 = [1.570, -1.570,  0.000,  0.000,  0,  0.000]

    return base_3d_obj, links, q0
