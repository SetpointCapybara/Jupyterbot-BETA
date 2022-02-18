from simulation import *
import numpy as np
from box import *
from utils import *
from meshmaterial import *
from pointlight import *
import robot as rb
import sys



def _control_demo_1():
    # Create simulation and add objects to the scene
    robot = rb.Robot.create_kukakr5(np.identity(4), "robo")
    mesh_board = MeshMaterial(roughness=1, metalness=0.9)
    board = Box(htm=Utils.trn([0.6, 0, 0.5]), width=0.05, depth=0.9, height=0.8, color="white",
                 mesh_material=mesh_board)
    mesh_ground = MeshMaterial(texture_map='https://i.imgur.com/oduDrf3.png', roughness=1, metalness=1)
    ground = Box(name="ground", width=3, depth=3, height=0.01, htm=Utils.trn([0, 0, 0.005]), mesh_material=mesh_ground)

    light1 = PointLight(name="light1", color="white", intensity=3, position=[1, 0, 1.5])
    light2 = PointLight(name="light2", color="white", intensity=3, position=[-1, 0, 1.5])
    light3 = PointLight(name="light3", color="white", intensity=3, position=[0, 0, 1.5])
    light4 = PointLight(name="light4", color="white", intensity=3, position=[0, 0.5, 0.5])

    sim = Simulation([robot, board, ground, light1, light2, light3, light4], ambient_light_intensity=3)



    # Create curve
    theta = np.linspace(0, 2 * np.pi, num=300)
    curve = []
    for t in theta:
        curve.append([0.575, 0.2 * cos(t), 0.2 * sin(t) + 0.5])

    # Create vector field
    vf = rb.Robot.vector_field(curve, 10, 0.3)

    # Parameters
    dt = 0.01
    time_max = 20
    K = 1
    imax = round(time_max / dt)

    # Initializations
    hist_time = []
    hist_qdot = []
    hist_q = []
    hist_error_ori = []

    x_des = np.array([0, 0, 1]).reshape((3, 1))
    y_des = np.array([0, -1, 0]).reshape((3, 1))
    z_des = np.array([1, 0, 0]).reshape((3, 1))

    # Main loop
    for i in range(imax):

        if i % 50 == 0 or i == imax - 1:
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * round(20 * i / (imax - 1)), round(100 * i / (imax - 1))))
            sys.stdout.flush()

        jac_eef, htm_eef = robot.jac_geo()

        p_eef = htm_eef[0:3, 3]
        x_eef = htm_eef[0:3, 0]
        y_eef = htm_eef[0:3, 1]
        z_eef = htm_eef[0:3, 2]

        target = np.zeros((6,))
        target[0:3] = vf(p_eef)
        target[3] = -K * sqrt(max(1 - np.transpose(x_des) @ x_eef, 0))
        target[4] = -K * sqrt(max(1 - np.transpose(y_des) @ y_eef, 0))
        target[5] = -K * sqrt(max(1 - np.transpose(z_des) @ z_eef, 0))

        jac_target = np.zeros((6, 6))
        jac_target[0:3, :] = jac_eef[0:3, :]
        jac_target[3, :] = np.transpose(x_des) @ Utils.S(x_eef) @ jac_eef[3:6, :]
        jac_target[4, :] = np.transpose(y_des) @ Utils.S(y_eef) @ jac_eef[3:6, :]
        jac_target[5, :] = np.transpose(z_des) @ Utils.S(z_eef) @ jac_eef[3:6, :]

        qdot = Utils.dp_inv(jac_target, 0.002) @ target
        q_prox = np.array(robot.q).reshape((6,)) + qdot * dt

        robot.add_ani_frame(i * dt, q_prox)

        hist_time.append((i - 1) * dt)
        hist_q.append(robot.q)
        hist_error_ori.append([(180 / (np.pi)) * acos(1 - num * num / (K * K)) for num in target[3:6]])
        hist_qdot.append(qdot)

    # Run simulation
    sim.run()

    # Plot graphs
    Utils.plot(hist_time, np.transpose(hist_q), "", "Time (s)", "Joint configuration  (rad)", "q")
    Utils.plot(hist_time, np.transpose(hist_qdot), "", "Time (s)", "Joint speed (rad/s)", "u")
    Utils.plot(hist_time, np.transpose(hist_error_ori), "", "Time (s)", "Relative orientation error (degrees)", ['x', 'y', 'z'])

    return sim
