from simulation import *
import numpy as np
from box import *
from utils import *
import robot as rb
from scipy.linalg import null_space
from meshmaterial import *
from pointlight import *
import sys


def _control_demo_2():
    # Create simulation and add objects to the scene

    mesh_ground = MeshMaterial(texture_map='https://i.imgur.com/oduDrf3.png', roughness=1, metalness=1)
    mesh_box = MeshMaterial(texture_map='https://i.imgur.com/oduDrf3.png', roughness=1, metalness=1)
    mesh_plate = MeshMaterial(metalness=0.9, roughness=0.05, env_map_intensity=0.9, clearcoat=1, opacity=1, \
                            reflectivity=0.2, refraction_ratio=0.985, ior=1.52, specular_intensity=0.1,
                              specular_color="white", transmission=1, side="BackSide")

    robot_a = rb.Robot.create_kukakr5(name="robot_a", htm=Utils.trn([0, -0.9, 0.3]), color="#df6c25")
    box_a = Box(name="box_a", width=0.3, depth=0.3, height=0.3, htm=Utils.trn([0, -0.9, 0.15]), mesh_material=mesh_box)

    robot_b = rb.Robot.create_kukakr5(name="robot_b", htm=Utils.rotz(np.pi) @ Utils.trn([0, -0.9, 0.3]), color="#4d6fc4")
    box_b = Box(name="box_b", width=0.3, depth=0.3, height=0.3, htm=Utils.trn([0, 0.9, 0.15]), mesh_material=mesh_box)

    box_left = Box(name="box_left", width=0.3, depth=0.3, height=0.366, htm=Utils.trn([0.3, 0, 0.183]), mesh_material=mesh_box)
    box_right = Box(name="box_right", width=0.3, depth=0.3, height=0.366, htm=Utils.trn([-0.3, 0, 0.183]), mesh_material=mesh_box)

    ground = Box(name="ground", width=3, depth=3, height=0.01, htm=Utils.trn([0, 0, 0.005]), mesh_material=mesh_ground)
    plate = Box(name="plate", width=0.2, depth=0.489, height=0.05, htm=Utils.trn([-0.3, 0, 0.391]), mesh_material=mesh_plate)

    light1 = PointLight(name="light1", color="white", intensity=3, position=[1, 0, 1.5])
    light2 = PointLight(name="light2", color="white", intensity=3, position=[-1, 0, 1.5])
    light3 = PointLight(name="light3", color="white", intensity=3, position=[0, 0, 1.5])
    light4 = PointLight(name="light4", color="white", intensity=3, position=[0, 0.5, 0.5])

    sim = Simulation([robot_a, robot_b, box_a, box_b, box_left, box_right, ground, plate, light1, light2, light3, light4], ambient_light_intensity=3)


    # Use inverse kinematics to set the robot to the starting pose
    q_a = robot_a.ikm(Utils.trn([-0.3, 0.1, -0.7]) @ robot_a.fkm())
    q_b = robot_b.ikm(Utils.trn([-0.3, -0.1, -0.7]) @ robot_b.fkm())

    robot_a.add_ani_frame(0, q_a)
    robot_b.add_ani_frame(0, q_b)

    mth_a_b_des = np.linalg.inv(robot_a.fkm()) @ robot_b.fkm()
    mth_plate = np.linalg.inv(robot_a.fkm()) @ plate.htm

    # Simulation parameters
    dt = 0.005
    time_max_1 = 3 #1
    time_max_2 = 4 #1.5

    jmax = round((2 * time_max_1 + time_max_2) / dt)
    K = 1

    # Initializations
    hist_time = []
    hist_qdot_a = []
    hist_q_a = []
    hist_qdot_b = []
    hist_q_b = []
    hist_error_rel_pos = []
    hist_error_rel_ori = []
    hist_time.append(0)

    # Main loop
    # First part: lift plate
    mth_a_des = Utils.trn([0, 0, 0.1]) @ robot_a.fkm()
    j = 0
    for i in range(round(time_max_1 / dt)):

        if j % 50 == 0 or j == jmax - 1:
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * round(20 * j / (jmax - 1)), round(100 * j / (jmax - 1))))
            sys.stdout.flush()

        r, jac_r = rb.Robot.coop_task_function(robot_a, robot_b, mth_a_des, mth_a_b_des, q_a, q_b)

        r[3:6] = np.sqrt(r[3:6])
        r[9:12] = np.sqrt(r[9:12])

        u1 = Utils.dp_inv(jac_r[0:6, :], 0.001) @ (-K * r[0:6])
        null_jac_r = null_space(jac_r[0:6, :])

        u2 = Utils.dp_inv(jac_r[6:12, :] @ null_jac_r, 0.001) @ (-K * r[6:12] - jac_r[6:12, :] @ u1)

        u = u1 + null_jac_r @ u2
        q_a += u[0:6] * dt
        q_b += u[6:12] * dt

        robot_a.add_ani_frame(j * dt, q_a)
        robot_b.add_ani_frame(j * dt, q_b)
        plate.add_ani_frame(j * dt, robot_a.fkm() @ mth_plate)

        hist_time.append(hist_time[-1] + dt)
        hist_q_a.append(robot_a.q)
        hist_qdot_a.append(u[0:6])
        hist_q_b.append(robot_b.q)
        hist_qdot_b.append(u[6:12])
        hist_error_rel_pos.append(r[0:3])
        hist_error_rel_ori.append([(180 / (np.pi)) * acos(1 - num * num) for num in r[3:6]])

        j += 1

        # Second part: move to the left (negative x axis)
    mth_a_des = Utils.trn([0.6, 0, 0]) @ robot_a.fkm()
    for i in range(round(time_max_2 / dt)):

        if j % 50 == 0 or j == jmax - 1:
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * round(20 * j / (jmax - 1)), round(100 * j / (jmax - 1))))
            sys.stdout.flush()

        r, jac_r = rb.Robot.coop_task_function(robot_a, robot_b, mth_a_des, mth_a_b_des, q_a, q_b)

        r[3:6] = np.sqrt(r[3:6])
        r[9:12] = np.sqrt(r[9:12])

        u1 = Utils.dp_inv(jac_r[0:6, :], 0.001) @ (-K * r[0:6])
        null_jac_r = null_space(jac_r[0:6, :])

        u2 = Utils.dp_inv(jac_r[6:12, :] @ null_jac_r, 0.001) @ (-K * r[6:12] - jac_r[6:12, :] @ u1)

        u = u1 + null_jac_r @ u2
        q_a += u[0:6] * dt
        q_b += u[6:12] * dt

        robot_a.add_ani_frame(j * dt, q_a)
        robot_b.add_ani_frame(j * dt, q_b)
        plate.add_ani_frame(j * dt, robot_a.fkm() @ mth_plate)

        hist_time.append(hist_time[-1] + dt)
        hist_q_a.append(robot_a.q)
        hist_qdot_a.append(u[0:6])
        hist_q_b.append(robot_b.q)
        hist_qdot_b.append(u[6:12])
        hist_error_rel_pos.append(r[0:3])
        hist_error_rel_ori.append([(180 / (np.pi)) * acos(1 - num * num) for num in r[3:6]])

        j += 1

        # Third part: move down the plate
    mth_a_des = Utils.trn([0, 0, -0.1]) @ robot_a.fkm()
    for i in range(round(time_max_1 / dt)):

        if j % 50 == 0 or j == jmax - 1:
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * round(20 * j / (jmax - 1)), round(100 * j / (jmax - 1))))
            sys.stdout.flush()

        r, jac_r = rb.Robot.coop_task_function(robot_a, robot_b, mth_a_des, mth_a_b_des, q_a, q_b)

        r[3:6] = np.sqrt(r[3:6])
        r[9:12] = np.sqrt(r[9:12])

        u1 = Utils.dp_inv(jac_r[0:6, :], 0.001) @ (-K * r[0:6])
        null_jac_r = null_space(jac_r[0:6, :])

        u2 = Utils.dp_inv(jac_r[6:12, :] @ null_jac_r, 0.001) @ (-K * r[6:12] - jac_r[6:12, :] @ u1)

        u = u1 + null_jac_r @ u2
        q_a += u[0:6] * dt
        q_b += u[6:12] * dt

        robot_a.add_ani_frame(j * dt, q_a)
        robot_b.add_ani_frame(j * dt, q_b)
        plate.add_ani_frame(j * dt, robot_a.fkm() @ mth_plate)

        hist_time.append(hist_time[-1] + dt)
        hist_q_a.append(robot_a.q)
        hist_qdot_a.append(u[0:6])
        hist_q_b.append(robot_b.q)
        hist_qdot_b.append(u[6:12])
        hist_error_rel_pos.append(r[0:3])
        hist_error_rel_ori.append([(180 / (np.pi)) * acos(1 - num * num) for num in r[3:6]])

        j += 1

        # Run simulation
    sim.run()

    # Plot graphs
    hist_time.pop()
    Utils.plot(hist_time, np.transpose(hist_q_a), "", "Time (s)", "Joint configuration robot A (rad)", "q")
    Utils.plot(hist_time, np.transpose(hist_qdot_a), "", "Time (s)", "Joint speed robot A (rad/s)", "u")
    Utils.plot(hist_time, np.transpose(hist_q_b), "", "Time (s)", "Joint configuration robot B (rad)", "q")
    Utils.plot(hist_time, np.transpose(hist_qdot_b), "", "Time (s)", "Joint speed robot B (rad/s)", "u")
    Utils.plot(hist_time, np.transpose(hist_error_rel_pos), "", "Time (s)", "Relative position error (m)", ['x', 'y', 'z'])
    Utils.plot(hist_time, np.transpose(hist_error_rel_ori), "", "Time (s)", "Relative orientation error (degrees)", ['x', 'y', 'z'])

    return sim
