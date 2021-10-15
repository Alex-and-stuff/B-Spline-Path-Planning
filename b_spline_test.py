# Code implemented from: https://opensourc.es/blog/b-spline/
# There are some adjustments needed to transfer from Julia to Python
# In line 47, this causes some odd problems where the starting point always starts from the origin (0,0)
import matplotlib.pyplot as plt
import math

def b_spline_basis(i,k,t,u):
    u_i = u[i]
    u_i1 = u[i+1]
    u_ir = u[i+k]
    u_i1r = u[i+k+1]

    if k == 0:
        if u_i <= t and t < u_i1:
            return 1
        else:
            return 0
    else:
        N = b_spline_basis
        left = 0
        if (u_ir-u_i) != 0:
            left = (t-u_i)/(u_ir-u_i) *b_spline_basis(i,k-1,t,u)

        right = 0
        if (u_i1r-u_i1) != 0:
            right = (u_i1r-t)/(u_i1r-u_i1) * b_spline_basis(i+1,k-1,t,u)

        return left + right 

def b_spline(dBx,dBy, k, knot_vector, steps=200):
    m = len(dBx)-1
    if len(knot_vector) != m+1+k+1:
        print("The knot vector has the wrong length it should have m+1+n+1 = ", m+1+k+1, " components.")
        print("But it has ", len(knot_vector), " components")
        return
    u_min = knot_vector[k]
    u_max = knot_vector[-k-1]
    step_size = (u_max-u_min)/(steps-1)
    curve, curve_x, curve_y = [], [], []
    c = 0
    # start = 0
    t = u_min
    while t < u_max:
        pos_x, pos_y = 0,0
        a = 0

        for i in reversed(range(len(knot_vector))):
            if t >= knot_vector[i]:
                start = i+1
                # print(t,knot_vector[i],start)
                break

        # here since we know that it ranges through n+1 sub polynomials every timestep, we could give a limit rather than iterate through the whole m+1 range
        # for i in range(m+1):
        # print(start)
        for i in range(start-k-1, start):
            # if knot_vector[i+1] != knot_vector[i+n+1]: # This line seems like the biggest bug
                pos_x += dBx[i]*b_spline_basis(i,k,t,knot_vector)
                pos_y += dBy[i]*b_spline_basis(i,k,t,knot_vector)

                # print(dBx[i],b_spline_basis(i,k,t,knot_vector))

        curve += [c,[pos_x,pos_y]]
        curve_x.append(pos_x)
        curve_y.append(pos_y)
        c += 1
        # print(t, t+step_size)
        t+= step_size

    # since cannot achieve t = u_max 
    curve_x.append(dBx[-1])
    curve_y.append(dBy[-1])

    plt.scatter(dBx,dBy,c='r')
    plt.plot(dBx,dBy, 'r-')
    plt.scatter(curve_x, curve_y, c='b')
    plt.plot(curve_x, curve_y, c='b')

def generate_knot_vector(dBx, dBy, k): # (Point_x_list, point_y_list, degree_of_spline)
    length_x = len(dBx)
    length_y = len(dBy)
    if length_x != length_y:
        print("path length inconsistent, check if both lengths match")
        return
    length = length_x
    knot_vector = []
    duplacants = k+1
    for i in range(duplacants):
        knot_vector.append(0)
    # inner = -duplacants + length
    inner = length-1-k
    for i in range(inner):
        knot_vector.append(i+1)
    end_duplacants = knot_vector[-1]+1
    for i in range(duplacants):
        knot_vector.append(end_duplacants)
    return knot_vector


if __name__ == '__main__':

    # dBx = [1,-1,5,7, 9,13]
    # dBy = [0, 3,5,4,-2, 1]

    dBx = [0,0,1,3,5,6,6,5,3,1,0,0]
    dBy = [3,1,0,0,1,1,5,6,6,6,5,3]
    degree_of_spline = 3
    knot_vector = generate_knot_vector(dBx,dBy,3)


    # num_of_control_points = len(dBx)
    # n = num_of_control_points - 1
    # k = 2  # degree of curve 
    # knot_vector = []
    # for i in range(n+k+2):
    #     knot_vector.append()
    # [0,0,0,1,2,3,4,5,6,7,7,7]
    # [0,0,0,0,1,2,3,4,5,6,6,6,6]
    # [0,0,0,1,2,3,4,5,6,7,8,8,8]
    # [0,1,2,3,4,5,6,7,8,9,10,11,12]

    # good!!!
    # b_spline(dBx, dBy, 2, [0,0,0,1,2,3,4,4,4])

    # good!!!
    # b_spline(dBx, dBy, 2, [0,0,0,1,2,3,4,5,6,7,7,7])

    # good!!!
    # b_spline(dBx, dBy, 2, [0,1,2,3,4,5,6,7,8,9,10,11])
    
    b_spline(dBx, dBy, degree_of_spline, knot_vector)

    print(generate_knot_vector(dBx,dBy,3))
    plt.show()
