import numpy as np
import matplotlib.pyplot as plt
import math

def  stangebyModel(r, fd, potential, ti, te, ne, B, background_amu):
    # Stangeby model for sheath potential
    me = 9.10938356e-31
    sheath_fac = abs(0.5*np.log((2*np.pi*me/background_amu)*(1+ti/te)))

    norm = np.acos(np.pow(math.exp(1), -sheath_fac))
    fd = 1.0 + np.log(np.cos(angle / 90.0 * norm)) / sheath_fac

    potential = sheath_fac*te 

    debyeLength = np.sqrt(k * te / (ne * e * e))
    larmorRadius = np.sqrt(2.0 * k * ti / (background_amu * me)) / (e * B)

    # Standby model for sheath potential
    potential_DS = fd * np.exp(-r / (2.0 * debyeLength))
    potential_CD = (1.0 - fd) * np.exp(-r / larmorRadius)
    return potential * (potential_DS + potential_CD)



def  stangebyModel(r, fd, potential,debyeLength, larmorRadius):
    # Stangeby model for sheath potential

    # Standby model for sheath potential
    potential_DS = fd * np.exp(-r / (2.0 * debyeLength))
    potential_CD = (1.0 - fd) * np.exp(-r / larmorRadius)
    return -potential * (potential_DS + potential_CD)

    # # Standby model for sheath potential
    # E_component_DS = fd / (2.0 * debyeLength) * np.exp(-r / (2.0 * debyeLength))
    # E_component_CD = (1.0 - fd) / larmorRadius * np.exp(-r / larmorRadius)
    # Emag = potential * (E_component_DS + E_component_CD)
    # return abs(Emag)

def  couletteManfrediModel(r, alpha):
    # Coulette Manfredi model for sheath potential
    C1_alpha = -0.00281407 * alpha - 2.31655435
    C2_alpha = 0.00640402 * alpha + 0.01023915
    # k = 1.38064852e-23
    # e = 1.60217662e-19
    # kte = k * te * 11605.0/ e
    # debyeLength=1
    # result = C1_alpha * C2_alpha *  np.exp(-C2_alpha * r) * kte  
    result = C1_alpha *  np.exp(-C2_alpha * r )   
    # derivative = abs(-C1_alpha * C2_alpha * np.exp(-C2_alpha * r /DebyeLength) * kte) /DebyeLength
    return result


# minium distance 0.012990 

k = 1.38064852e-23
e = 1.60217662e-19
te  = 20.000000
DebyeLength= 0.000011 
kte = k * te * 11605.0/ e
pot = 57.630254 / kte
fd = 0.780993 
larmorRadius = 0.000322
# angle = 180/2-59.999966 
r =  0.012990 
# print("Stangeby: ", stangebyModel(r, fd, pot, DebyeLength, larmorRadius))
# print("Coulette: ", couletteManfrediModel(r/DebyeLength, te, DebyeLength, angle))
ne = 1.0E+19
_angle =59.999966 # radians
angle0 = _angle  #* 180 / np.pi # degrees

# for angle in [1.8,5, 10, 20, 40, 60, 87]:
for angle in [2,20,angle0]:
    r = np.linspace(0,200,1000)
    print("angle: ", angle)
    # print("Stangeby: ", stangebyModel(r, fd, pot, DebyeLength, larmorRadius))
    # print("Coulette: ", couletteManfrediModel(r/DebyeLength, te, DebyeLength, angle))

    # plt.plot(r*1e3, couletteManfrediModel(r/DebyeLength, te, DebyeLength, angle)[1], label='angle = '+str(angle))

    # plt.plot(r*1e3, stangebyModel(r, fd, pot, DebyeLength, larmorRadius), label='Stangeby')

    # plt.plot(r*1e3, couletteManfrediModel(r/DebyeLength, te, DebyeLength, angle)[0], label='angle = '+str(angle))

    # plt.plot(r*1e3, stangebyModel(r, fd, pot, DebyeLength, larmorRadius), label='Stangeby')
    plt.plot(r, couletteManfrediModel(r, angle), label='PIC')
    plt.plot(r, stangebyModel(r*DebyeLength, fd, pot, DebyeLength, larmorRadius), label='Stangeby')

    # plt.ylim(0, 1e5) 
    # plt.xlim(0.,10)
    plt.legend()
    plt.xlabel('r [mm]' )
    plt.ylabel('Potential')
plt.title('Coulette Manfredi Model')
plt.grid(True)
plt.show()