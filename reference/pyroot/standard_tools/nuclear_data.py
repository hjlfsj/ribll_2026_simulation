





def get_mass_data():
    '''
    Docstring for get_mass_data
    get the atomic mass data from nndc
    '''
# 931.49410242 MeV/c^2
# Mass Excess经过校对，与nndc数据一致，第一项不一定准确
    ATOMIC_MASS_DB = {
    '1H': (938.783, 7.288971),
    '2H': (1875.61294257, 13.135722),
    '3H': (2808.92113298, 14.949810),
    '4H': (3728.4, 24.6),
    
    '3He': (2808.391607, 14.931213),
    '4He': (3727.3794066, 2.424915),
    '5He': (4668.0, 11.231),
    '6He': (5605.577, 17.592),
    '8He': (7472.0, 31.61),
    
    '6Li': (5601.518, 14.086),
    '7Li': (6533.834, 14.908),
    '8Li': (7471.636, 20.946),
    '9Li': (8401.762, 24.954),
    '11Li': (10270.0, 40.728),
    
    '7Be': (6534.184, 15.769),
    '9Be': (8392.748, 11.348),
    '10Be': (9326.104, 12.607),
    '11Be': (10263.0, 20.17),
    
    '10B': (9324.436, 12.051),
    '11B': (10252.542, 8.668),
    '12B': (11178.5, 13.37),

    '10C': (0, 15.699),
    '11C': (0, 10.649),
    '12C': (11174.862, 0.0),
    '13C': (12109.482, 3.125),
    '14C': (13040.988, 3.020),
    '15C': (0, 9.873),
    
    '13N': (12110.6, 5.345),
    '14N': (13043.2, 2.863),
    '15N': (13971.8, 0.101),
    
    '14O': (13036.0, 8.008),
    '15O': (13971.2, 2.855),
    '16O': (14895.079, -4.737),
    '17O': (15828.0, -0.809),
    '18O': (16763.0, -0.782),
    }
    return ATOMIC_MASS_DB



def get_particle_mass(particle):
    '''
    get the mass of a particle
    iuput formation : "6Li"
    return : mass in MeV/c^2
    '''
    ATOMIC_MASS_DB = get_mass_data()
    try:
        import re
        match = re.match(r'(\d+)([A-Za-z]+)', particle)
        if match:
            isotope = match.group(1) + match.group(2)
            if isotope in ATOMIC_MASS_DB:
                return ATOMIC_MASS_DB[isotope][1] + int(match.group(1)) * 931.494
            else:
                A = int(match.group(1))
                u = 931.494
                return A * u
    except:
        pass
    
    raise ValueError(f"未知粒子符号: {particle}")





if __name__ == '__main__':
    print(get_particle_mass('6Li'))
    print(get_particle_mass('16C')/ 16)
    print(get_particle_mass('N'))
