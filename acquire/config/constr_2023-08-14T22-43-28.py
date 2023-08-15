
#Constraint file for experiment: exp_2023-08-14T22-43-28

def constraint(raw_point):
    O8Cl2Ca_mLs, LiF6P_mLs, Cl2Ca_mLs, H6C4O3_mLs = raw_point
    return (LiF6P_mLs==0 and O8Cl2Ca_mLs==0) and (Cl2Ca_mLs+H6C4O3_mLs>=2.0) and (Cl2Ca_mLs+H6C4O3_mLs<=100.0)
