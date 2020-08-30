import numpy as np
import easyaccess as ea

fileName = "Source_Locations.txt"
names = np.loadtxt("Source_Locations.txt", dtype={'names':('ID', 'RA', 'DEC'), 'formats':(np.int, np.float, np.float)})
tol=0.00027

connection = ea.connect()

outLoc = "query/DESDR_"

for i in range(len(names)):
    print("Getting photometry for "+str(names['ID'][i]))

    q_pt1="select COADD_OBJECT_ID, TILENAME, RA, DEC, MAG_AUTO_G, MAG_AUTO_R,MAG_AUTO_I, MAGERR_AUTO_G, MAGERR_AUTO_R, MAGERR_AUTO_I, "
    q_pt2="FLAGS_G, FLAGS_R, FLAGS_I from DR1_MAIN where ra between "
    q_pt3="dec between " 
    q_and=" and "
    ra=names['RA'][i]
    dec=names['DEC'][i]
    query = q_pt1+q_pt2+str(ra-tol)+q_and+str(ra+tol)+q_and+q_pt3+str(dec-tol)+q_and+str(dec+tol)    

    connection.query_and_save(query,outLoc+str(names['ID'][i])+'.tab')
