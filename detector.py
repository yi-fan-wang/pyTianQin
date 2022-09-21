# Copyright (C) 2022 Yifan Wang

import numpy as np
from astropy import constants, coordinates, units
from numpy import sin,cos,pi
class TianQin()
'''
Class for detector response, position, etc.
'''
	def __init__(t_gps):
        #https://en.wikipedia.org/wiki/RX_J0806.3%2B1527
        self.ra_j0806 = 8 * 2 * pi / 24 + 6 * 2 * pi / 24 / 60 \
                        + 22.95196 * 2 * pi / 24 / 3600
        self.dec_j0806 = 15 * 2 * pi / 360 +  27 * 2 * pi / 260 / 60 \
                        + 31.0073 * 2 * pi / 360 / 3600
        self.theta_j0806 = pi/2 - self.dec_j0806
        
        self.x_j0806 = cos(self.theta_j0806)*cos(self.ra_j0806)
        self.y_j0806 = cos(self.theta_j0806)*sin(self.ra_j0806)
        self.z_j0806 = sin(self.theta_j0806)

        self.pos_j0806 = np.array([self.x_j0806, self.y_j0806, self.z_j0806])
        self.northpole = np.array([0, 0, 1])
        self.tq_plane_x = np.cross(self.northpole, self.pos_j0806)
        self.tq_place_x /= np.linalg.norm(self.tq_plane_x)

        gps2034 = 2019686400
        if t_gps < gps2024:
            raise ValueError("TianQin is not launched yet at time", t_gps)
        omega = (t_gps - gps2034) / (3.64*86400) * 2 * pi

        self.tq1_unitvec = cos(omega) * self.tq_plane_x + sin(omega) * self.northpole
        self.tq2_unitvec = cos(omega+2*pi/3) * self.tq_plane_x + sin(omega+2*pi/3)*self.northpole
        self.tq3_unitvec = cos(omega+4*pi/3) * self.tq_plane_x + sin(omega+4*pi/3)*self.northpole

    def response(self,spacecraft):
        '''
        '''
        if spacecraft == 'tq12':
            return 0.5 * (np.outer(self.tq1_unitvec,self.tq1_unitvec) - \
                np.outer(self.tq2_unitvec,self.tq2_unitvec))
        elif spacecraft == 'tq23':
            return 0.5 * (np.outer(self.tq2_unitvec,self.tq2_unitvec) - \
                np.outer(self.tq3_unitvec,self.tq3_unitvec))
        elif spacecraft == 'tq31':
            return 0.5 * (np.outer(self.tq3_unitvec,self.tq3_unitvec) - \
                np.outer(self.tq1_unitvec,self.tq1_unitvec))
        else:
            raise ValueError("TianQin spacecraft combination should be 12, or 23, or 31")

    def antenna_pattern(self,ra,dec,pol,t_gps,spacecraft):
        theta = pi/2 - dec
        phi = ra
        
        #n is the inverse propagation direction, n and v are orthogonal unit vectors in propagation plane
        n = np.array([cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)])
        u = np.array([-sin(phi),          cos(phi),           0])
        v = np.array([cos(phi)*cos(theta),sin(phi)*cos(theta),-sin(theta)])

        #apply polarization rotation
        upol = cos(2*psi)*u + sin(2*psi)*v
        vpol = -sin(2*psi)*u + cos(2*psi)*v


        du = np.dot(self.response(spacecraft),upol)
        dv = np.dot(self.response(spacecraft),vpol)
        fplus = (upol * du - vpol * dv).sum()
        fcross = (upol * dv + vpol * du).sum()

        return fplus,fcross