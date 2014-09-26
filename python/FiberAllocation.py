import numpy as np
import matplotlib.pyplot as plt
import struct 
import healpy as hp

class ObjectCatalog(object):
    """
    Data Structure that defines a galaxy catalog.
    """
    def __init__(self,
    random=False, #is it a random catalog?
    filein=False,
    n_points=0):
    
        self.ra = np.array([0.0]) # right ascension
        self.dec = np.array([0.0]) # declination
        self.zz = np.array([0.0]) # redshift
        self.tag = np.array([-1]) # galaxy type
        self.ID = np.array([-1]) # galaxy unique ID
        self.priority = np.array([-1]) # galaxy priority
        self.maxnumobs = np.array([-1]) #maximum number of observations.
        self.healpix_pixels = np.array([-1])
        self.healpix_n_side = -1
    
        if((filein!=False) & (random==False)):
            f = open(filein, "rb")
            sizeint = struct.calcsize("i")
            sizefloat = struct.calcsize("f")
            n_p = f.read(sizeint)
            n_p = (struct.unpack('i', n_p))[0]
            
            format_f = str(n_p)+'f'
            format_i = str(n_p)+'i'
            
            array_data = f.read(n_p*sizefloat)
            self.ra = np.array(struct.unpack(format_f, array_data))
            
            array_data = f.read(n_p*sizefloat)
            self.dec = np.array(struct.unpack(format_f, array_data))
            
            array_data = f.read(n_p*sizefloat)
            self.zz = np.array(struct.unpack(format_f, array_data))
            
            array_data = f.read(n_p*sizeint)
            self.tag = np.array(struct.unpack(format_i, array_data))
            
            array_data = f.read(n_p*sizefloat)
            self.priority = np.array(struct.unpack(format_f, array_data))
            
            array_data = f.read(n_p*sizeint)
            self.maxnumobs = np.array(struct.unpack(format_i, array_data))
            
            self.ID = np.arange(n_p)
            
            print n_p
            f.close()
        if((random==True)&(n_points>0)):
            self.ra = np.random.random(n_points)*360.0
            self.dec = 90.0 - np.arccos((np.random.random(n_points)*2.0) - 1.0)*180.0/np.pi
            self.zz = np.random.random(n_points)
            self.tag = np.int_(np.random.random(n_points)*5)
            self.ID = np.arange(n_points)
            self.priority = np.int_(np.random.random(n_points)*5)
            self.maxnumobs = np.int_(np.random.random(n_points)*3)
        
        #make the initization of new arrays describing the geometry of the problem
        self.theta = (90.0-self.dec)*np.pi/180.0
        self.phi = self.ra*np.pi/180.0
        
        self.nhat0 = np.sin(self.theta)*np.cos(self.phi)
        self.nhat1 = np.sin(self.theta)*np.sin(self.phi)
        self.nhat2 = np.cos(self.theta)
        
    def healpixelize(self, n_side):
        """
        For each galaxy, finds the corresponding pixel number in healpix
        pixelization of n_side.
        """
        self.healpix_n_side = n_side
        self.healpix_pixels = hp.vec2pix(self.healpix_n_side,self.nhat0, self.nhat1, self.nhat2)
        
    def select(self, select_id):
        if(np.size(select_id)>0):
            new_objects = ObjectCatalog(n_points=np.size(select_id), random=True)
            new_objects.ra = self.ra[select_id]
            new_objects.dec = self.dec[select_id]
            new_objects.zz = self.zz[select_id]
            new_objects.tag = self.tag[select_id]
            new_objects.ID = self.ID[select_id]
            new_objects.maxnumobs = self.maxnumobs[select_id]
            new_objects.theta = self.theta[select_id]
            new_objects.phi = self.phi[select_id]
            new_objects.nhat0 = self.nhat0[select_id]
            new_objects.nhat1 = self.nhat1[select_id]
            new_objects.nhat2 = self.nhat2[select_id]
        return new_objects
            
            
        
    def writebinary(self,
                        filein="catalog.dat"):
        f=open(filein, "wb")
        n_p = np.size(self.ra)
        f.write(struct.pack('i', n_p))
        sizeint = struct.calcsize("i")
        sizefloat = struct.calcsize("f")
        
        format_f = str(n_p)+'f'
        format_i = str(n_p)+'i'
            
        bytedata = struct.pack(format_f, *(list(self.ra)))
        f.write(bytedata)
        
        bytedata = struct.pack(format_f, *(list(self.dec)))
        f.write(bytedata)
        
        bytedata = struct.pack(format_f, *(list(self.zz)))
        f.write(bytedata)
        
        bytedata = struct.pack(format_i, *(list(self.tag)))
        f.write(bytedata)
        
        bytedata = struct.pack(format_f, *(list(self.priority)))
        f.write(bytedata)
        
        bytedata = struct.pack(format_i, *(list(self.maxnumobs)))
        f.write(bytedata)

        
        f.close()
    
class FiberSetup(object):
    """
    Data Structure that defines the fiber setup
    """
    def __init__(self,
    random=False, #is it a random catalog?
    filein=False,
    plate_radius=1.65, # in degrees
    patrol_radius_min=0.0, #in mm
    patrol_radius_max=6.0, #in mm
    n_fibers=0):
        if((filein!=False)&(random==False)):
            data = np.loadtxt(filein)
            self.fiber_ID = data[:,0]
            self.positioner_ID = data[:,1]
            self.spectro_ID = data[:,2]
            self.x = data[:,3]
            self.y = data[:,4]
            self.z = data[:,5]
        self.plate_radius = plate_radius*np.pi/180.0
        self.patrol_radius_min = patrol_radius_min
        self.patrol_radius_max = patrol_radius_max
    
class TilingSetup(object):
    def __init__(self,
                 random=False,
                 filein=False,
                 n_tiles=0):
        if((filein!=False)&(random==False)):
            f = open(filein, "rb")
            sizeint = struct.calcsize("i")
            sizefloat = struct.calcsize("f")
            n_p = f.read(sizeint)
            n_p = (struct.unpack('i', n_p))[0]
            format_f = str(n_p)+'f'
            format_i = str(n_p)+'i'
            print n_p

            array_data = f.read(n_p*sizefloat)
            self.ra = np.array(struct.unpack(format_f, array_data))
            
            array_data = f.read(n_p*sizefloat)
            self.dec = np.array(struct.unpack(format_f, array_data))
            
            array_data = f.read(n_p*sizeint)
            self.passnumber = np.array(struct.unpack(format_i, array_data))
            
            f.close()
            
        #make the initization of new arrays describing the geometry of the problem
        self.theta = (90.0-self.dec)*np.pi/180.0
        self.phi = self.ra*np.pi/180.0
        
        self.nhat0 = np.sin(self.theta)*np.cos(self.phi)
        self.nhat1 = np.sin(self.theta)*np.sin(self.phi)
        self.nhat2 = np.cos(self.theta)
        self.healpix_n_side = -1
        self.healpix_pixels = np.array([-1])
        
    def healpixelize(self, n_side):
        """
        For each file, finds the corresponding pixel number in healpix
        pixelization of n_side.
        """
        self.healpix_n_side = n_side
        self.healpix_pixels = hp.vec2pix(self.healpix_n_side,self.nhat0, self.nhat1, self.nhat2)
        
   

def plate_dist(theta):
        """
        Returns the radial distance on the plate (mm) given the angle (radians).
        This is a fit to the provided data
        """
        p = np.array([8.297E5,-1750.0,1.394E4,0.0])
        radius = 0.0
        for i in range(4):
            radius = theta*radius + p[i]
        return radius


def get_rotation_matrix(theta_final=0.0, phi_final=0.0):
    """
    Computes the rotation matrix if the initial position is theta=0.0, phi=0.0
    to the final position over the sphere theta_final, phi_final
    """
    id_matrix = np.zeros((3,3))
    id_matrix[0][0] = 1.0; id_matrix[1][1] = 1.0; id_matrix[2][2]=1.0
    id_matrix = np.matrix(id_matrix)
    if(theta_final>0.0):
        original_center = np.array([0.0,0.0])
        new_center  = np.array([theta_final, phi_final])

        #vector positions
        vec_original = hp.rotator.dir2vec(original_center)
        vec_new  = hp.rotator.dir2vec(new_center)

        # direction and angle of the rotation
        vec_rot_axis = np.cross(vec_original, vec_new)    
        norm_rot_axis = np.linalg.norm(vec_rot_axis)
    
        if(((norm_rot_axis-1.0)>0) & ((norm_rot_axis-1.0)<1E-6)):
            norm_rot_axis = 1.0
    
        rot_angle = np.arccos(np.inner(vec_original, vec_new))
        vec_rot_axis = vec_rot_axis/norm_rot_axis

    
        #Rodrigues' rotation formula
        
        tensor_matrix = np.zeros((3,3))
        tensor_matrix = np.matrix(vec_rot_axis).T*np.matrix(vec_rot_axis)

        cross_matrix = np.zeros((3,3))
        cross_matrix[0][1] = -vec_rot_axis[2]; cross_matrix[0][2] = vec_rot_axis[1]; 
        cross_matrix[1][0] = vec_rot_axis[2]; cross_matrix[1][2] = -vec_rot_axis[0]; 
        cross_matrix[2][0] = -vec_rot_axis[1]; cross_matrix[2][1] = vec_rot_axis[0]; 
        cross_matrix = np.matrix(cross_matrix)

        rot_matrix = np.cos(rot_angle)*id_matrix + np.sin(rot_angle)*cross_matrix + (1.0-np.cos(rot_angle))*tensor_matrix

        rotated_theta, rotated_phi = hp.rotator.rotateDirection(rot_matrix, original_center[0], phi=original_center[1])
        rotated_phi = rotated_phi%(2.0*np.pi)

        delta = (rotated_theta - new_center[0])**2 + (rotated_phi - new_center[1])**2
        vec_final = hp.rotator.dir2vec(np.array([rotated_theta, rotated_phi]))
        if((delta > 1E-6) & (abs(vec_final[2])<0.99)):
            print " ", i_tiling
            print "rotated", rotated_theta, rotated_phi
            print "wanted", new_center
            print "vec wanted", vec_new 
            print "vec original", vec_original
            print "vec done", vec_final
            print "rot angle", rot_angle
            print "vec rotation axis", vec_rot_axis
            print "error", delta
            status = 1
        status = 0
    else:
        rot_matrix = id_matrix
        status = 0
    return status, rot_matrix


def radec2xy(object_ra, object_dec, tile_ra, tile_dec):
    """
    Returns the x,y coordinats of an object on the plate.
    It takes as an input the ra,dec coordinates ob the object 
    and the ra,dec coordinates of the plate's center.
    """
    object_theta = (90.0 - object_dec)*np.pi/180.0
    object_phi = object_ra*np.pi/180.0
    o_hat0 = np.sin(object_theta)*np.cos(object_phi)
    o_hat1 = np.sin(object_theta)*np.sin(object_phi)
    o_hat2 = np.cos(object_theta)
    
    tile_theta = (90.0 - tile_dec)*np.pi/180.0
    tile_phi = tile_ra*np.pi/180.0
    t_hat0 = np.sin(tile_theta)*np.cos(tile_phi)
    t_hat1 = np.sin(tile_theta)*np.sin(tile_phi)
    t_hat2 = np.cos(tile_theta)
    
    
    #we make a rotation on o_hat, so that t_hat ends up aligned with 
    #the unit vector along z. This is composed by a first rotation around z
    #of an angle pi/2 - phi and a second rotation around x by an angle theta, 
    #where theta and phi are the angles describin t_hat.
    
    costheta = t_hat2
    sintheta = np.sqrt(1.0-costheta*costheta) + 1E-10
    cosphi = t_hat0/sintheta
    sinphi = t_hat1/sintheta
    
    #First rotation, taking into account that cos(pi/2 -phi) = sin(phi) and sin(pi/2-phi)=cos(phi)
    n_hat0 = sinphi*o_hat0 - cosphi*o_hat1
    n_hat1 = cosphi*o_hat0 + sinphi*o_hat1
    n_hat2 = o_hat2
    
    #Second rotation
    nn_hat0 = n_hat0
    nn_hat1 = costheta*n_hat1 - sintheta*n_hat2
    nn_hat2 = sintheta*n_hat1 + costheta*n_hat2
    
    #Now find the radius on the plate
    theta = np.sqrt(nn_hat0*nn_hat0 + nn_hat1*nn_hat1)
    radius = plate_dist(theta)
    x = radius * nn_hat0/theta
    y = radius * nn_hat1/theta
    
    return x,y
    
    
def find_available_galaxies(fiber_set, tile_set, object_set, 
                            tile_ID, filename="available_galaxies"):
    
    assert tile_set.healpix_n_side==object_set.healpix_n_side, "healpix n_side is different for tiles and objects"
    assert tile_set.healpix_n_side>0, "healpix n_side is negative"
    
    tile_theta = tile_set.theta[tile_ID]
    tile_phi = tile_set.phi[tile_ID]  
    tile_center = np.array([tile_theta,tile_phi])
    tile_vector = hp.rotator.dir2vec(tile_center)
    pix = hp.query_disc(tile_set.healpix_n_side, tile_vector, fiber_set.plate_radius, inclusive=True) 
    
    # stores the galaxies that fall into the tile 
    objects_in_id = np.empty((0))
    
    for i_pix in pix:
        tmp_id_in = np.where(object_set.healpix_pixels==i_pix)
        tmp_id_in = tmp_id_in[0]
        objects_in_id = np.append(objects_in_id, tmp_id_in)
        
    objects_in_id = np.int_(objects_in_id)
    #if(not(i_tiling%(n_tilings/20))):
        
    selected_objects = object_set.select(objects_in_id)

    #these selected objects must have new a x,y coordinates in the focal plane
    selected_x, selected_y = radec2xy(selected_objects.ra, selected_objects.dec, tile_set.ra[tile_ID], tile_set.dec[tile_ID])
    
    n_fibers = np.size(fiber_set.x)
    fiber_list = np.arange(n_fibers)
    #np.random.shuffle(fiber_list)
    #print fiber_list

    out = open("tile_%d_%s.dat"%(tile_ID, filename), "w")
    
    for fiber_i in fiber_list:
        fiber_x = fiber_set.x[fiber_i]
        fiber_y = fiber_set.y[fiber_i]
        radius = np.sqrt((selected_x - fiber_x)**2 + (selected_y-fiber_y)**2)
        inside = np.where((radius>fiber_set.patrol_radius_min) & (radius<fiber_set.patrol_radius_max))
        inside = inside[0]
        n_available = np.size(inside)
        if(n_available):
            out.write("%d %d %d "%(tile_ID, fiber_i, n_available))
            id_available = ''.join('%d ' % i for i in selected_objects.ID[inside])
            out.write("%s\n"%(id_available))
        else:
            out.write("%d %d %d \n"%(tile_ID, fiber_i, n_available))
    out.close()
    return selected_x, selected_y
