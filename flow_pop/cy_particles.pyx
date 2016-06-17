import numpy as np
cimport numpy as np

from cython_gsl cimport *

cdef float tolerance = 10.**-9.

cdef class Simulation_2d(object):

    cdef:
        unsigned long int seed

        public float phys_Lx
        public float phys_Ly
        public float phys_z
        public float phys_N
        public float droplet_density
        public float R
        public float time_prefactor
        public float phys_muc
        public float[:] phys_mu_list
        public float phys_Dc
        public float[:] phys_D_list
        public float phys_D_nutrient

        public float Lc
        public float Tc
        public float[:] dim_Di_list
        public float dim_D_nutrient
        public float[:] dim_Gi_list
        public float[:] dim_Dgi_list

        public float phys_delta
        public float dim_delta
        public float N_delta
        public float N_L

        public float[:] micro_Gi_list

        public float dim_dt
        public float dim_dx

        public float dim_Lx
        public float dim_Ly

        public float[:] dim_L_array

        public int num_populations
        public int nutrient_id
        public int num_fields

        public int num_bins_x
        public int num_bins_y

        public dict particle_dict
        public int[:, :, :] grid

        public int[:, :] total_growth_grid

        gsl_rng *random_generator

    def __cinit__(self, unsigned long int seed = 0, **kwargs):
        self.seed = seed
        cdef gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937)
        gsl_rng_set(r, self.seed)
        self.random_generator = r

    def __dealloc__(self):
        gsl_rng_free(self.random_generator)

    def __init__(self, float Lx=1., float Ly=1., float z = .1,
                 float N = 10., float R = 4., float time_prefactor = 0.1,
                 float droplet_density=1.0,
                 float mu_c = 1.0, mu_list = None,
                 float Dc = 1.0, D_list = None,
                 float D_nutrient = 1.0, **kwargs):

        self.phys_Lx = Lx
        self.phys_Ly = Ly
        self.phys_z = z

        self.phys_N = N # Number of particles per unit area
        self.droplet_density = droplet_density # Number of particles per unit area in the droplet
        self.R = R # Resolution: Number of interaction lengths a deme of size Lc is divided into
        self.time_prefactor = time_prefactor

        self.phys_muc = mu_c
        self.phys_mu_list = np.array(mu_list, dtype=np.float32)
        self.phys_Dc = Dc
        self.phys_D_list = np.array(D_list, dtype=np.float32)

        self.phys_D_nutrient = D_nutrient

        #### Define Characteristic Length and Time Scales ####
        self.Lc = 2*np.sqrt(self.phys_Dc/self.phys_muc)
        print 'Lc (effective deme size, physical units):', self.Lc
        self.Tc = 1./self.phys_muc
        print 'Tc (characteristic time scale, physical units):', self.Tc

        #### Define Dimensionless Parameters ####
        self.dim_Di_list = D_list/(4*self.phys_Dc)
        print 'dim_Di:', np.asarray(self.dim_Di_list)

        self.dim_D_nutrient = self.phys_D_nutrient/(4*self.phys_Dc)
        print 'dim_D_nutrient:', np.asarray(self.dim_D_nutrient)

        np_Gi_list = mu_list/self.phys_muc
        self.dim_Gi_list = np_Gi_list
        print 'dim_Gi:', np.asarray(self.dim_Gi_list)

        self.dim_Dgi_list = np_Gi_list/(self.phys_N*self.Lc**2) # Two-dimensional
        print 'dim_Dgi:', np.asarray(self.dim_Dgi_list)

        #### Define Simulation Parameters ####
        self.phys_delta = self.Lc/self.R # The interaction length. R should be larger than or equal to 1, always.
        self.dim_delta = self.phys_delta / self.Lc

        self.N_delta = self.phys_N * self.phys_delta**2 # Average number of particles inside the interaction radius
        self.N_L = self.phys_N * self.Lc**2 # Average number of particles inside a deme. Controls stochasticity.

        self.micro_Gi_list = np_Gi_list / self.N_delta # The microscopic reaction rates. Dimensionless.
        print 'Microscopic Gi:', np.asarray(self.micro_Gi_list)
        self.dim_dt = time_prefactor * (1./np.max(self.micro_Gi_list))
        self.dim_dx = 1./self.R # The reaction radius spacing in dimensionless units
        print 'Time step (to resolve microscopic reaction rates):', self.dim_dt

        ##### Initialize Particles and Grid #####
        self.dim_Lx = self.phys_Lx / self.Lc
        self.dim_Ly = self.phys_Ly / self.Lc

        self.num_populations = len(self.phys_mu_list)
        self.nutrient_id = self.num_populations # The ID corresponding to the nutrient field

        self.num_fields = self.num_populations + 1 # We need the concentration field

        self.num_bins_x = np.int32(self.phys_Lx/self.phys_delta) + 1
        self.num_bins_y = np.int32(self.phys_Ly/self.phys_delta) + 1

        #### Create the particle dictionary ####
        particle_id_num = 0
        self.particle_dict = {}

        #### Inoculate Nutrients ####
        # Inoculate nutrient particles first. phys_N particles per deme, roughly. The carrying capacity, basically.
        total_num_nutrients = np.int32(self.phys_N * self.phys_Lx * self.phys_Ly)

        print 'Total number of nutrients:', total_num_nutrients

        for i in range(total_num_nutrients):
            # Scatter randomly in space throughout the system. Positions are stored in NON-DIMENSIONAL SPACE
            particle_id = particle_id_num

            cur_x = gsl_rng_uniform(self.random_generator) * self.dim_Lx
            cur_y = gsl_rng_uniform(self.random_generator) * self.dim_Ly

            cur_gridx = np.int32(cur_x / self.dim_delta)
            cur_gridy = np.int32(cur_y / self.dim_delta)

            new_particle = Particle(self, self.nutrient_id, cur_x, cur_y, cur_gridx, cur_gridy,
                                    D=self.dim_D_nutrient, k=0) # k is zero as nutrient decay is correleated with other growth

            self.particle_dict[particle_id] = new_particle

            particle_id_num += 1


        #### Inoculate Populations ####

        # Inoculate them in a circle of width phys_N. We can do this by drawing a random R and theta
        # over a specified range

        # Inoculate a density of phys_N particles all over the circle...
        total_num_population = np.int32(self.droplet_density * np.pi*self.phys_z**2)

        print 'Number of particles in initial droplet:', total_num_population

        x0 = self.dim_Lx/2.
        y0 = self.dim_Ly/2.

        for _ in range(total_num_population):

            r = gsl_rng_uniform(self.random_generator)*self.phys_z/self.Lc
            theta = gsl_rng_uniform(self.random_generator)*2*np.pi

            x = x0 + r*gsl_sf_cos(theta)
            y = y0 + r*gsl_sf_sin(theta)

            xgrid = int(x/ self.dim_delta)
            ygrid = int(y / self.dim_delta)

            pop_type = gsl_rng_uniform_int(self.random_generator, self.num_populations)
            new_particle = Particle(self, pop_type, x, y, xgrid, ygrid,
                                    D = self.dim_Di_list[pop_type],
                                    k = self.micro_Gi_list[pop_type])

            self.particle_dict[particle_id_num] = new_particle

            particle_id_num += 1

        #### Setup the grid ####

        self.grid = np.zeros((self.num_bins_x, self.num_bins_y, self.num_fields), dtype=np.int32)
        for cur_particle in self.particle_dict.values():
            gridx = cur_particle.gridx
            gridy = cur_particle.gridy

            pop_num = cur_particle.pop_type

            self.grid[gridx, gridy, pop_num] += 1

        self.total_growth_grid = np.zeros((self.num_bins_x, self.num_bins_y), dtype=np.int32)

    def react(self):
        """Right now, simple concentration-based growth"""

        ##### REACT POPULATION PARTICLES ####

        particles_to_add = []
        positions_to_increase = []

        concentration_index = self.num_fields - 1

        cdef int cur_key
        cdef Particle cur_particle, new_particle

        cdef int gridx, gridy
        cdef int num_c
        cdef float prob, rand

        for cur_key in self.particle_dict:
            cur_particle = self.particle_dict[cur_key]
            gridx = cur_particle.gridx
            gridy = cur_particle.gridy

            if cur_particle.pop_type != concentration_index: # Last type is the concentration field

                num_c = self.grid[gridx, gridy, concentration_index]

                prob = num_c * cur_particle.k * self.dim_dt
                rand = gsl_rng_uniform(self.random_generator)

                if rand < prob: # React!
                    new_particle = cur_particle.birth()
                    particles_to_add.append(new_particle)
                    positions_to_increase.append([gridx, gridy, new_particle.pop_type])
                    self.total_growth_grid[gridx, gridy] += 1

        #### ADJUST THE NUTRIENT FIELD APPROPRIATELY ####

        keys_to_delete = []
        positions_to_decrease = []

        # Adjust the nutrient field appropriately, based on the growth of the others
        for cur_key in self.particle_dict:
            cur_particle = self.particle_dict[cur_key]
            x = cur_particle.gridx
            y = cur_particle.gridy

            if cur_particle.pop_type == concentration_index: # Last type is the concentration field

                if self.total_growth_grid[x, y] > 0:
                    positions_to_decrease.append([x, y, cur_particle.pop_type])
                    keys_to_delete.append(cur_key)
                    self.total_growth_grid[x, y] -= 1

        #### UPDATE THE GRIDS AND PARTICLE DICTIONARY ####

        # Remove particles that died (nutrients)
        for cur_key in keys_to_delete:
            del self.particle_dict[cur_key]

        # Add new particles to the dictionary that were born
        max_key = np.max(self.particle_dict.keys())
        count = 1

        for cur_particle in particles_to_add:
            self.particle_dict[max_key + count] = cur_particle
            count += 1

        # Update the grid based on populations
        for xyc in positions_to_increase:
            self.grid[xyc[0], xyc[1], xyc[2]] += 1
        for xyc in positions_to_decrease:
            self.grid[xyc[0], xyc[1], xyc[2]] -= 1

        # Reset the growth grid
        self.total_growth_grid[:, :] = 0

    def move(self):
        for cur_particle in self.particle_dict.values():
            cur_particle.move()

    def run(self, num_iterations):
        for i in range(num_iterations):
            self.move()
            self.react()

cdef class Particle(object):
    cdef:
        public Simulation_2d sim
        public int pop_type
        public float x
        public float y
        public int gridx
        public int gridy
        public float D
        public float k

    def __init__(self, Simulation_2d simulation, int pop_type, float x, float y, int gridx, int gridy,
                 float D=1.0, float k = 1.0):

        self.sim = simulation

        self.pop_type = np.int32(pop_type)

        self.x = x
        self.y = y

        self.gridx = gridx
        self.gridy = gridy

        self.D = D
        self.k = k

    def move(self):

        cdef Simulation_2d sim = self.sim

        sim.grid[self.gridx, self.gridy, self.pop_type] -= 1

        cdef float x = self.x
        cdef float y = self.y

        cdef float Lx, Ly, dx, dy
        cdef int gridx, gridy

        cdef float D = self.D
        cdef float dt = sim.dim_dt

        cdef gsl_rng *r = sim.random_generator

        with nogil:

            x += sqrt(2 * D * dt) * gsl_ran_gaussian(r, 1)
            y += sqrt(2 * D * dt) * gsl_ran_gaussian(r, 1)

            # Deal with moving out of the system...bounce back
            # The issue with bounceback is that if you move farther that the system size twice,
            # due to the randn draw, you can run into trouble...
            Lx = sim.dim_Lx
            Ly = sim.dim_Ly

            if (x < 0):
                dx = (-x) % Lx
                x = dx + tolerance
            elif (x > Lx):
                dx = (x - Lx) % Lx # Just to avoid super bounces
                x = Lx - dx - tolerance
            if (y < 0):
                dy = (-y) % Ly
                y = dy + tolerance
            elif (y > Ly):
                dy = (y - Ly) % Ly  # Just to avoid super bounces
                y = Ly - dy - tolerance

            gridx = int(x / sim.dim_delta)
            gridy = int(y / sim.dim_delta)

        self.x = x
        self.y = y

        self.gridx = gridx
        self.gridy = gridy

        xout = (self.gridx < 0) or (self.gridx > self.sim.num_bins_x - 1)
        yout = (self.gridy < 0) or (self.gridy > self.sim.num_bins_y - 1)

        if xout or yout:
            print 'out of bounds, wtf'
            print 'position:', self.x, self.y
            print

        sim.grid[self.gridx, self.gridy, self.pop_type] += 1

    cdef Particle birth(Particle self):
        return Particle(self.sim, self.pop_type, self.x, self.y, self.gridx, self.gridy,
                        D=self.D, k=self.k)