import numpy as np

tolerance = 10.**-9.

class Simulation_2d(object):

    def __init__(self, Lx=1., Ly=1., z = .1,
                 N = 10., R = 4., time_prefactor = 0.1,
                 mu_c = 1.0, mu_list = None,
                 Dc = 1.0, D_list = None,
                 D_nutrient = 1.0):

        self.phys_Lx = Lx
        self.phys_Ly = Ly
        self.phys_z = z

        self.N = N # Number of particles per unit area
        self.R = R # Resolution: Number of interaction lengths a deme of size Lc is divided into
        self.time_prefactor = time_prefactor

        self.phys_muc = mu_c
        self.phys_mu_list = np.array(mu_list, dtype=np.float32)
        self.phys_Dc = Dc
        self.phys_D_list = np.array(D_list, dtype=np.float32)

        self.phys_D_nutrient = D_nutrient

        #### Define Characteristic Length and Time Scales ####
        self.Lc = 2*np.sqrt(self.phys_Dc/self.phys_muc)
        self.Tc = 1./self.phys_muc

        #### Define Dimensionless Parameters ####
        self.dim_Di_list = self.phys_D_list/(4*self.phys_Dc)
        print 'dim_Di:', self.dim_Di_list

        self.dim_D_nutrient = self.phys_D_nutrient/(4*self.phys_Dc)
        print 'dim_D_nutrient:', self.dim_D_nutrient

        self.dim_Gi_list = self.phys_mu_list/self.phys_muc
        print 'dim_Gi:', self.dim_Gi_list

        self.dim_Dgi_list = self.dim_Gi_list/(self.N*self.Lc**2) # Two-dimensional
        print 'dim_Dgi:', self.dim_Dgi_list

        #### Define Simulation Parameters ####
        self.phys_delta = self.Lc/self.R # The interaction length. R should be larger than or equal to 1, always.
        self.dim_delta = self.phys_delta / self.Lc

        self.N_delta = self.N * self.phys_delta**2 # Average number of particles inside the interaction radius
        self.N_L = self.N * self.Lc**2 # Average number of particles inside a deme. Controls stochasticity.

        self.micro_Gi_list = self.dim_Gi_list / self.N_delta # The microscopic reaction rates. Dimensionless.
        print 'Microscopic Gi:', self.micro_Gi_list
        self.dim_dt = time_prefactor * (1./np.max(self.micro_Gi_list))
        self.dim_dx = 1./self.R # The reaction radius spacing in dimensionless units
        print 'Time step (to resolve microscopic reaction rates):', self.dim_dt

        ##### Initialize Particles and Grid #####
        self.dim_Lx = self.phys_Lx / self.Lc
        self.dim_Ly = self.phys_Ly / self.Lc

        self.dim_L_array = np.array(self.dim_Lx, self.dim_Ly)

        self.num_populations = len(self.phys_mu_list)
        self.nutrient_id = self.num_populations # The ID corresponding to the nutrient field

        self.num_fields = self.num_populations + 1 # We need the concentration field

        self.num_bins_x = np.int32(self.phys_Lx/self.phys_delta)
        self.num_bins_y = np.int32(self.phys_Ly/self.phys_delta)

        #### Create the particle dictionary ####
        particle_id_num = 0
        self.particle_dict = {}

        #### Inoculate Nutrients ####
        # Inoculate nutrient particles first. N particles per deme, roughly. The carrying capacity, basically.
        total_num_nutrients = np.int32(self.N * self.phys_Lx * self.phys_Ly)

        print 'Total number of nutrients:', total_num_nutrients

        for i in range(total_num_nutrients):
            # Scatter randomly in space throughout the system. Positions are stored in NON-DIMENSIONAL SPACE
            particle_id = particle_id_num
            cur_position = np.random.rand(2) * self.dim_L_array

            cur_grid = np.int32(cur_position / self.dim_delta)

            new_particle = Particle(self, self.nutrient_id, cur_position, cur_grid,
                                    D=self.dim_D_nutrient, k=0) # k is zero as nutrient decay is correleated with other growth

            self.particle_dict[particle_id] = new_particle

            particle_id_num += 1


        #### Inoculate Populations ####

        # Inoculate them in a circle of width N. We can do this by drawing a random R and theta
        # over a specified range

        # Inoculate a density of N particles all over the circle...
        total_num_population = np.int32(self.N * np.pi*self.phys_z**2)

        print 'Number of particles in initial droplet:', total_num_population

        x0 = self.dim_Lx/2.
        y0 = self.dim_Ly/2.

        for _ in range(total_num_population):

            r = np.random.uniform(0, self.phys_z/self.Lc)
            theta = np.random.uniform(0, 2*np.pi)

            x = x0 + r*np.cos(theta)
            y = y0 + r*np.sin(theta)

            cur_position = np.array([x, y], dtype=np.float32)

            cur_grid = np.int32(cur_position / self.dim_delta)

            pop_type = np.random.randint(self.num_populations)
            new_particle = Particle(self, pop_type, cur_position, cur_grid,
                                    D = self.dim_Di_list[pop_type],
                                    k = self.micro_Gi_list[pop_type])

            self.particle_dict[particle_id_num] = new_particle

            particle_id_num += 1

        #### Setup the grid ####

        self.grid = np.zeros((self.num_bins_x, self.num_bins_y, self.num_fields), dtype=np.int32)
        for cur_particle in self.particle_dict.values():
            xy = cur_particle.grid_point
            pop_num = cur_particle.pop_type

            self.grid[xy[0], xy[1], pop_num] += 1

        self.total_growth_grid = np.zeros((self.num_bins_x, self.num_bins_y), dtype=np.int32)

    def react(self):
        """Right now, simple concentration-based growth"""

        ##### REACT POPULATION PARTICLES ####

        particles_to_add = []
        positions_to_increase = []

        concentration_index = self.num_fields - 1

        for cur_key in self.particle_dict:
            cur_particle = self.particle_dict[cur_key]
            x = cur_particle.grid_point[0]
            y = cur_particle.grid_point[1]

            if cur_particle.pop_type != concentration_index: # Last type is the concentration field

                num_c = self.grid[x, y, concentration_index]

                prob = num_c * cur_particle.k * self.dim_dt
                rand = np.random.rand()

                if rand < prob: # React!
                    new_particle = cur_particle.birth()
                    particles_to_add.append(new_particle)
                    positions_to_increase.append([x, y, new_particle.pop_type])
                    self.total_growth_grid[x, y] += 1

        #### ADJUST THE NUTRIENT FIELD APPROPRIATELY ####

        keys_to_delete = []
        positions_to_decrease = []

        # Adjust the nutrient field appropriately, based on the growth of the others
        for cur_key in self.particle_dict:
            cur_particle = self.particle_dict[cur_key]
            x = cur_particle.grid_point[0]
            y = cur_particle.grid_point[1]

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

class Particle(object):
    def __init__(self, simulation, pop_type, position, grid_point, D=1.0, k = 1.0):

        self.sim = simulation

        self.pop_type = np.int32(pop_type)
        self.position = np.float32(position)
        self.grid_point = np.int32(grid_point)

        self.D = D
        self.k = k

    def move(self):

        sim = self.sim

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] -= 1

        rand2 = np.random.randn(2)

        self.position += np.sqrt(2 * self.D * sim.dim_dt) * rand2

        # Deal with moving out of the system...bounce back
        # The issue with bounceback is that if you move farther that the system size twice,
        # due to the randn draw, you can run into trouble...
        Lx = sim.dim_Lx
        Ly = sim.dim_Ly

        x = self.position[0]
        y = self.position[1]

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

        self.position[0] = x
        self.position[1] = y

        gridx = np.int32(x / sim.dim_delta)
        gridy = np.int32(y / sim.dim_delta)

        self.grid_point[0] = gridx
        self.grid_point[1] = gridy

        xout = (self.grid_point[0] < 0) or (self.grid_point[0] > self.sim.num_bins_x - 1)
        yout = (self.grid_point[1] < 0) or (self.grid_point[1] > self.sim.num_bins_y - 1)

        if xout or yout:
            print 'out of bounds, wtf'
            print 'position:', self.position
            print 'random:', rand2
            print

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] += 1

    def birth(self):
        return Particle(self.sim, self.pop_type, self.position, self.grid_point,
                        D=self.D, k=self.k)