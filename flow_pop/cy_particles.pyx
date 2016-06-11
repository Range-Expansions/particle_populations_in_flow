import numpy as np
cimport numpy as np

cdef float tolerance = 10.**-9.

class Simulation_2d(object):

    def __init__(self, Lx=1., Ly=1., interaction_length=0.02, num_particles=300,
                 num_nutrients = 10**4,
                 num_populations = 2, dt = 0.1):

        self.Lx = Lx
        self.Ly = Ly
        self.L = np.array([Lx, Ly])
        self.dt = dt

        self.interaction_length = interaction_length

        self.num_particles = num_particles
        self.num_nutrients = num_nutrients

        self.num_populations = num_populations
        self.num_fields = num_populations + 1 # We need the concentration field

        self.num_bins = np.int32(self.L/interaction_length)

        # Create particles randomly for now
        starting_num = num_particles / num_populations

        self.particle_dict = {}
        for cur_pop in range(self.num_populations):
            for i in range(starting_num):
                particle_id = cur_pop*starting_num + i
                cur_position = np.random.rand(2) * self.L

                cur_grid = np.int32(cur_position / interaction_length)

                new_particle = Particle(self, cur_pop, cur_position, cur_grid)

                self.particle_dict[particle_id] = new_particle

        # Create nutrients randomly for now as well...
        for i in range(num_nutrients):
            particle_id = self.num_populations * starting_num + i
            cur_position = np.random.rand(2) * self.L

            cur_grid = np.int32(cur_position / interaction_length)

            new_particle = Particle(self, self.num_populations, cur_position, cur_grid)

            self.particle_dict[particle_id] = new_particle

        # Setup the grid
        self.grid = np.zeros((self.num_bins[0], self.num_bins[1], self.num_fields), dtype=np.int32)
        for cur_particle in self.particle_dict.values():
            xy = cur_particle.grid_point
            pop_num = cur_particle.pop_type

            self.grid[xy[0], xy[1], pop_num] += 1

    def react(self):
        """Right now, simple concentration-based growth"""


        particles_to_add = []
        positions_to_increase = []

        keys_to_delete = []
        positions_to_decrease = []

        concentration_index = self.num_fields - 1

        for cur_key in self.particle_dict:
            cur_particle = self.particle_dict[cur_key]
            x = cur_particle.grid_point[0]
            y = cur_particle.grid_point[1]

            if cur_particle.pop_type != concentration_index: # Last type is the concentration field

                num_c = self.grid[x, y, concentration_index]

                prob = num_c * cur_particle.k * self.dt
                rand = np.random.rand()

                if rand < prob: # React!
                    new_particle = cur_particle.birth()
                    particles_to_add.append(new_particle)
                    positions_to_increase.append([x, y, new_particle.pop_type])

            else:
                total_n = 0
                for i in range(self.num_populations): # Loop over all possible 2-pair interactions
                    total_n += self.grid[x, y, i]

                prob = total_n * cur_particle.k * self.dt
                rand = np.random.rand()

                if rand < prob: # Die
                    positions_to_decrease.append([x, y, cur_particle.pop_type])
                    keys_to_delete.append(cur_key)

        # Delete dead particles
        for cur_key in keys_to_delete:
            del self.particle_dict[cur_key]

        # Add new particles to the dictionary
        max_key = np.max(self.particle_dict.keys())
        count = 1
        for cur_particle in particles_to_add:
            self.particle_dict[max_key + count] = cur_particle
            count += 1

        # Update the grid
        for xyc in positions_to_increase:
            self.grid[xyc[0], xyc[1], xyc[2]] += 1
        for xyc in positions_to_decrease:
            self.grid[xyc[0], xyc[1], xyc[2]] -= 1

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

        self.position += np.sqrt(2 * self.D * sim.dt) * rand2

        # Deal with moving out of the system...bounce back
        # The issue with bounceback is that if you move farther that the system size twice,
        # due to the randn draw, you can run into trouble...
        Lx = sim.L[0]
        Ly = sim.L[1]

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

        gridx = np.int32(x / sim.interaction_length)
        gridy = np.int32(y / sim.interaction_length)

        self.grid_point[0] = gridx
        self.grid_point[1] = gridy

        xout = (self.grid_point[0] < 0) or (self.grid_point[0] > self.sim.num_bins[0] - 1)
        yout = (self.grid_point[1] < 0) or (self.grid_point[1] > self.sim.num_bins[1] - 1)

        if xout or yout:
            print 'out of bounds, wtf'
            print 'position:', self.position
            print 'random:', rand2
            print

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] += 1

    def birth(self):
        return Particle(self.sim, self.pop_type, self.position, self.grid_point,
                        D=self.D, k=self.k)