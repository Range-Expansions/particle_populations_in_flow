import numpy as np
import weakref


class Simulation_2d(object):

    def __init__(self, Lx=1., Ly=1., interaction_length=0.02, num_particles=300,
                 num_populations = 2, dt = 0.1):

        self.Lx = Lx
        self.Ly = Ly
        self.L = np.array([Lx, Ly])
        self.dt = dt

        self.interaction_length = interaction_length

        self.num_particles = num_particles
        self.num_populations = num_populations

        self.num_bins = np.int32(self.L/interaction_length)

        # Create particles randomly for now
        starting_num = num_particles / num_populations

        self.particle_list = []
        for cur_pop in range(num_populations):
            for i in range(starting_num):
                cur_position = np.random.rand(2) * self.L

                cur_grid = np.int32(cur_position / interaction_length)

                new_particle = Particle(self, cur_pop, cur_position, cur_grid)

                self.particle_list.append(new_particle)

        # Setup the grid
        self.grid = np.zeros((self.num_bins[0], self.num_bins[1], self.num_populations), dtype=np.int32)
        for cur_particle in self.particle_list:
            xy = cur_particle.grid_point
            pop_num = cur_particle.pop_type

            self.grid[xy[0], xy[1], pop_num] += 1

    def react(self):
        """Right now, simple concentration-based growth"""

        for cur_particle in self.particle_list:
            if cur_particle.pop_type != self.num_populations: # Last type is the concentration field
                x = cur_particle.grid_point[0]
                y = cur_particle.grid_point[1]

                num_c = self.grid[x, y, self.num_populations]
                rand = np.random.rand()

                prob = num_c * cur_particle.k * self.dt

                if rand < prob: # react...the issue is that we actually need the ID of each particle...
                    pass



    def run(self, num_iterations):
        for i in range(num_iterations):
            # Move
            for cur_particle in self.particle_list:
                cur_particle.move()


class Particle(object):
    def __init__(self, simulation, pop_type, position, grid_point, D=1.0, k = 1.0):

        self.sim = weakref.proxy(simulation)

        self.pop_type = np.int32(pop_type)
        self.position = np.float32(position)
        self.grid_point = np.int32(grid_point)

        self.D = D
        self.k = k

    def move(self):

        sim = self.sim

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] -= 1

        self.position += np.sqrt(2 * self.D * sim.dt) * np.random.randn(2)

        # Deal with moving out of the system...bounce back
        Lx = sim.L[0]
        Ly = sim.L[1]

        x = self.position[0]
        y = self.position[1]

        if (x < 0):
            x = np.abs(x)
        if (x > Lx):
            x = Lx - (x - Lx)

        if (y < 0):
            y = np.abs(y)
        if (y > Ly):
            y = Ly - (y - Ly)

        self.position[0] = x
        self.position[1] = y

        self.grid_point = np.int32(self.position / sim.interaction_length)

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] += 1